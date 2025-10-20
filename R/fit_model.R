#' Quick GeMLR Fit Without Cross-Validation
#'
#' @description
#' Fit a generative mixture of logistic regression (GeMLR) quickly (no CV),
#' mainly for a user-friendly preview before cross-validation.
#' The model builds a GMM over continuous features for soft clustering, and
#' fits one logistic regression per cluster. Final class probability is a
#' responsibility-weighted mixture over clusters.
#'
#' @param x Numeric matrix/data.frame of continuous predictors (n x p).
#' @param y Binary response vector (0/1) of length n.
#' @param Indi Optional indicator (dummy) matrix (n x q). Appended to the
#' logistic-regression design matrix (no standardization).
#' @param K Integer, number of clusters (default 2).
#' @param num_gmm Integer, number of continuous columns (by descending variance)
#' to use in the GMM. If \code{vargmm} is provided, it takes precedence.
#' If both are missing, all continuous columns of \code{x} are used.
#' @param vargmm Integer vector of column indices in \code{x} to use for GMM
#' (overrides \code{num_gmm}).
#' @param standardize Logical; whether to standardize \code{x} (default TRUE).
#' @param nseeds Number of k-means seeds to try in EM initialization (default 3).
#' @param rangeSeed Upper bound of the random seed sampling range (default 500).
#' @param vlasso Optional \code{glmnet} lambda value; if NULL, it is estimated
#' by \code{glmnet::cv.glmnet}.
#' @param alphaLasso Elastic-net \code{alpha} for \code{glmnet} (default 0.5).
#' @param verbose Integer verbosity flag (0 = silent).
#'
#' @return An object of class \code{"GeMLR_fit"}:
#' \itemize{
#'   \item \code{beta}: coefficient matrix (intercept + features) per cluster
#'   \item \code{a2, mu2, sigma2}: GMM weights, means, covariances
#'   \item \code{p}: predicted \(P(Y=1|X)\) on the training data
#'   \item \code{pij}: responsibilities (n x K)
#'   \item \code{clusterid}: hard assignments via \code{max.col(pij)}
#'   \item \code{metrics}: list with \code{accuracy}, \code{auc} (if \pkg{pROC} present), \code{logloss}
#'   \item \code{scaler}: list with \code{mu}, \code{sd} used for standardization
#'   \item \code{K}, \code{vargmm}, \code{call}
#' }
#'
#' @details
#' This function relies on internal GeMLR routines: \code{init_MLMoption()},
#' \code{estimateBestSD()}, \code{GMMFormatConvert()}, and \code{MLMclassify()}.
#' If the EM/GMM path fails (e.g., due to numerical issues), the function
#' automatically degrades to a single logistic regression (\code{K=1}) fitted
#' via \code{glmnet}, returning a structure compatible with \code{GeMLR_fit}.
#'
#' @export
#' @importFrom stats predict
fit_model <- function(x, y, Indi = NULL,
                      K = 2,
                      num_gmm = 0,
                      vargmm = NULL,
                      standardize = TRUE,
                      nseeds = 3, rangeSeed = 500,
                      vlasso = NULL,
                      alphaLasso = 0.5,
                      verbose = 0) {
  X <- as.matrix(x); Y <- as.numeric(y)
  stopifnot(all(Y %in% c(0,1)))
  if (!is.null(Indi)) Indi <- as.matrix(Indi)
  
  # Standardize continuous features
  scaler <- NULL
  if (standardize) {
    mu <- colMeans(X); sdv <- apply(X, 2, sd); sdv[sdv == 0] <- 1
    scaler <- list(mu=mu, sd=sdv)
    Xs <- scale(X, center=mu, scale=sdv)
  } else {
    Xs <- X
  }
  
  # Logistic design matrix
  Xlogit <- if (is.null(Indi)) Xs else cbind(Xs, Indi)
  storage.mode(Xlogit) <- "double"
  colnames(Xlogit) <- make.names(colnames(Xlogit), unique=TRUE)
  
  # Select GMM features
  if (!is.null(vargmm)) {
    idx_gmm <- as.integer(vargmm)
  } else if (is.numeric(num_gmm) && num_gmm > 0) {
    vars <- apply(Xs, 2, var)
    ord <- order(vars, decreasing=TRUE)
    idx_gmm <- ord[seq_len(min(num_gmm, ncol(Xs)))]
  } else {
    idx_gmm <- seq_len(ncol(Xs))
  }
  if (length(idx_gmm) == 0L || !all(idx_gmm %in% seq_len(ncol(Xs))))
    stop("vargmm is empty or out of range; please check.")
  X_gmm <- Xs[, idx_gmm, drop=FALSE]
  dimgmm <- ncol(X_gmm)
  
  # Basic checks
  if (any(!is.finite(X_gmm)) || any(!is.finite(Xlogit)) || any(!is.finite(Y)))
    stop("X/Indi/Y contains NA/Inf; please clean the data first.")
  
  # Lasso lambda (if not provided)
  if (is.null(vlasso)) {
    if (!requireNamespace("glmnet", quietly=TRUE))
      stop("Package 'glmnet' is required. Please install.packages('glmnet').")
    suppressWarnings({
      cv_fit <- glmnet::cv.glmnet(data.matrix(Xlogit), Y, alpha=1, family="binomial", nfolds=5)
    })
    vlasso <- cv_fit$lambda.min
  }
  
  # Build MLMoption
  MLMoption <- init_MLMoption(
    alphaLasso = alphaLasso,
    vlasso     = 1,       # placeholder; set lambdaLasso below
    numcmp     = K,
    verbose    = verbose,
    DISTR      = 'binomial',
    AUC        = 1,
    NOEM       = 0
  )
  MLMoption$lambdaLasso <- rep(vlasso, K)
  
  # Seeds to try for k-means init
  seedlist <- sample.int(max(10, rangeSeed), nseeds, replace=FALSE)
  
  # Try estimation; if it fails, degrade to K=1 LR
  est <- NULL
  ok <- FALSE
  for (seed_i in seedlist) {
    MLMoption$kmseed <- seed_i
    one <- NULL
    res <- try({
      suppressWarnings(suppressMessages(
        invisible(capture.output({
          one <- estimateBestSD(X_gmm, Xlogit, Y, MLMoption, seedlist = seed_i)
        }, type="output"))
      ))
    }, silent=TRUE)
    if (!inherits(res, "try-error") && !is.null(one) && !is.null(one$c) && !is.null(one$c$supp)) {
      est <- one; ok <- TRUE; break
    }
  }
  
  # Degrade to single LR if EM/GMM fails
  if (!ok) {
    if (!requireNamespace("glmnet", quietly=TRUE))
      stop("Estimation failed and fallback requires 'glmnet'. Please install it.")
    suppressWarnings({
      cv_lr <- glmnet::cv.glmnet(data.matrix(Xlogit), Y, family="binomial", alpha=1, nfolds=5)
    })
    p_lr <- as.numeric(stats::predict(cv_lr, newx = data.matrix(Xlogit),
                                      s = "lambda.min", type = "response"))
    # Safe logloss
    eps <- 1e-12
    p_safe <- p_lr
    p_safe[!is.finite(p_safe)] <- NA
    p_safe <- pmin(pmax(p_safe, eps), 1 - eps)
    metrics <- list(
      accuracy = mean((p_safe>=0.5)==Y, na.rm = TRUE),
      auc      = if (requireNamespace("pROC", quietly=TRUE)) as.numeric(pROC::auc(Y, p_safe)) else NA_real_,
      logloss  = -mean(Y*log(p_safe) + (1-Y)*log(1-p_safe), na.rm = TRUE)
    )
    # Coefficients at lambda.min
    lam_idx <- which(cv_lr$glmnet.fit$lambda == cv_lr$lambda.min)
    b_intercept <- cv_lr$glmnet.fit$a0[lam_idx]
    b_coef <- cv_lr$glmnet.fit$beta[, lam_idx, drop = FALSE]
    beta <- rbind("(Intercept)" = b_intercept, as.matrix(b_coef))
    colnames(beta) <- "Cluster 1"
    a2 <- 1
    mu2 <- matrix(0, nrow=dimgmm, ncol=1)
    sigma2 <- array(diag(dimgmm), dim=c(dimgmm,dimgmm,1))
    pij <- matrix(1, nrow=nrow(X), ncol=1)
    return(structure(list(
      beta=beta, a2=a2, mu2=mu2, sigma2=sigma2,
      pij=pij, p=p_lr, clusterid=rep(1, nrow(X)),
      scaler=scaler, K=1, vargmm=idx_gmm, metrics=metrics, call=match.call()
    ), class="GeMLR_fit"))
  }
  
  # Normal path: convert + classify
  gmm <- GMMFormatConvert(dimgmm, est$c)
  a2 <- gmm$a; mu2 <- gmm$mu; sigma2 <- gmm$sigma
  beta <- est$beta
  rn <- c("Intercept", colnames(Xlogit))
  rownames(beta) <- rn[seq_len(nrow(beta))]
  colnames(beta) <- paste0("Cluster ", seq_len(ncol(beta)))
  
  cls <- MLMclassify(a2, mu2, sigma2, beta, X_gmm, Xlogit)
  pyi <- as.numeric(cls$pyi); pij <- cls$pij
  clusterid <- max.col(pij, ties.method="first")
  
  # Safe metrics
  eps <- 1e-12
  p_safe <- pyi
  p_safe[!is.finite(p_safe)] <- NA
  p_safe <- pmin(pmax(p_safe, eps), 1 - eps)
  metrics <- list(
    accuracy = mean((p_safe>=0.5)==Y, na.rm = TRUE),
    auc      = if (requireNamespace("pROC", quietly=TRUE)) as.numeric(pROC::auc(Y, p_safe)) else NA_real_,
    logloss  = -mean(Y*log(p_safe) + (1-Y)*log(1-p_safe), na.rm = TRUE)
  )
  
  structure(list(
    beta=beta, a2=a2, mu2=mu2, sigma2=sigma2,
    pij=pij, p=pyi, clusterid=clusterid,
    scaler=scaler, K=K, vargmm=idx_gmm, metrics=metrics, call=match.call()
  ), class="GeMLR_fit")
}

#' Predict Class Labels from a GeMLR Fit
#'
#' @description
#' Predict 0/1 labels for new data given a \code{GeMLR_fit} object. The function
#' automatically standardizes \code{newx} using the scaler stored in \code{fit}.
#' If the model degraded to a single logistic regression (\code{K=1}), this
#' routine returns the thresholded logistic probabilities.
#'
#' @param fit An object returned by \code{\link{fit_model}}.
#' @param newx Numeric matrix/data.frame of new continuous features (n x p).
#' @param Indi Optional indicator matrix for new data (n x q).
#' @param threshold Probability cutoff for class 1 (default 0.5).
#'
#' @return Integer vector of predicted labels (0/1).
#' @export
predict_class <- function(fit, newx, Indi = NULL, threshold = 0.5) {
  stopifnot(inherits(fit, "GeMLR_fit"))
  X <- as.matrix(newx)
  if (!is.null(fit$scaler)) {
    mu <- fit$scaler$mu; sdv <- fit$scaler$sd; sdv[sdv==0] <- 1
    X <- scale(X, center=mu, scale=sdv)
  }
  Xlogit <- if (is.null(Indi)) X else cbind(X, as.matrix(Indi))
  colnames(Xlogit) <- make.names(colnames(Xlogit), unique=TRUE)
  X_gmm <- X[, fit$vargmm, drop=FALSE]
  
  # K=1 fallback path: plain logistic regression
  if (length(fit$a2)==1 && ncol(fit$beta)==1) {
    b <- fit$beta[,1,drop=TRUE]
    eta <- as.numeric(b[1] + Xlogit %*% b[-1])
    p <- 1/(1+exp(-eta))
    return(as.integer(p >= threshold))
  }
  
  cls <- MLMclassify(fit$a2, fit$mu2, fit$sigma2, fit$beta, X_gmm, Xlogit)
  as.integer(as.numeric(cls$pyi) >= threshold)
}
