#' Quick GeMLR Fit Without Cross-Validation
#'
#' @param x     Numeric matrix/data.frame of predictors (n x p), ALREADY standardized (Xs from read_data()).
#' @param y     Binary response vector (0/1), length n.
#' @param Indi  Optional indicator (dummy) matrix/data.frame (n x q), appended to logistic design.
#' @param K     Number of clusters (default 2).
#' @param num_gmm If NULL, use all columns for GMM; if 0, must provide vargmm;
#'                if positive integer, pick top-variance num_gmm columns (default 5).
#' @param vargmm Column names or indices to use for GMM (takes precedence over num_gmm).
#' @param nseeds Number of kmeans seeds to try (default 3).
#' @param rangeSeed Upper bound for random seed sampling (default 500).
#' @param vlasso Optional glmnet lambda; if NULL, estimated by cv.glmnet.
#' @param alphaLasso Elastic-net alpha for glmnet (default 0.5).
#' @param verbose Verbosity flag passed to GeMLR internals (0 = silent).
#'
#' @return An object of class "GeMLR_fit".
#' @export
fit_model <- function(x, y, Indi = NULL,
                      K = 2,
                      num_gmm = 5,
                      vargmm = NULL,
                      nseeds = 3, rangeSeed = 500,
                      vlasso = NULL,
                      alphaLasso = 0.5,
                      verbose = 0) {

  # -- inputs (assume x is already standardized: Xs)
  Xs <- as.matrix(x)
  Y  <- as.numeric(y)
  stopifnot(all(Y %in% c(0, 1)))
  if (!is.null(Indi)) Indi <- as.matrix(Indi)

  # -- logistic design matrix
  Xlogit <- if (is.null(Indi)) Xs else cbind(Xs, Indi)
  storage.mode(Xlogit) <- "double"
  colnames(Xlogit) <- make.names(colnames(Xlogit), unique = TRUE)

  # ---------- GMM feature selection on Xs ----------
  p <- ncol(Xs)
  if (!is.null(vargmm) && length(vargmm) > 0L) {
    # use vargmm (names or indices)
    if (is.character(vargmm)) {
      if (is.null(colnames(Xs))) stop("vargmm uses names but Xs has no column names.")
      if (!all(vargmm %in% colnames(Xs))) stop("Some vargmm names not found in Xs.")
      idx_gmm <- match(vargmm, colnames(Xs))
    } else {
      idx_gmm <- as.integer(vargmm)
    }
  } else {
    # derive by num_gmm
    if (is.null(num_gmm)) {
      idx_gmm <- seq_len(p)  # use all
    } else if (identical(num_gmm, 0L) || identical(num_gmm, 0)) {
      stop("num_gmm = 0 requires providing vargmm (names or indices).")
    } else if (is.numeric(num_gmm) && num_gmm > 0 && floor(num_gmm) == num_gmm) {
      X_var <- apply(Xs, 2, var)
      ord   <- order(X_var, decreasing = TRUE)
      k     <- min(num_gmm, p)
      idx_gmm <- ord[seq_len(k)]
    } else {
      stop("num_gmm must be NULL, 0, or a positive integer.")
    }
  }
  if (length(idx_gmm) == 0L || !all(idx_gmm %in% seq_len(p))) {
    stop("Resolved vargmm is empty or out of range.")
  }
  idx_gmm <- sort(unique(idx_gmm))
  X_gmm   <- Xs[, idx_gmm, drop = FALSE]
  dimgmm  <- ncol(X_gmm)

  # -- sanity
  if (any(!is.finite(X_gmm)) || any(!is.finite(Xlogit)) || any(!is.finite(Y))) {
    stop("X/Indi/Y contains NA/Inf.")
  }

  # -- lambda for glmnet (if not provided)
  if (is.null(vlasso)) {
    if (!requireNamespace("glmnet", quietly = TRUE))
      stop("Package 'glmnet' is required.")
    suppressWarnings({
      cv_fit <- glmnet::cv.glmnet(data.matrix(Xlogit), Y, alpha = 1, family = "binomial", nfolds = 5)
    })
    vlasso <- cv_fit$lambda.min
  }

  # -- GeMLR options
  MLMoption <- init_MLMoption(
    alphaLasso = alphaLasso,
    vlasso     = 1,      # placeholder (we set lambdaLasso below)
    numcmp     = K,
    verbose    = verbose,
    DISTR      = "binomial",
    AUC        = 1,
    NOEM       = 0
  )
  MLMoption$lambdaLasso <- rep(vlasso, K)

  # -- seeds for initialization
  seedlist <- sample.int(max(10, rangeSeed), nseeds, replace = FALSE)

  # -- try EM/GMM estimation
  est <- NULL; ok <- FALSE
  for (seed_i in seedlist) {
    MLMoption$kmseed <- seed_i
    one <- NULL
    res <- try({
      suppressWarnings(suppressMessages(
        invisible(capture.output({
          one <- estimateBestSD(X_gmm, Xlogit, Y, MLMoption, seedlist = c(seed_i))
        }, type = "output"))
      ))
    }, silent = TRUE)
    if (!inherits(res, "try-error") && !is.null(one) && !is.null(one$c) && !is.null(one$c$supp)) {
      est <- one; ok <- TRUE; break
    }
  }

  # -- fallback: single LR if EM/GMM fails
  if (!ok) {
    if (!requireNamespace("glmnet", quietly = TRUE))
      stop("Estimation failed and 'glmnet' missing for fallback.")
    suppressWarnings({
      cv_lr <- glmnet::cv.glmnet(data.matrix(Xlogit), Y, family = "binomial", alpha = 1, nfolds = 5)
    })
    p_lr <- as.numeric(stats::predict(cv_lr, newx = data.matrix(Xlogit), s = "lambda.min", type = "response"))
    eps <- 1e-12
    p_safe <- p_lr
    p_safe[!is.finite(p_safe)] <- NA
    p_safe <- pmin(pmax(p_safe, eps), 1 - eps)
    metrics <- list(
      accuracy = mean((p_safe >= 0.5) == Y, na.rm = TRUE),
      auc      = if (requireNamespace("pROC", quietly = TRUE)) as.numeric(pROC::auc(Y, p_safe)) else NA_real_,
      logloss  = -mean(Y * log(p_safe) + (1 - Y) * log(1 - p_safe), na.rm = TRUE)
    )
    lam_idx <- which(cv_lr$glmnet.fit$lambda == cv_lr$lambda.min)
    b_intercept <- cv_lr$glmnet.fit$a0[lam_idx]
    b_coef <- cv_lr$glmnet.fit$beta[, lam_idx, drop = FALSE]
    beta <- rbind("(Intercept)" = b_intercept, as.matrix(b_coef))
    colnames(beta) <- "Cluster 1"
    a2 <- 1
    mu2 <- matrix(0, nrow = dimgmm, ncol = 1)
    sigma2 <- array(diag(dimgmm), dim = c(dimgmm, dimgmm, 1))
    pij <- matrix(1, nrow = nrow(Xs), ncol = 1)
    return(structure(list(
      beta = beta, a2 = a2, mu2 = mu2, sigma2 = sigma2,
      pij = pij, p = p_lr, clusterid = rep(1, nrow(Xs)),
      scaler = NULL, K = 1, vargmm = idx_gmm, metrics = metrics, call = match.call()
    ), class = "GeMLR_fit"))
  }

  # -- normal path: convert and classify
  gmm <- GMMFormatConvert(dimgmm, est$c)
  a2 <- gmm$a; mu2 <- gmm$mu; sigma2 <- gmm$sigma
  beta <- est$beta
  rn <- c("Intercept", colnames(Xlogit))
  rownames(beta) <- rn[seq_len(nrow(beta))]
  colnames(beta) <- paste0("Cluster ", seq_len(ncol(beta)))

  cls <- MLMclassify(a2, mu2, sigma2, beta, X_gmm, Xlogit)
  pyi <- as.numeric(cls$pyi); pij <- cls$pij
  clusterid <- max.col(pij, ties.method = "first")

  # -- safe metrics
  eps <- 1e-12
  p_safe <- pyi
  p_safe[!is.finite(p_safe)] <- NA
  p_safe <- pmin(pmax(p_safe, eps), 1 - eps)
  metrics <- list(
    accuracy = mean((p_safe >= 0.5) == Y, na.rm = TRUE),
    auc      = if (requireNamespace("pROC", quietly = TRUE)) as.numeric(pROC::auc(Y, p_safe)) else NA_real_,
    logloss  = -mean(Y * log(p_safe) + (1 - p_safe), na.rm = TRUE)
  )

  structure(list(
    beta = beta, a2 = a2, mu2 = mu2, sigma2 = sigma2,
    pij = pij, p = pyi, clusterid = clusterid,
    scaler = NULL, K = K, vargmm = idx_gmm, metrics = metrics, call = match.call()
  ), class = "GeMLR_fit")
}

#' Predict Class Labels (expects newx standardized same as training Xs)
#' @export
predict_class <- function(fit, newx, Indi = NULL, threshold = 0.5) {
  stopifnot(inherits(fit, "GeMLR_fit"))
  Xs_new <- as.matrix(newx)  # must already be standardized, same as training Xs
  Xlogit <- if (is.null(Indi)) Xs_new else cbind(Xs_new, as.matrix(Indi))
  colnames(Xlogit) <- make.names(colnames(Xlogit), unique = TRUE)
  X_gmm <- Xs_new[, fit$vargmm, drop = FALSE]

  if (length(fit$a2) == 1 && ncol(fit$beta) == 1) {
    b <- fit$beta[, 1, drop = TRUE]
    eta <- as.numeric(b[1] + Xlogit %*% b[-1])
    p <- 1 / (1 + exp(-eta))
    return(as.integer(p >= threshold))
  }

  cls <- MLMclassify(fit$a2, fit$mu2, fit$sigma2, fit$beta, X_gmm, Xlogit)
  as.integer(as.numeric(cls$pyi) >= threshold)
}





