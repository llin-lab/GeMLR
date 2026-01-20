#' Quick GeMLR Fit Without Cross-Validation
#'
#' @param X Raw (unstandardized) predictor matrix (n x p).
#' @param Xs Standardized predictor matrix (n x p).
#' @param Y Binary response vector (0/1).
#' @param Indi Optional indicator matrix/data.frame.
#' @param K Number of clusters (default 2).
#' @param vargmm Column names or indices defining the POOL for GMM. NULL = all features.
#' @param VS Variable Selection for GMM. NA = use all in pool, Integer = select top VS by variance.
#' @param varreg Column names or indices for LR. NULL = all features.
#' @param vlasso Lambda for Lasso. NULL = auto-compute.
#' @param nseeds Number of kmeans seeds (default 3).
#' @param rangeSeed Upper bound for seed sampling (default 500).
#' @param alphaLasso Elastic-net alpha (default 0.5).
#' @param verbose Verbosity (default 0).
#' @export
fit_model <- function(X, Xs, Y, Indi = NULL,
                      K = 2,
                      vargmm = NULL,
                      VS = NA,
                      varreg = NULL,
                      vlasso = NULL,
                      nseeds = 3, 
                      rangeSeed = 500,
                      alphaLasso = 0.5,
                      verbose = 0) {
  
  X_raw <- as.matrix(X)
  Xs <- as.matrix(Xs)
  Y  <- as.numeric(Y)
  stopifnot(all(Y %in% c(0, 1)))
  if (!is.null(Indi)) Indi <- as.matrix(Indi)
  
  p <- ncol(Xs)
  
  # STEP 1: Define vargmm pool
  if (is.null(vargmm)) {
    vargmm_pool_idx <- seq_len(p)
  } else if (is.character(vargmm)) {
    if (is.null(colnames(Xs))) {
      stop("vargmm uses column names but Xs has no column names.")
    }
    if (!all(vargmm %in% colnames(Xs))) {
      stop("Some vargmm names not found in Xs: ", 
           paste(setdiff(vargmm, colnames(Xs)), collapse = ", "))
    }
    vargmm_pool_idx <- match(vargmm, colnames(Xs))
  } else {
    vargmm_pool_idx <- as.integer(vargmm)
    if (any(vargmm_pool_idx < 1 | vargmm_pool_idx > p)) {
      stop("Some vargmm indices are out of range [1,", p, "]")
    }
  }
  
  # STEP 2: Apply Variable Selection (VS)
  if (is.na(VS)) {
    idx_gmm <- vargmm_pool_idx
  } else {
    VS <- as.integer(VS)
    if (VS < 1) {
      stop("VS must be a positive integer or NA")
    }
    
    pool_data <- Xs[, vargmm_pool_idx, drop = FALSE]
    pool_vars <- apply(pool_data, 2, var, na.rm = TRUE)
    
    n_select <- min(VS, length(vargmm_pool_idx))
    top_positions <- order(pool_vars, decreasing = TRUE)[1:n_select]
    idx_gmm <- vargmm_pool_idx[top_positions]
  }
  
  idx_gmm <- sort(unique(idx_gmm))
  X_gmm   <- Xs[, idx_gmm, drop = FALSE]
  dimgmm  <- ncol(X_gmm)
  # ========================================
  # AUTO-COMPUTE vlasso if not provided
  # CRITICAL: Use raw X for consistency with old version
  # ========================================
  if (is.null(vlasso)) {
    if (verbose > 0) cat("vlasso not specified. Computing via cross-validation using raw X...\n")
    if (!requireNamespace("glmnet", quietly = TRUE)) {
      stop("glmnet package required for automatic vlasso computation. Please install it or provide vlasso manually.")
    }
    # Use RAW X for vlasso computation (consistent with compute_gemlr_params and finalModel)
    set.seed(1)
    Xlogit_temp <- if (is.null(Indi)) X_raw else cbind(X_raw, Indi)
    cv_fit <- suppressWarnings(
      glmnet::cv.glmnet(data.matrix(Xlogit_temp), Y, alpha = 1, family = "binomial", nfolds = 5)
    )
    vlasso <- cv_fit$lambda.min
    if (verbose > 0) cat("Computed vlasso =", vlasso, "\n")
  }
  
  # ========================================
  # Handle varreg (variables for LR)
  # ========================================
  if (is.null(varreg)) {
    # Default: use all variables
    X_reg <- Xs
    idx_reg <- NULL
  } else {
    # User-specified
    if (is.character(varreg)) {
      if (is.null(colnames(Xs))) stop("varreg uses names but Xs has no column names.")
      if (!all(varreg %in% colnames(Xs))) stop("Some varreg names not found in Xs.")
      idx_reg <- match(varreg, colnames(Xs))
    } else {
      idx_reg <- as.integer(varreg)
    }
    
    if (length(idx_reg) == 0L || !all(idx_reg %in% seq_len(p))) {
      stop("Resolved varreg is empty or out of range.")
    }
    X_reg <- Xs[, idx_reg, drop = FALSE]
  }
  
  # -- Build logistic design matrix using standardized features
  Xlogit <- if (is.null(Indi)) X_reg else cbind(X_reg, Indi)
  storage.mode(Xlogit) <- "double"
  colnames(Xlogit) <- make.names(colnames(Xlogit), unique = TRUE)
  
  # -- Sanity check
  if (any(!is.finite(X_gmm)) || any(!is.finite(Xlogit)) || any(!is.finite(Y))) {
    stop("X/Indi/Y contains NA/Inf.")
  }
  
  # -- Initialize GeMLR options
  MLMoption <- init_MLMoption(
    alphaLasso = alphaLasso,
    vlasso     = 1,
    numcmp     = K,
    verbose    = verbose,
    DISTR      = "binomial",
    AUC        = 1,
    NOEM       = 0
  )
  MLMoption$lambdaLasso <- rep(vlasso, K)
  
  set.seed(9)
  # -- Generate seeds for initialization
  seedlist <- sample.int(max(10, rangeSeed), nseeds, replace = FALSE)
  
  # -- Attempt EM/GMM estimation
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
  
  # -- Fallback: single LR if EM/GMM fails
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
    beta <- cbind(beta, beta[, 1, drop = FALSE])
    colnames(beta)[2] <- "LR"
    a2 <- 1
    mu2 <- matrix(0, nrow = dimgmm, ncol = 1)
    sigma2 <- array(diag(dimgmm), dim = c(dimgmm, dimgmm, 1))
    pij <- matrix(1, nrow = nrow(Xs), ncol = 1)
    return(structure(list(
      beta = beta, a2 = a2, mu2 = mu2, sigma2 = sigma2,
      pij = pij, p = p_lr, clusterid = rep(1, nrow(Xs)),
      scaler = NULL, K = 1, 
      vargmm = idx_gmm, varreg = idx_reg, vlasso = vlasso,
      metrics = metrics, call = match.call()
    ), class = "GeMLR_fit"))
  }
  
  # -- Normal path: convert GMM results and classify
  gmm <- GMMFormatConvert(dimgmm, est$c)
  a2 <- gmm$a; mu2 <- gmm$mu; sigma2 <- gmm$sigma
  beta <- est$beta
  rn <- c("Intercept", colnames(Xlogit))
  rownames(beta) <- rn[seq_len(nrow(beta))]
  colnames(beta) <- paste0("Cluster ", seq_len(ncol(beta)))
  
  cls <- MLMclassify(a2, mu2, sigma2, beta, X_gmm, Xlogit)
  pyi <- as.numeric(cls$pyi); pij <- cls$pij
  clusterid <- max.col(pij, ties.method = "first")
  
  # -- Add baseline LR model
  if (requireNamespace("glmnet", quietly = TRUE)) {
    tryCatch({
      set.seed(1)
      suppressWarnings({
        cv_lr <- glmnet::cv.glmnet(
          x = data.matrix(Xlogit), 
          y = Y, 
          alpha = 1, 
          family = "binomial", 
          nfolds = 5
        )
      })
      B_lr <- coef(cv_lr, s = "lambda.min")
      if (nrow(beta) == nrow(B_lr)) {
        beta <- cbind(beta, as.matrix(B_lr))
        colnames(beta)[ncol(beta)] <- "LR"
      }
    }, error = function(e) {
      warning("Failed to add LR baseline: ", e$message)
    })
  }
  
  # -- Compute safe metrics
  eps <- 1e-12
  p_safe <- pyi
  p_safe[!is.finite(p_safe)] <- NA
  p_safe <- pmin(pmax(p_safe, eps), 1 - eps)
  metrics <- list(
    accuracy = mean((p_safe >= 0.5) == Y, na.rm = TRUE),
    auc      = if (requireNamespace("pROC", quietly = TRUE)) as.numeric(pROC::auc(Y, p_safe)) else NA_real_,
    logloss  = -mean(Y * log(p_safe) + (1 - Y) * log(1 - p_safe), na.rm = TRUE)
  )
  
  structure(list(
    beta = beta, a2 = a2, mu2 = mu2, sigma2 = sigma2,
    pij = pij, p = p_safe, clusterid = clusterid,
    scaler = NULL, K = K, 
    vargmm = idx_gmm, varreg = idx_reg, vlasso = vlasso,
    metrics = metrics, call = match.call()
  ), class = "GeMLR_fit")
}

#' Predict Class Labels
#' 
#' @param fit A GeMLR_fit object from fit_model()
#' @param newx Standardized predictor matrix (same standardization as training Xs)
#' @param Indi Optional indicator variables for new data
#' @param threshold Classification threshold (default 0.5)
#' 
#' @return Binary predictions (0/1)
#' @export
predict_class <- function(fit, newx, Indi = NULL, threshold = 0.5) {
  stopifnot(inherits(fit, "GeMLR_fit"))
  Xs_new <- as.matrix(newx)
  
  # Apply varreg if it was used in training
  if (!is.null(fit$varreg)) {
    X_reg <- Xs_new[, fit$varreg, drop = FALSE]
  } else {
    X_reg <- Xs_new
  }
  
  Xlogit <- if (is.null(Indi)) X_reg else cbind(X_reg, as.matrix(Indi))
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