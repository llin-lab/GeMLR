#' final Model
#' @description
#' use AUC to decide the model
#'
#' @param cvAUCfinal the AUC results in cross-verified
#' @param ncmp the number of clusters
#' @param nseeds the number of random seeds
#' @param rangeSeed the largest number of random seeds
#' @param vargmm the index of variables that are used in gmm model
#' @param vlasso Lambda for Lasso penalty
#' @param Y the response variable
#' @param Xs the standardized independent variables
#' @param X the raw input of independent variables
#' @param Indi the dummy variables
#' @param varreg (optional) Variables for logistic regression. If NULL, uses vargmm
#' @param alphaLasso Elastic net parameter (default 0.8)
#' @param stopratio Convergence threshold (default 1.0e-5)
#' @param verbose Verbosity flag (default 1)
#' @param constrain Covariance constraint (default "DIAS")
#' @param diagshrink Diagonal shrinkage (default 0.9)
#' @param kappa Instance weighting (default -1)
#' @param AUC Use AUC for selection (default 1)
#' @param DISTR Distribution family (default "binomial")
#'
#' @return a list of clustering results
#' @export

finalModel <- function(cvAUCfinal, ncmp, nseeds, rangeSeed, 
                       vargmm, vlasso, Y, Xs, X, Indi, 
                       varreg = NULL,
                       alphaLasso = 0.8,
                       stopratio = 1.0e-5,
                       verbose = 1,
                       constrain = "DIAS",
                       diagshrink = 0.9,
                       kappa = -1,
                       AUC = 1,
                       DISTR = "binomial") {
  dimgmm <- length(vargmm)
  
  maxvf <- max(colMeans(cvAUCfinal))
  maxcmp <- which.max(colMeans(cvAUCfinal))
  
  # Determine which variables to use for logistic regression
  # If varreg is NULL, use vargmm (same variables for GMM and LR)
  if (is.null(varreg)) {
    varreg_to_use <- vargmm
    if (verbose > 0) {
      cat("varreg not specified. Using vargmm for logistic regression.\n")
    }
  } else {
    varreg_to_use <- varreg
  }
  
  # Prepare logistic regression data
  if (is.null(Indi)) {
    Xlogit <- Xs[, varreg_to_use, drop = FALSE]
    X_raw_logit <- X[, varreg_to_use, drop = FALSE]
  } else {
    Xlogit <- cbind(Xs[, varreg_to_use, drop = FALSE], Indi)
    X_raw_logit <- cbind(X[, varreg_to_use, drop = FALSE], Indi)
  }
  
  # Create MLMoption internally
  MLMoption <- init_MLMoption(
    alphaLasso = alphaLasso,
    vlasso = vlasso,
    numcmp = ncmp[maxcmp],
    stopratio = stopratio,
    verbose = verbose,
    minloop = 3,
    maxloop = 50,
    constrain = constrain,
    diagshrink = diagshrink,
    kmseed = 0,
    algorithm = 1,
    kappa = kappa,
    AUC = AUC,
    DISTR = DISTR,
    NOEM = 0,
    Yalpha = 1.0
  )
  MLMoption$lambdaLasso <- rep(vlasso, ncmp[maxcmp])
  
  set.seed(9)
  rseeds <- sample(1:rangeSeed, nseeds, replace = FALSE)
  
  # GMM clustering uses vargmm, LR uses varreg_to_use
  est_result <- estimateBestSD(Xs[, vargmm, drop = FALSE], Xlogit, Y, MLMoption, rseeds)
  
  beta_rownames <- c('Intercept', colnames(Xlogit))
  beta <- est_result$beta
  if (nrow(beta) == length(beta_rownames)) {
    rownames(beta) <- beta_rownames
  } else if (nrow(beta) < length(beta_rownames)) {
    rownames(beta) <- beta_rownames[1:nrow(beta)]
  } else {
    warning("beta has more rows than beta_rownames. Row names not assigned.")
  }
  colnames(beta) <- paste0("Cluster ", 1:ncol(beta))
  
  # Add pure LR model using raw X
  set.seed(1)
  cv_fit <- suppressWarnings(cv.glmnet(data.matrix(X_raw_logit), Y, alpha = 1, family = "binomial", nfolds = 5))
  B <- coef(cv_fit, s = "lambda.min")
  
  beta <- cbind(beta, B)
  colnames(beta)[ncol(beta)] <- "LR"
  
  gmm_result <- GMMFormatConvert(dimgmm, est_result$c)
  a2 <- gmm_result$a
  mu2 <- gmm_result$mu
  sigma2 <- gmm_result$sigma
  
  classify_result <- MLMclassify(a2, mu2, sigma2, beta, Xs[, vargmm, drop = FALSE], Xlogit)
  pij <- classify_result$pij
  clusterid <- apply(pij, 1, which.max)
  
  return(list(beta = as.matrix(beta), clusterid = clusterid, a2 = a2, mu2 = mu2, sigma2 = sigma2))
}