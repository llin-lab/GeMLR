#' final Model
#' @description
#' use AUC to decide the model
#'
#' @param cvAUCfinal the AUC results in cross-verified
#' @param ncmp the number of clusters
#' @param nseeds the number of random seeds
#' @param rangeSeed the largest number of random seeds
#' @param vargmm the index of variables that are used in gmm model
#' @param Y the response variable
#' @param Xs the standardized independent variables
#' @param Indi the dummy variables
#' @param MLMoption all necessary variables
#'
#' @return a list of clustering results
#' @export

finalModel <- function(cvAUCfinal, ncmp, nseeds, rangeSeed, vargmm, Y, Xs, Indi, MLMoption) {
  dimgmm <- length(vargmm)

  maxvf <- max(colMeans(cvAUCfinal))
  maxcmp <- which.max(colMeans(cvAUCfinal))
  MLMoption$numcmp <- ncmp[maxcmp]
  MLMoption$lambdaLasso <- rep(MLMoption$lambdaLasso[1], ncmp[maxcmp])
  set.seed(9)
  rseeds <- sample(1:rangeSeed, nseeds, replace = FALSE)

  est_result <- estimateBestSD(Xs[, vargmm], cbind(Xs, Indi), Y, MLMoption, rseeds)

  beta_rownames <- c('Intercept',colnames(cbind(Xs, Indi)))
  beta <- est_result$beta
  if (nrow(beta) == length(beta_rownames)) {
    rownames(beta) <- beta_rownames
  } else if (nrow(beta) < length(beta_rownames)) {
    rownames(beta) <- beta_rownames[1:nrow(beta)]
  } else {
    warning("beta has more rows than beta_rownames. Row names not assigned.")
  }
  colnames(beta) <- paste0("Cluster ", 1:ncol(beta))

  gmm_result <- GMMFormatConvert(dimgmm, est_result$c)
  a2 <- gmm_result$a
  mu2 <- gmm_result$mu
  sigma2 <- gmm_result$sigma

  classify_result <- MLMclassify(a2, mu2, sigma2, beta, Xs[, vargmm], cbind(Xs, Indi))
  pij <- classify_result$pij
  clusterid <- apply(pij, 1, which.max)

  return(list(beta = beta, clusterid = clusterid, a2 = a2, mu2 = mu2, sigma2 = sigma2))
}
