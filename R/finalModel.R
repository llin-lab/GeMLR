#' final Model
#' @description
#' use AUC to decide the model
#'
#' @param cvAUCfinal the AUC results in cross-verified
#' @param ncmp the number of clusters
#' @param nseeds the number of random seeds
#' @param rangeSeed the largest number of random seeds
#' @param vargmm the number of variables that are used in gmm model
#' @param Y1 the response variable
#' @param X1s the standardized features
#' @param Indi the dummy variables
#' @param MLMoption all necessary variables
#'
#' @return a list of clustering results
#' @export

finalModel <- function(cvAUCfinal, ncmp, nseeds, rangeSeed, vargmm, Y1, X1s, Indi, MLMoption) {
  dimgmm <- length(vargmm)

  maxvf <- max(colMeans(cvAUCfinal))
  maxcmp <- which.max(colMeans(cvAUCfinal))
  MLMoption$numcmp <- ncmp[maxcmp]
  MLMoption$lambdaLasso <- rep(MLMoption$lambdaLasso[1], ncmp[maxcmp])
  set.seed(9)
  rseeds <- sample(1:rangeSeed, nseeds, replace = FALSE)

  est_result <- estimateBestSD(X1s[, vargmm], cbind(X1s, Indi), Y1, MLMoption, rseeds)
  beta <- est_result$beta

  gmm_result <- GMMFormatConvert(dimgmm, est_result$c)
  a2 <- gmm_result$a
  mu2 <- gmm_result$mu
  sigma2 <- gmm_result$sigma

  classify_result <- MLMclassify(a2, mu2, sigma2, beta, X1s[, vargmm], cbind(X1s, Indi))
  pij <- classify_result$pij
  clusterid <- apply(pij, 1, which.max)

  return(list(beta = beta, clusterid = clusterid, a2 = a2, mu2 = mu2, sigma2 = sigma2))
}
