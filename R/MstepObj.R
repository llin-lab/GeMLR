#' Compute the Q(theta'|theta) function in the EM algorithm
#'
#' @param dim the number of features
#' @param numdata the number of samples
#' @param supp mu and sigma informations about the GMM model
#' @param ww the weights of old and new model
#' @param X the features of each sample
#' @param pij the posterior matrix pij should be pre-computed based on previous model, not the updated model
#'
#' @return the likelihood
#' @export
#'

MstepObj <- function(dim, numdata, supp, ww, X, pij) {

  numcmp <- length(ww)
  mu <- supp[1:dim, 1:numcmp]
  a <- ww

  sigma <- array(0, dim = c(dim, dim, numcmp))
  for (k2 in 1:numcmp) {
    sigma[, , k2] <- matrix(supp[(dim + 1):(dim + dim * dim), k2], nrow = dim, ncol = dim)
  }

  sigmainv <- array(0, dim = c(dim, dim, numcmp))
  sigmadetsqrt <- numeric(numcmp)

  for (j in 1:numcmp) {
    sigmainv[, , j] <- solve(sigma[, , j])
    sigmadetsqrt[j] <- sqrt(det(sigma[, , j]))
  }

  likelihood <- 0.0

  for (i in 1:numdata) {
    v1 <- 0.0
    v2 <- 0.0

    for (j in 1:numcmp) {
      v1 <- v1 + pij[i, j] * log(a[j])
      v3 <- -0.5 * t(X[, i] - mu[, j]) %*% sigmainv[, , j] %*% (X[, i] - mu[, j])
      v3 <- v3 - log(sigmadetsqrt[j]) - (dim / 2) * log(2 * pi)
      v2 <- v2 + pij[i, j] * v3
    }

    likelihood <- likelihood + v1 + v2
  }

  return(likelihood)
}
