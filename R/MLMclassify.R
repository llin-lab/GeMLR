#' MLM classify
#' @description
#' classify based on a given MLM model
#'
#' @param a the weight of each cluster
#' @param mu the mean of each cluster
#' @param sigma the sigma matrix of each cluster
#' @param beta the logistic coefficients of each cluster
#' @param X all non-dummy feature variables
#' @param Xlogit all feature variables
#'
#' @return P(Y=1|given X) and P(Z=c cluster|given X)
#' @export
MLMclassify <- function(a, mu, sigma, beta, X, Xlogit) {
  numdata <- nrow(X)
  dim <- ncol(X)
  dimlogit <- ncol(Xlogit)
  numcmp <- length(a)

  sigmainv <- array(0, dim = c(dim, dim, numcmp))
  sigmadetsqrt <- numeric(numcmp)

  for (j in 1:numcmp) {
    sigmainv[,,j] <- solve(sigma[,,j])
    sigmadetsqrt[j] <- sqrt(det(sigma[,,j]))
  }

  pij <- matrix(0, nrow = numdata, ncol = numcmp)
  pyij <- matrix(0, nrow = numdata, ncol = numcmp)
  pyi <- numeric(numdata)

  for (i in 1:numdata) {
    tmp <- 0.0
    for (j in 1:numcmp) {
      pij[i, j] <- a[j] / sigmadetsqrt[j] * exp(-0.5 * t(as.numeric(X[i,] - mu[,j])) %*% sigmainv[,,j] %*% as.numeric(X[i,] - mu[,j]))
      if (pij[i, j] >= 0) {
        tmp <- tmp + pij[i, j]
      } else {
        warning(sprintf('Numerical error when computing pij: pij(%d, %d)=%f', i, j, pij[i, j]))
        warning(sprintf('sigmadetsqrt(%d)=%f', j, sigmadetsqrt[j]))
        print(sprintf('Data point %d:', i))
        print(X[i,])
        stop('MLMclassify: Joint density of X and the component should be nonnegative')
      }
    }

    for (j in 1:numcmp) {
      if (tmp > 0) {
        pij[i, j] <- pij[i, j] / tmp
      } else {
        pij[i, j] <- 1 / numcmp
      }
    }

    for (j in 1:numcmp) {
      v1 <- exp(sum(Xlogit[i,] * beta[2:(dimlogit+1), j]) + beta[1, j])
      if (v1 > 1.0e+10) {
        pyij[i, j] <- 1
      } else {
        pyij[i, j] <- v1 / (1.0 + v1)
      }
    }

    pyi[i] <- sum(pij[i,] * pyij[i,])
  }

  return(list(pyi = pyi, pij = pij))
}
