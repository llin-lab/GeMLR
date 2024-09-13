#' classify based on a given MLM model
#'
#'
#' @return classification
#' @export
#'
KMclassify <- function(a, mu, sigma, beta, X, Xlogit) {
  # get the data and dimension
  numdata <- nrow(X)
  dim <- ncol(X)
  dimlogit <- ncol(Xlogit)
  numcmp <- length(a)

  # initialize
  dij <- matrix(0, nrow = numdata, ncol = numcmp)
  pyi <- numeric(numdata)
  pij <- matrix(0, nrow = numdata, ncol = numcmp)
  loglike <- 0.0

  # calculate the classification
  for (i in 1:numdata) {
    for (j in 1:numcmp) {
      dij[i, j] <- sum((X[i, ] - mu[, j])^2)
    }

    mv <- min(dij[i, ])
    jj <- which.min(dij[i, ])

    pij[i, jj] <- 1

    v1 <- exp(Xlogit[i, ] %*% beta[2:(dimlogit + 1), jj] + beta[1, jj])
    if (v1 > 1.0e+10) {
      pyi[i] <- 1
    } else {
      pyi[i] <- v1 / (1.0 + v1)
    }
  }

  return(list(pyi = pyi, pij = pij))
}
