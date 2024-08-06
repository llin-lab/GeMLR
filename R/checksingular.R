#' check singular
#' @description
#' Check whether covariance matrix is ill-conditioned or nearly singular. If so, regularize covariance by changing the diagonal elements.
#'
#' @param sigma the sigma matrix
#' @param Xvar the gmm features variables
#' @param shrinkrate the rate of shrinking
#'
#' @return the new singular sigma matrix
#' @export
#'
checksingular <- function(sigma, Xvar, shrinkrate) {
  dim <- nrow(sigma)
  sigmanew <- sigma

  thred <- mean(Xvar) * 1.0e-4
  t <- mean(diag(sigma))

  if (rcond(sigma) < 1.0e-8 || t < thred) {
    v1 <- sum(diag(sigma)) / dim
    for (i in 1:dim) {
      sigma[i, i] <- sigma[i, i] + shrinkrate * v1
    }
  } else {
    return(sigmanew)
  }
  sigmanew <- sigma

  t <- mean(diag(sigma))

  if (rcond(sigma) < 1.0e-8 || t < thred) {
    v1 <- mean(Xvar)
    for (i in 1:dim) {
      sigma[i, i] <- sigma[i, i] + shrinkrate * v1
    }
  } else {
    return(sigmanew)
  }
  sigmanew <- sigma

  t <- mean(diag(sigma))

  if (rcond(sigma) < 1.0e-8 || t < thred) {
    warning('Unable to solve singular covariance matrix, recommend to modify data')
    print(sigma)
    stop('checksingular: failed to modify sigma to be away from singular')
  }

  return(sigmanew)
}
