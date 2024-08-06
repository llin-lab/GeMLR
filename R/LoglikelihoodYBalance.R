#' LoglikehoodYBalance
#' @description
#'  Compute the conditional log likelihood of Y given X.As the classes may be very unbalanced, adjust the likelihood by assuming the class proportions are given in input w. If w is set to the empirical proportions of the classes, the output is simply the conditional log likelihood of the observed Y.
#'
#'
#' @param pyi P(Y=1|X,Z)
#' @param Y the resonse variables
#' @param w the weigh
#'
#' @return the loglike Y
#' @export
#'
LoglikehoodYBalance <- function(pyi, Y, w) {

  numdata <- length(pyi)

  loglike <- numeric(2)
  for (i in 1:numdata) {
    if (Y[i] == 1) {
      if (pyi[i] > 0.0001) {
        loglike[1] <- loglike[1] + log(pyi[i])
      } else {
        loglike[1] <- loglike[1] + log(0.0001)
      }
    } else {
      if (1 - pyi[i] > 0.0001) {
        loglike[2] <- loglike[2] + log(1 - pyi[i])
      } else {
        loglike[2] <- loglike[2] + log(0.0001)
      }
    }
  }

  n1 <- sum(Y == 1)
  n <- length(Y)

  if (sum(w) <= 0) {
    stop(sprintf("LoglikehoodYBalance: Input weight w has non-positive sum, [%f, %f]", w[1], w[2]))
  }

  w <- w / sum(w)
  if (n1 > 0 & n - n1 > 0) {
    loglikeY <- n * (loglike[1] / n1 * w[1] + loglike[2] / (n - n1) * w[2])
  } else {
    loglikeY <- sum(loglike)
    warning("Warning LoglikehoodYBalance: Input Y has only one class label, loglikelhood of Y not balanced")
  }

  return(loglikeY)
}
