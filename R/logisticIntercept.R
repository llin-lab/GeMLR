#' Compute the logistic model's intercept
#'
#' @param Y the response variable
#' @param wt weights
#' @param dim the dimension of dataset(number of features)
#'
#' @return the beta coefficients
#' @export
#'
logisticIntercept <- function(Y, wt, dim) {
  vec <- which(Y == 1)

  if (sum(wt) > 0) {
    probclass1 <- sum(wt[vec]) / sum(wt)
  } else {
    probclass1 <- 0.5
  }

  if (probclass1 >= 0.9999) {
    probclass1 <- 0.9999
  } else if (probclass1 <= 0.0001) {
    probclass1 <- 0.0001
  }

  beta <- rep(0, dim + 1)
  beta[1] <- log(probclass1 / (1 - probclass1))  # Intercept

  return(beta)
}
