#' glm Intercept
#' @description
#' Fit a linear model with only intercept term. This is used when the number of data points for estimation is extremely small.
#'
#' @param Y the respinse variable
#' @param wt the weight
#' @param dim the number of features
#' @param DISTR default = binominal
#'
#' @return the intercept coefficient of the elastic net regression
#' @export
glmIntercept <- function(Y, wt, dim, DISTR) {
  if (DISTR == 'binomial') {
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

    beta <- matrix(0,nrow=dim + 1,ncol=1)
    beta[1,1] <- log(probclass1 / (1 - probclass1))
    return(beta[,1])
  }

  stop('In glmIntercept: DISTR improper\n')
}
