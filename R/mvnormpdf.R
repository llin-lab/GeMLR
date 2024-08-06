#' mvnormpdf
#'
#' @param x p.n array of n values of p-dim MV normal
#' @param mu column p vector mean
#' @param Sigma p.p variance matrix
#'
#' @return the value of pdf
#'
mvnormpdf <- function(x, mu, Sigma) {
  p <- nrow(x)
  n <- ncol(x)
  C <- chol(Sigma)
  e <- solve(t(C)) %*% (x - matrix(mu, nrow = length(mu), ncol = n, byrow = TRUE))
  pdf <- exp(-rowSums(e * e) / 2) / (prod(diag(C)) * (2 * pi)^(p / 2))
  return(pdf)
}
