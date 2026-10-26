#' mvnormpdf
#'
#' @param x p x n array, one sample per COLUMN
#' @param mu length-p vector mean (or p x 1 matrix)
#' @param Sigma p.p variance matrix
#'
#' @return the value of pdf
#'
mvnormpdf <- function(x, mu, Sigma) {
  p <- nrow(x)
  n <- ncol(x)

  C <- chol(Sigma)

  e <- solve(t(C)) %*% (x - matrix(mu, nrow = length(mu), ncol = n))
  pdf <- exp(-colSums(e * e) / 2) / (prod(diag(C)) * (2 * pi)^(p / 2))

  return(pdf)
}
