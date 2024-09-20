#' Use eigen value decomposition to solve squre root
#' @param A matrix
#'
#' @return the square root matrix
#' @export
#'
sqrtm_eig <- function(A) {
  # Get matrix dimension
  dim_A <- nrow(A)

  # Force symmetry
  A <- (A + t(A)) / 2

  # Eigen decomposition
  eig_decomp <- eigen(A)
  V <- eig_decomp$vectors
  D <- diag(eig_decomp$values)

  # Defend against precision error causing negative eigenvalues
  for (k in 1:dim_A) {
    if (D[k, k] < 0) {
      D[k, k] <- 0
    }
  }

  # Return the square root matrix
  B <- V %*% sqrt(D) %*% t(V)
  return(B)
}
