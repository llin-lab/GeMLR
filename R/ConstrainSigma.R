#' ConstrainSigma
#'
#' @return sigmaout
#' @export
ConstrainSigma <- function(a, sigma, dim, numcmp, Xvar, MLMoption) {
  code <- MLMoption$constrain
  ncode <- switch(code,
                  'N' = 0,
                  'EI' = 1,
                  'VI' = 2,
                  'EEE' = 3,
                  'VVV' = 4,
                  'DIA' = 5,
                  'DIAE' = 6,
                  'DIAS' = 7,
                  'EEV' = 8,
                  'VEV' = 9,
                  0  # default case
  )
  
  if (ncode == 0 || ncode == 4) {
    return(sigma)
  }
  
  if (ncode == 5) {
    sigmaout <- array(0, dim = dim(sigma))
    for (j in 1:numcmp) {
      # FIX: Handle 1x1 matrix case
      if (dim == 1) {
        sigmaout[,,j] <- matrix(sigma[1,1,j], 1, 1)
      } else {
        sigmaout[,,j] <- diag(diag(sigma[,,j]))
      }
    }
    return(sigmaout)
  }
  
  if (ncode == 7) {
    v1 <- MLMoption$diagshrink
    v1 <- pmin(pmax(v1, 0.0), 1.0)
    sigmaout <- array(0, dim = dim(sigma))
    for (j in 1:numcmp) {
      # FIX: Handle 1x1 matrix case
      if (dim == 1) {
        sigmaout[,,j] <- matrix(sigma[1,1,j], 1, 1)
      } else {
        sigmaout[,,j] <- sigma[,,j] * (1 - v1) + v1 * diag(diag(sigma[,,j]))
      }
    }
    return(sigmaout)
  }
  
  sigmaave <- matrix(0, nrow = dim, ncol = dim)
  for (j in 1:numcmp) {
    sigmaave <- sigmaave + a[j] * sigma[,,j]
  }
  
  if (ncode == 6) {
    sigmaout <- array(0, dim = dim(sigma))
    for (j in 1:numcmp) {
      # FIX: Handle 1x1 matrix case
      if (dim == 1) {
        sigmaout[,,j] <- matrix(sigmaave[1,1], 1, 1)
      } else {
        sigmaout[,,j] <- diag(diag(sigmaave))
      }
    }
    return(sigmaout)
  }
  
  if (ncode == 1) {
    sigmaout <- array(0, dim = dim(sigma))
    for (j in 1:numcmp) {
      sigmaout[,,j] <- mean(diag(sigmaave)) * diag(rep(1, dim))
    }
    return(sigmaout)
  }
  
  if (ncode == 2) {
    sigmaout <- array(0, dim = dim(sigma))
    for (j in 1:numcmp) {
        sigmaout[,,j] <- mean(diag(matrix(sigma[,,j], dim, dim))) * diag(dim)
    }
    return(sigmaout)
  }
  
  if (ncode == 3) {
    sigmaout <- array(0, dim = dim(sigma))
    for (j in 1:numcmp) {
      sigmaout[,,j] <- sigmaave
    }
    return(sigmaout)
  }
  
  vscale <- numeric(numcmp)
  V <- array(0, dim = c(dim, dim, numcmp + 1))
  D <- array(0, dim = c(dim, dim, numcmp + 1))
  
  for (j in 1:numcmp) {
    eig_result <- eigen(matrix(sigma[,,j], dim, dim), symmetric = TRUE)
    V[,,j] <- eig_result$vectors
    D[,,j] <- diag(eig_result$values, nrow = dim)
    vscalelog <- sum(log(diag(matrix(D[,,j], dim, dim))))
    if (any(diag(matrix(D[,,j], dim, dim)) <= 0)) {
      stop("ConstrainSigma: Error\n")
    }
    vscale[j] <- exp(vscalelog / dim)
    D[,,j] <- D[,,j] / vscale[j]
  }
  
  eig_result <- eigen(sigmaave, symmetric = TRUE)
  V[,,numcmp + 1] <- eig_result$vectors
  D[,,numcmp + 1] <- diag(eig_result$values, nrow = dim)
  vscalelog <- sum(log(diag(matrix(D[,,numcmp + 1], dim, dim))))
  if (any(diag(matrix(D[,,numcmp + 1], dim, dim)) <= 0)) {
    stop("ConstrainSigma: Error\n")
  }
  vscaleave <- exp(vscalelog / dim)
  D[,,numcmp + 1] <- D[,,numcmp + 1] / vscaleave
  
  sigmaout <- array(0, dim = dim(sigma))
  if (ncode == 8) {
    for (j in 1:numcmp) {
        sigmaout[,,j] <- vscaleave * matrix(V[,,j], dim, dim) %*% matrix(D[,,numcmp + 1], dim, dim) %*% t(matrix(V[,,j], dim, dim))
    }
  }
  
  if (ncode == 9) {
    for (j in 1:numcmp) {
        sigmaout[,,j] <- vscale[j] * matrix(V[,,j], dim, dim) %*% matrix(D[,,numcmp + 1], dim, dim) %*% t(matrix(V[,,j], dim, dim))
    }
  }
  
  for (j in 1:numcmp) {
    sigmaout_j <- as.matrix(sigmaout[,,j])  # Prevent dimension drop
    sigmaout[,,j] <- checksingular(sigmaout_j, Xvar, 0.05)
  }
  
  return(sigmaout)
}
