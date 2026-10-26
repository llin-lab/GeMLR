#' accuracy
#' @description
#' Compute accuracy for the selected dimensions (index) for a specific component comp0,in normal mixture model.
#'
#'
accuracy <- function(pi, mu, Sigma, index) {
  num <- nrow(index)
  k <- ncol(mu)
  D <- array(0, dim = c(k, k, num))

  pi <- as.numeric(pi)
  Acc <- matrix(0, k, num)
  delta_ncnc <- matrix(0, k, num)

  for (s in seq_len(num)) {
    index1 <- index[s, ]
    index1 <- index1[index1 != 0]

    for (i in 1:k) {
      for (j in i:k) {
        D[i, j, s] <- mvnormpdf(mu[index1, j, drop = FALSE], mu[index1, i, drop = FALSE], Sigma[index1, index1, i] + Sigma[index1, index1, j])
        D[j, i, s] <- D[i, j, s]
      }
    }

    Ds <- matrix(D[, , s], k, k)
    E <- Ds - diag(diag(Ds), nrow = k)
    dcc <- as.numeric(E %*% pi) / (1 - pi)
    Delta_c <- dcc / diag(Ds)

    for (i in 1:k) {
      Dsub <- Ds[-i, -i, drop = FALSE]
      pisub <- pi[-i]

      delta_ncnc[i, s] <- (1 / (1 - pi[i])^2) * sum(Dsub * outer(pisub, pisub))
    }

    Delta_nc <- dcc / delta_ncnc[, s]
    Acc[, s] <- (pi * pi) / (pi + (1 - pi) * Delta_c) + (1 - pi) * (1 - pi) / (1 - pi + pi * Delta_nc)
  }

  list(Acc = Acc, D = D)
}
