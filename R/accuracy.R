#' accuracy
#' @description
#' Compute accuracy for the selected dimensions (index) for a specific component comp0,in normal mixture model.
#'
#'
accuracy <- function(pi, mu, Sigma, index) {
  num <- nrow(index)
  k <- ncol(mu)
  D <- array(0, dim = c(k, k, num))

  pi <- t(pi)
  Acc <- matrix(0, k, num)
  delta_ncnc <- matrix(0, k, num)

  for (s in 1:num) {
    index1 <- index[s, ]
    index1 <- index1[index1 != 0]

    for (i in 1:k) {
      for (j in i:k) {
        D[i, j, s] <- mvnormpdf(mu[index1, j, drop = FALSE], mu[index1, i, drop = FALSE], Sigma[index1, index1, i] + Sigma[index1, index1, j])
        D[j, i, s] <- D[i, j, s]
      }
    }

    E <- D[, , s] - diag(diag(D[, , s]))
    dcc <- 1 / (1 - pi) * (E %*% pi)
    Delta_c <- dcc / diag(D[, , s])

    for (i in 1:k) {
      Dsub <- D[, , s]
      Dsub <- Dsub[-i, -i]
      pisub <- pi
      pisub[i] <- 0

      delta_ncnc[i, s] <- (1 / (1 - pi[i])^2) * sum(Dsub * outer(pisub, pisub))
    }

    Delta_nc <- dcc / delta_ncnc[, s]
    Acc[, s] <- (pi * pi) / (pi + (1 - pi) * Delta_c) + (1 - pi) * (1 - pi) / (1 - pi + pi * Delta_nc)
  }

  list(Acc = Acc, D = D)
}
