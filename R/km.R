#' km
#' @description
#' use k-means algorithm to initiate the Gussian mixture model
#'
#' @param dat the input data with every sample stored in one row
#' @param nkm the number of clusters prespecified
#' @param thred the threshold for determining when to stop the iteration in kmeans. Suggested value for thred: 1.0e-4 to 1.0e-6
#' @param kmseed a random number seed that specifies a way to assign random numbers to each sample
#'
#' @return the center of each cluster and the identity of each sample
#' @export
#'
km <- function(dat, nkm, thred, kmseed) {
  len <- nrow(dat)
  dim <- ncol(dat)

  set.seed(kmseed)
  b <- runif(len)
  ind2 <- order(b)

  cdbk <- matrix(0, nkm, dim)

  for (i in 1:nkm) {
    cdbk[i, ] <- as.matrix(dat[ind2[i], ])
  }

  ind <- rep(0, len)
  aa <- rep(0, nkm)

  stdlist <- apply(dat, 2, sd)
  mvlist <- colMeans(dat)
  dist <- sum(stdlist^2) * len * 10

  done <- 0
  distlist <- numeric()

  while (done != -1) {
    newcdbk <- matrix(0, nkm, dim)
    newct <- rep(0, nkm)
    newdist <- 0

    for (i in 1:len) {
      for (j in 1:nkm) {
        aa[j] <- sum((dat[i, ] - cdbk[j, ])^2)
      }
      minv <- min(aa)
      ind[i] <- which.min(aa)
      newcdbk[ind[i], ] <- as.matrix(newcdbk[ind[i], ] + dat[i, ])
      newct[ind[i]] <- newct[ind[i]] + 1
      newdist <- newdist + minv
    }

    for (j in 1:nkm) {
      if (newct[j] > 0) {
        newcdbk[j, ] <- as.matrix(newcdbk[j, ] / newct[j])
      } else {
        newcdbk[j, ] <- as.matrix(cdbk[j, ])
      }
    }

    done <- done + 1
    distlist[done] <- newdist

    if ((dist - newdist) / dist < thred) {
      done <- -1
    }

    dist <- newdist
    cdbk <- newcdbk
  }

  ind <- rep(1, len)
  aa <- rep(0, nkm)

  for (i in 1:len) {
    for (j in 1:nkm) {
      aa[j] <- sum((dat[i, ] - cdbk[j, ])^2)
    }
    minv <- min(aa)
    ind[i] <- which.min(aa)
  }

  return(list(cdbk = cdbk, distlist = distlist, ind = ind))
}
