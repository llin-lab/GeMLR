#' init emMLM
#' @description use EM algorithm to initiate a model.
#'
#' @param X all non-dummy features that are used in logistic model
#' @param Xlogit include all features and dummy feature
#' @param Y the response variable
#' @param MLMoption some necessary variables
#'
#' @return a list of initialization informations of each cluster
#' @export
#'
initemMLM <- function(X, Xlogit, Y, MLMoption) {
  numdata <- nrow(X)
  dim <- ncol(X)
  dimlogit <- ncol(Xlogit)
  numcmp <- MLMoption$numcmp[1]

  if (numcmp > numdata) {
    stop(sprintf('Abort: number of components=%d is larger than data size %d\n', numcmp, numdata))
  }

  kmseed <- MLMoption$kmseed
  Xvar <- apply(X, 2, var)
  if (mean(Xvar) < 1.0e-6) {
    cat(sprintf('Warning: average variance is very small: %f, may lead to singular matrix\n', mean(Xvar)))
  }
  for (j in 1:dim) {
    if (Xvar[j] < 1.0e-6) {
      cat(sprintf('Warning: variance of dimension %d is very small: %f, may lead to singular matrix\n', j, Xvar[j]))
    }
  }

  if (nrow(MLMoption$InitCluster) == 0) {
    nloops <- 10
  } else {
    nloops <- 0
  }

  for (ii in 1:nloops) {
    km_result <- km(X, numcmp, 1.0e-4, kmseed)
    cdbk <- km_result$cdbk
    cd1 <- km_result$ind
    ndatpercls <- rep(0, numcmp)
    for (i in 1:numdata) {
      ndatpercls[cd1[i]] <- ndatpercls[cd1[i]] + 1
    }

    k <- sum(ndatpercls == 0)
    if (k == 0) {
      break
    }

    cat(sprintf('Warning initemMLM: kmeans generated empty cluster. Rerun kmeans using seed: %d\n', kmseed + ii))

    if (ii == nloops) {
      cat(sprintf('Warning initemMLM: kmeans generated empty cluster %d times.\n', nloops))
      cat(sprintf('Warning: failure to generate kmeans with NO empty clusters\n'))
      cat(sprintf('Split clusters to generate %d clusters from %d data points\n', numcmp, numdata))
      km_result <- kmRmEmpty(X, numcmp, cdbk)
      cdbk <- km_result$cdbk
      cd1 <- km_result$ind
    }

    kmseed <- kmseed + ii
  }

  if (nloops == 0) {
    cd1 <- MLMoption$InitCluster[, 1]

    if (max(cd1) > numcmp) {
      cat(sprintf('Warning: the number of components %d in MLMoption is smaller than %d appeared in the initial cluster labels\n', numcmp, max(cd1)))
      cat('Abort initemMLMcls\n')
      return(NULL)
    }

    ndatpercls <- rep(0, numcmp)
    for (i in 1:numdata) {
      ndatpercls[cd1[i]] <- ndatpercls[cd1[i]] + 1
    }

    k <- sum(ndatpercls == 0)
    if (k > 0) {
      cat('Warning: there are one or more empty clusters in the input initial cluster labels\n')
      cat('Abort initemMLMcls\n')
      return(NULL)
    }
  }

  rsigma2 <- rep(0, numcmp)
  ainit <- ndatpercls / numdata
  muinit <- matrix(0, nrow = dim, ncol = numcmp)
  for (i in 1:numdata) {
    muinit[, cd1[i]] <- as.matrix(muinit[, cd1[i]] + X[i, ])
  }

  for (j in 1:numcmp) {
    muinit[, j] <- as.matrix(muinit[, j] / ndatpercls[j])
  }

  sigmainit <- array(0, dim = c(dim, dim, numcmp))
  for (i in 1:numdata) {
    sigmainit[, , cd1[i]] <- sigmainit[,,cd1[i]] + as.numeric((X[i, ] - muinit[, cd1[i]]))%*%t(as.numeric((X[i, ] - muinit[, cd1[i]])))
  }

  for (j in 1:numcmp) {
    if (ndatpercls[j] > 0) {
      sigmainit[, , j] <- as.matrix(sigmainit[, , j] / ndatpercls[j])
    }
    sigmainit[, , j] <- checksingular(sigmainit[, , j], Xvar, 0.05)
  }

  beta <- matrix(0, nrow = dimlogit + 1, ncol = numcmp)
  for (j in 1:numcmp) {
    numdata_small <- sum(cd1 == j)
    Ysmall <- Y[cd1 == j]
    Xsmall <- Xlogit[cd1 == j, ]
    counts <- table(Ysmall)
    count_0 <- ifelse("0" %in% names(counts), counts["0"], 0)
    count_1 <- ifelse("1" %in% names(counts), counts["1"], 0)
    min_count <- min(count_0, count_1)
    if (min_count >= 2) {
      lasso_result <- glmnet(Xsmall, Ysmall, MLMoption$DISTR, lambda = MLMoption$lambdaLasso[j], alpha = MLMoption$AlphaLasso)
      beta_nointercept <- as.matrix(lasso_result$beta)
      beta_intercept <- lasso_result$a0
      beta[, j] <- rbind(beta_intercept, beta_nointercept)
    } else {
      beta[, j] <- glmIntercept(Ysmall, rep(1, length(Ysmall)), dimlogit, MLMoption$DISTR)
    }

    if (MLMoption$DISTR == 'normal') {
      rsigma2[j] <- sum((Ysmall - (Xsmall %*% beta[2:(dimlogit + 1), j] + beta[1, j]))^2)
    }
  }

  v1 <- sum(rsigma2)
  rsigma2 <- rep(v1 / numdata, length(rsigma2))

  return(list(muinit = muinit, sigmainit = sigmainit, ainit = ainit, beta = beta, rsigma2 = rsigma2))
}
