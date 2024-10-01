#' em_MLM_maw
#' @description
#' Use maw distance to advance the iteration in em_MLM.
#'
#' @param X all non-dummy features that are used in the logistic model
#' @param Xlogit include the dummy variables
#' @param Y the response variable
#' @param cinit some information about the Gaussian model
#' @param betainit the coefficients in logistics model
#' @param MLMoption some necessary variabls
#' @param cprior the given template GMM model
#'
#' @return the optimized parameters
#' @export
#'
em_MLM_maw <- function(X, Xlogit, Y, cinit, betainit, MLMoption,cprior) {
  library(glmnet)

  alpha <- MLMoption$alpha
  kappa <- MLMoption$kappa
  numdata <- nrow(X)
  dim <- ncol(X)
  dimlogit <- ncol(Xlogit)
  numcmp <- length(cinit$w)

  lambdaLasso <- rep(0, numcmp)
  m <- length(MLMoption$lambdaLasso)
  k <- 1
  for (i in 1:numcmp) {
    lambdaLasso[i] <- MLMoption$lambdaLasso[k]
    k <- k + 1
    if (k > m){
      k <- 1
    }
  }

  Wi <- rep(1, numdata)

  # normalize to be sure
  cprior.w <- cprior.w / sum(cprior.w)
  stride <- c(numcmp, ncompprior)
  N <- length(stride)
  instanceW <- rep(1, N)

  Xvar <- apply(X, 2, var)
  if (mean(Xvar) < 1.0e-6) {
    cat('Warning: average variance is very small:', mean(Xvar), ', may lead to singular matrix\n')
  }
  for (j in 1:dim) {
    if (Xvar[j] < 1.0e-6) {
      cat('Warning: variance of dimension', j, 'is very small:', Xvar[j], ', may lead to singular matrix\n')
    }
  }

  X <- t(X)
  Xlogit <- t(Xlogit)
  mu <- cinit$supp[1:dim,]
  start_ind <- dim+1
  end_ind <- dim+dim*dim
  sigma <- array(cinit$supp[start_ind:end_ind,], dim = c(dim, dim, numcmp))
  a <- cinit$w
  beta <- betainit

  minloop <- max(c(MLMoption$minloop, 2))
  # maxloop <- max(c(MLMoption$maxloop, 5))
  maxloop <- 200

  oldloglike <- -1.0e+30
  oldloglikepen <- oldloglike

  mawdist <- 0

  sigmainv <- array(0, dim = c(dim, dim, numcmp))
  sigmadetsqrt <- rep(0, numcmp)
  for (j in 1:numcmp) {
    sigmainv[,,j] <- solve(sigma[,,j])
    sigmadetsqrt[j] <- sqrt(det(sigma[,,j]))
  }

  loop <- 1
  while (loop < maxloop) {
    pij <- matrix(0, nrow = numdata, ncol = numcmp)
    pyij <- matrix(0, nrow = numdata, ncol = numcmp)
    loglike <- 0.0
    for (i in 1:numdata) {
      tmp <- 0.0
      for (j in 1:numcmp) {
        v1 <- exp(sum(beta[2:(dimlogit+1),j] * Xlogit[,i]) + beta[1,j])
        if (v1 > 1.0e+10) {
          pyij[i,j] <- 1
        } else {
          pyij[i,j] <- v1 / (1.0 + v1)
        }


        pij[i,j] <- a[j] / sigmadetsqrt[j] * exp(-0.5 * t(X[,i] - mu[,j]) %*% sigmainv[,,j] %*% (X[,i] - mu[,j]))
        v5 <- Y[i] * pyij[i,j] + (1 - Y[i]) * (1 - pyij[i,j])
        pij[i,j] <- pij[i,j] * v5


        if (pij[i,j] >= 0) {
          tmp <- tmp + pij[i,j]
        } else {
          cat('Warning: numerical error when computing pij: pij(', i, ',', j, ') =', pij[i,j], '\n')
          cat('sigmadetsqrt(', j, ') =', sigmadetsqrt[j], '\n')
          cat('beta(:,', j, ')\n')
          print(beta[,j])
          cat('pyij(', i, ',', j, ') =', pyij[i,j], '\n')
          cat('Data point', i, ':\n')
          print(t(X[,i]))
          stop('em_MLM_maw: joint density of X, Y, and component should be nonnegative')
        }
      }

      for (j in 1:numcmp) {
        if (tmp > 0) {
          pij[i,j] <- pij[i,j] / tmp
        } else {
          pij[i,j] <- 1 / numcmp
        }
      }

      for (j in 1:numcmp) {
        if (!(pij[i,j] >= 0 || pij[i,j] < 0)) {
          cat('pij(', i, ',', j, ') =', pij[i,j], '\n')
          stop('em_MLM_maw: Numerical error with computing posterior pij')
        }
      }

      loglike <- loglike + (log(tmp) - dim / 2 * log(2 * pi)) * Wi[i]
    }

    pij_wt <- pij
    for (i in 1:numdata) {
      pij_wt[i,] <- pij[i,] * Wi[i]
    }

    if (kappa > 0) {
      Wiunit <- Wi / sum(Wi)
      Wientropy <- 0.0
      for (i in 1:numdata) {
        if (Wiunit[i] > 0) {
          Wientropy <- Wientropy + Wiunit[i] * log(Wiunit[i])
        }
      }
    } else {
      Wientropy <- 0
    }

    if (MLMoption$algorithm == 1) {
      penbeta <- MLMoption$AlphaLasso * sum(abs(beta[2:(dimlogit+1),])) + (1 - MLMoption$AlphaLasso) * sum(beta[2:(dimlogit+1),]^2)
      loglikepen <- loglike - alpha*mawdist - sum(lambdaLasso * penbeta) - kappa * Wientropy
    } else {
      loglikepen <- loglike - alpha*mawdist - kappa * Wientropy
    }

    if (abs((loglikepen - oldloglikepen) / oldloglikepen) < MLMoption$stopratio && loop > minloop) {
      break
    }

    if (loglikepen < oldloglikepen && loop > minloop) {
      #loop
      #print(c(loglike, oldloglike, loglikepen, oldloglikepen))
      break
    }

    oldloglike <- loglike
    oldloglikepen <- loglikepen

    pj <- colSums(pij_wt)
    a <- pj / sum(pj)
    numnonzero <- sum(!(a > 0))
    if (numnonzero == numcmp) {
      cat('All components have non-positive prior\n')
      cat('pj:\n')
      print(pj)
      cat('a:\n')
      print(a)
    }

    muprev <- mu
    sigmaprev <- sigma
    betaprev <- beta

    mu <- X %*% pij_wt
    for (j in 1:numcmp) {
      if (pj[j] > 0) {
        mu[,j] <- mu[,j] / pj[j]
      } else {
        mu[,j] <- muprev[,j]
        cat('Warning: zero prior for component', j, ', Component mean, Covariance, and Beta all set to the same as previous round.\n')
      }
    }

    for (j in 1:numcmp) {
      Phi <- matrix(0, nrow = dim, ncol = dim)
      for (i in 1:numdata) {
        Phi <- Phi + pij_wt[i,j] * (X[,i] - mu[,j]) %*% t(X[,i] - mu[,j])
      }
      if (pj[j] > 0) {
        sigma[,,j] <- Phi / sum(pj[j])
      } else {
        sigma[,,j] <- sigmaprev[,,j]
      }
      sigma[,,j] <- checksingular(sigma[,,j], Xvar, 0.05)
    }

    sigma <- ConstrainSigma(a, sigma, dim, numcmp, Xvar, MLMoption)

    ### Compute barycenters by changing weights assigned to the new model
    ### and the prior model.
    if (ncompprior > 0 & alpha > 0.0) {

      lambda <- seq(0.0, 0.9, by = 0.1)

      c0$w <- a
      c0$supp <- matrix(0, dim + dim * dim, numcmp)
      c0$supp[1:dim, ] <- mu

      for (kk in 1:numcmp) {
        c0$supp[(dim + 1):(dim + dim * dim), kk] <- as.vector(sigma[, , kk])
      }

      supp <- cbind(c0$supp, cprior$supp)
      ww <- c(c0$w, cprior$w)

      options$max_iters <- 500
      options$fixpoint_niter <- 20
      options$thred_obj <- 1e-4

      objlamb <- numeric(length(lambda))
      mawdist1 <- numeric(length(lambda))
      supplamb <- matrix(0, nrow = nrow(c0$supp), ncol = ncol(c0$supp) * length(lambda))
      wwlamb <- matrix(0, nrow = nrow(c0$w), ncol = ncol(c0$w) * length(lambda))

      for (kk in 1:length(lambda)) {
        instanceW <- c(1 - lambda[kk], lambda[kk])
        c <- centroid_basic_GMM(dim, stride, instanceW, supp, ww, c0, options)

        Xmatch2 <- Xmatch[, (stride[1] + 1):(stride[1] + stride[2])]
        mawdist1[kk] <- MAWwithMatch(dim, c$supp, cprior$supp, Xmatch2)
        objlamb[kk] <- MstepObj(dim, numdata, c$supp, c$w, X, pij_wt) - alpha * mawdist1[kk]

        supplamb[, ((kk - 1) * numcmp + 1):(kk * numcmp)] <- c$supp
        wwlamb[((kk - 1) * numcmp + 1):(kk * numcmp)] <- c$w
      }

      mvind <- which.max(objlamb)
      mawdist <- mawdist1[mvind]

      mu <- supplamb[1:dim, ((mvind - 1) * numcmp + 1):(mvind * numcmp)]
      sigma <- array(supplamb[(dim + 1):(dim + dim * dim), ((mvind - 1) * numcmp + 1):(mvind * numcmp)],
                     dim = c(dim, dim, numcmp))
      a <- wwlamb[((mvind - 1) * numcmp + 1):(mvind * numcmp)]

      for (i in 1:numcmp) {
        if (a[i] < 0) a[i] <- 0
      }

      a <- a / sum(a)

      for (j in 1:numcmp) {
        sigma[, , j] <- checksingular(sigma[, , j], Xvar, 0.05)
      }

    } else {
      mawdist <- 0.0
    }


    # follow are as the same as em_MLM

    for (j in 1:numcmp) {
      if (pj[j] > 0) {
        if (pj[j] >= 3) {
          if (MLMoption$algorithm == 1) {
            # cat('the iteration in em_MLM is:',j,'\n')
            fit <- glmnet(t(Xlogit), Y, family = MLMoption$DISTR,  weights = pij_wt[,j],alpha = MLMoption$AlphaLasso, lambda = lambdaLasso[j])
            beta_nonintercept = as.matrix(fit$beta)
            beta[,j] <- rbind(fit$a0, beta_nonintercept)
          } else {
            weights <- pij_wt[, j] / pj[j]
            Xlogit_transposed <- t(Xlogit)
            # Fit the logistic regression model
            model <- glm(Y ~ Xlogit_transposed, family = binomial(link = "logit"), weights = weights)
            # Extract the coefficients (including intercept)
            beta[, j] <- coef(model)
          }
        } else {
          beta[,j] <- logisticIntercept(Y, pij_wt[,j], dimlogit)
        }
      } else {
        beta[,j] <- betaprev[,j]
      }
    }

    for (j in 1:numcmp) {
      sigmainv[,,j] <- solve(sigma[,,j])
      sigmadetsqrt[j] <- sqrt(det(sigma[,,j]))
    }

    if (kappa > 0) {
      Li <- rep(0, numdata)
      for (i in 1:numdata) {
        for (j in 1:numcmp) {
          v1 <- sum(beta[2:(dimlogit+1),j] * Xlogit[,i]) + beta[1,j]
          if (exp(v1) == Inf) {
            v2 <- -v1 + Y[i] * v1 + log(a[j]) - log(sigmadetsqrt[j]) - 0.5 * t(X[,i] - mu[,j]) %*% sigmainv[,,j] %*% (X[,i] - mu[,j])
          } else {
            v2 <- -log(1 + exp(v1)) + Y[i] * v1 + log(a[j]) - log(sigmadetsqrt[j]) - 0.5 * t(X[,i] - mu[,j]) %*% sigmainv[,,j] %*% (X[,i] - mu[,j])
          }
          if(a[j] > 0){
            Li[i] <- Li[i] + pij[i,j] * v2
          }
        }
      }
      Wi <- exp(Li/kappa)
      if (sum(Wi)==0){
        Wi <- rep(1, numdata)
      }else{
        Wi = Wi/sum(Wi)*numdata
      }
      for (i in 1:numdata){
        if (!(Wi[i] >= 0 || Wi[i] < 0)) { # NaN appears
          cat(sprintf('Wi(%d)=%f\n', i, Wi[i]))
          stop('em_MLM_maw: Numerical error with computing weights Wi')
        }
      }
    }


    loop <- loop + 1
    #print(loglikepen)

  }

  c = list()
  c$w = a;
  c$supp = matrix(0,nrow = dim+dim*dim,ncol = numcmp)
  c$supp[1:dim,] = mu;
  for (i in 1:numcmp){
    start_ind = dim+1
    end_ind = dim+dim*dim
    c$supp[start_ind:end_ind, i] <- as.vector(sigma[,,i])
  }

  result <- list(
    c = c,
    beta = beta,
    Wi = Wi,
    loglike = loglike,
    loglikepen = loglikepen
  )

  return(result)
}
