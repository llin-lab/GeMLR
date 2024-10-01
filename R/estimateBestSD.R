#' estimate the model that have the max AUC
#' @description
#' Repeat EM algorithm for estimating MLM using multiple random seeds for kmeans initialization. Choose the best seed based on training accuracy or AUC.
#'
#' @param X all non-dummy feature variables
#' @param Xlogit all feature variables
#' @param Y the response variable
#' @param MLMoption all necessary variables
#' @param seedlist the list of random seeds
#'
#' @return all informations about the gmm and logistic model
#' @export
estimateBestSD <- function(X, Xlogit, Y, MLMoption, seedlist) {

  library(pROC)
  if (any(seedlist < 0)) {
    cat('Warning: seedlist contains negative value, NOT recommended, can confuse the meaning of output bestseed\n')
  }

  nseedkm <- length(seedlist)
  numdata <- nrow(X)
  dim <- ncol(X)
  dimlogit <- ncol(Xlogit)
  InitCluster <- MLMoption$InitCluster
  numinit <- ncol(InitCluster)
  nseed <- nseedkm + numinit

  if (nseed < 1) {
    cat('Error: number of initializations is', nseed, '\n')
    cat('Input correct seedlist or MLMoption.InitCluster\n')
    return(NULL)
  }

  accuracy <- numeric(nseed)
  AUC <- numeric(nseed)

  for (ii in 1:nseed) {
    if (ii <= nseedkm) {
      MLMoption$InitCluster <- matrix(0, 0, 0)
      MLMoption$kmseed <- seedlist[ii]
    } else {
      MLMoption$InitCluster <- InitCluster[, ii - nseedkm]
      MLMoption$kmseed <- 0
    }
    result <- estimate(X, Xlogit, Y, MLMoption)
    c <- result$c
    beta <- result$beta
    rsigma2 <- result$rsigma2
    loglike <- result$loglike
    loglikepen <- result$loglikepen
    loglikeY <- result$loglikeY
    Wi <- result$Wi

    format_result <- GMMFormatConvert(dim, c)
    numcmp2 <- format_result$numcmp
    a2 <- format_result$a
    mu2 <- format_result$mu
    sigma2 <- format_result$sigma

    if (MLMoption$DISTR == 'binomial') {
      classify_result <- MLMclassify(a2, mu2, sigma2, beta, X, Xlogit)
      pyi <- classify_result$pyi
      pij <- classify_result$pij
      pyihard <- pyi > 0.5
      accuracy[ii] <- sum(Y == pyihard) / length(Y)

      auc_result <- roc(Y, pyi, levels = c(0,1))
      AUC[ii] <- auc_result$auc
    }else{
      regress_result <- MLMregress(a2, mu2, sigma2, beta, X, Xlogit)
      accuracy[ii] <- sum((pyi - t(Y))^2) / length(Y)  # MSE for regression
    }
  }

  if (nseed == 1) {
    if (nseedkm > 0) {
      bestseed <- seedlist[1]
    } else {
      bestseed <- -1
    }
    return(list(c = c, beta = beta, rsigma2 = rsigma2, loglike = loglike, loglikepen = loglikepen, loglikeY = loglikeY, Wi = Wi, bestseed = bestseed))
  }

  if (MLMoption$AUC == 1 && MLMoption$DISTR == 'binomial') {
    bestID <- which.max(AUC)
  } else {
    if (MLMoption$DISTR == 'binomial') {
      bestID <- which.max(accuracy)
    }else{
      bestID = which.min(accuracy)
    }
  }

  if (bestID <= nseedkm) {
    bestseed <- seedlist[bestID]
    MLMoption$kmseed <- bestseed
    MLMoption$InitCluster <- matrix(0, 0, 0)
    result <- estimate(X, Xlogit, Y, MLMoption)
  } else {
    bestseed <- bestID - nseedkm
    MLMoption$InitCluster <- InitCluster[, bestseed]
    MLMoption$kmseed <- 0
    result <- estimate(X, Xlogit, Y, MLMoption)
    bestseed <- -bestseed
  }

  MLMoption$InitCluster <- InitCluster

  return(list(c = result$c, beta = result$beta, rsigma2 = result$rsigma2, loglike = result$loglike, loglikepen = result$loglikepen, loglikeY = result$loglikeY, Wi = result$Wi, bestseed = bestseed))
}
