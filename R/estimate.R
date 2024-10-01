#' estimate
#' @description
#' use EM algorithm for estimating MLM.
#'
#' @param X all non-dummy feature variables
#' @param Xlogit all feature variables
#' @param Y the response variable
#' @param MLMoption all necessary variables
#'
#' @return all informations about the gmm and logistic model
#' @export
#'
estimate <- function(X, Xlogit, Y, MLMoption) {
  library(pROC)
  dim <- ncol(X)
  numdata <- nrow(X)
  numcmp <- MLMoption$numcmp[1]

  init_params <- initemMLM(X, Xlogit, Y, MLMoption)
  muinit <- init_params$muinit
  sigmainit <- init_params$sigmainit
  ainit <- init_params$ainit
  betainit <- init_params$beta
  rsigma2 <- init_params$rsigma2

  if (MLMoption$verbose > 0) {
    if (MLMoption$DISTR == 'binomial') {
      classify_result <- MLMclassify(ainit, muinit, sigmainit, betainit, X, Xlogit)
      pyi <- classify_result$pyi
      pij <- classify_result$pij
      pyihard <- pyi > 0.5
      cat('Classification accuracy by initial model (within training):', sum(Y == pyihard) / length(Y), '\n')
      auc_result <- roc(Y, pyi)
      AUC <- auc_result$auc
      cat('Classification AUC by initial model (within training):', AUC, '\n')
    }else{
      regress_result <- MLMregress(ainit,muinit,sigmainit,betainit,X,Xlogit)
      v1 <- sum((pyi - t(Y))^2) / numdata
      cat(sprintf('Regression MSE by initial model (within training): %f\n', v1))
    }
  }

  cinit <- list()
  cinit$w <- ainit
  cinit$supp <- matrix(0, nrow = dim + dim * dim, ncol = numcmp)
  cinit$supp[1:dim, ] <- muinit
  for (k3 in 1:numcmp) {
    cinit$supp[(dim + 1):(dim + dim * dim), k3] <- as.vector(sigmainit[, , k3])
  }

  if (MLMoption$NOEM == 1) {
    c <- cinit
    beta <- betainit
    loglike <- 0
    loglikepen <- 0.0
    Wi <- rep(1, numdata)
  } else {
    if (MLMoption$DISTR == 'binomial') {
      em_result <- em_MLM(X, Xlogit, Y, cinit, betainit, MLMoption)
      c <- em_result$c
      beta <- em_result$beta
      Wi <- em_result$Wi
      loglike <- em_result$loglike
      loglikepen <- em_result$loglikepen
    }else{
      em_result <- em_MLM_regress(X, Xlogit, Y, cinit, betainit, rsigma2 , MLMoption)
      c <- em_result$c
      beta <- em_result$beta
      Wi <- em_result$Wi
      rsigma2 <- em_result$rsigma2
      loglike <- em_result$loglike
      loglikepen <- em_result$loglikepen
    }
  }

  format_result <- GMMFormatConvert(dim, c)
  numcmp2 <- format_result$numcmp
  a2 <- format_result$a
  mu2 <- format_result$mu
  sigma2 <- format_result$sigma

  if (MLMoption$DISTR == 'binomial') {
    classify_result <- MLMclassify(a2, mu2, sigma2, beta, X, Xlogit)
    pyi <- classify_result$pyi
    pij <- classify_result$pij
    w2 <- c(sum(Y == 1) / length(Y), sum(Y == 0) / length(Y))
    loglikeY <- LoglikehoodYBalance(pyi, Y, w2)
    if (MLMoption$verbose > 0) {
      pyihard <- pyi > 0.5
      cat('Classification accuracy after EM (within training):', sum(Y == pyihard) / length(Y), '\n')
      auc_result <- roc(Y, pyi) # default: 0 represents control group
      AUC <- auc_result$auc
      cat('Classification AUC after EM (within training):', AUC, '\n')
    }
  } else {
    regress_result <- MLMregress(a2, mu2, sigma2, beta, X, Xlogit)
    pyi <- regress_result$pyi
    pij <- regress_result$pij
    v1 <- sum((pyi - Y)^2) / length(Y)
    loglikeY <- v1
    if (MLMoption$verbose > 0) {
      cat('Regression MSE after EM (within training):', v1, '\n')
    }
  }

  cat('Seed:', MLMoption$kmseed, 'Done---------------------\n')

  return(list(c = c, beta = beta, rsigma2 = rsigma2, loglike = loglike, loglikepen = loglikepen, loglikeY = loglikeY, Wi = Wi))
}
