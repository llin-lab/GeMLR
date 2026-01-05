#' run CV
#' @description
#' Run Cross-Verified for GeM-LR with number of components being ncmp.
#'
#' @param k the number of folds of cross-validation
#' @param ncmp the number of clusters
#' @param nseeds the number of random seeds
#' @param rangeSeed the largest number among random seeds
#' @param vargmm the index of variables that are used in gmm
#' @param vlasso Lambda for Lasso penalty
#' @param Y the response variable
#' @param X all independent variables
#' @param Indi dummy variable
#' @param alphaLasso Elastic net parameter (default 0.8)
#' @param stopratio Convergence threshold (default 1.0e-5)
#' @param verbose Verbosity flag (default 1)
#' @param constrain Covariance constraint (default "DIAS")
#' @param diagshrink Diagonal shrinkage (default 0.9)
#' @param kappa Instance weighting (default -1)
#' @param AUC Use AUC for selection (default 1)
#' @param DISTR Distribution family (default "binomial")
#'
#' @return all AUC informations and the final classification result
#' @export
#'
runCV <- function(k=5, ncmp=c(2,3,4), nseeds=20, rangeSeed=30, 
                  vargmm, vlasso, Y, X, Indi, 
                  alphaLasso = 0.8,
                  stopratio = 1.0e-5,
                  verbose = 1,
                  constrain = "DIAS",
                  diagshrink = 0.9,
                  kappa = -1,
                  AUC = 1,
                  DISTR = "binomial") {
  library(caret)
  library(pROC)

  lcmp <- length(ncmp)
  dimgmm <- length(vargmm)
  labels <- vector("list", lcmp * k)
  guess <- vector("list", lcmp * k)

  cvAUCfinal <- matrix(0, nrow = k, ncol = lcmp)
  rownames(cvAUCfinal) <- paste(1:k, "fold", sep = " ")
  colnames(cvAUCfinal) <- paste("cluster=", ncmp, sep = "")

  bestseed <- matrix(0, nrow = k, ncol = lcmp)

  set.seed(9)
  tuningK2 <- createFolds(Y, k = k, list = TRUE)
  rseeds <- sample(1:rangeSeed, nseeds, replace = FALSE)
  dim <- ncol(X)

  for (ifold in 1:k) {
    test_index <- tuningK2[[ifold]]
    training_index <- setdiff(seq_along(Y), test_index)

    Xtraining <- scale(X[training_index, ])
    Ctest1 <- attr(Xtraining, "scaled:center")
    Stest1 <- attr(Xtraining, "scaled:scale")
    Ytraining <- Y[training_index]

    Xtt <- scale(X[test_index, ], center = Ctest1, scale = Stest1)
    Ytt <- Y[test_index]

    if (is.null(Indi)) {
      Xtrain_indi <- Xtraining
      Xtt_indi <- Xtt
    } else {
      Xtrain_indi <- cbind(Xtraining, as.data.frame(Indi)[training_index,])
      Xtt_indi <- cbind(Xtt, as.data.frame(Indi)[test_index,])
    }

    for (jj in 1:lcmp) {
      # Create MLMoption internally
      MLMoption <- init_MLMoption(
        alphaLasso = alphaLasso,
        vlasso = vlasso,
        numcmp = ncmp[jj],
        stopratio = stopratio,
        verbose = verbose,
        minloop = 3,
        maxloop = 50,
        constrain = constrain,
        diagshrink = diagshrink,
        kmseed = 0,
        algorithm = 1,
        kappa = kappa,
        AUC = AUC,
        DISTR = DISTR,
        NOEM = 0,
        Yalpha = 1.0
      )
      MLMoption$lambdaLasso <- rep(vlasso, ncmp[jj])

      est <- estimateBestSD(Xtraining[, vargmm], Xtrain_indi, Ytraining, MLMoption, rseeds)
      c <- est[[1]]
      beta <- est[[2]]
      bestseed[ifold, jj] <- est[[8]]

      MLMoption$kmseed <- bestseed[ifold, jj]
      gmm <- GMMFormatConvert(dimgmm, c)
      a2 <- gmm[[2]]
      mu2 <- gmm[[3]]
      sigma2 <- gmm[[4]]

      mlm <- MLMclassify(a2, mu2, sigma2, beta, Xtt[, vargmm], Xtt_indi)
      pyi <- mlm[[1]]

      labels[[jj + (ifold - 1) * lcmp]] <- Ytt
      guess[[jj + (ifold - 1) * lcmp]] <- pyi

      roc_obj <- roc(Ytt, pyi, levels = c(0,1))
      cvAUCfinal[ifold, jj] <- roc_obj$auc
    }
  }

  return(list(cvAUCfinal = cvAUCfinal, labels = labels, guess = guess))
}
