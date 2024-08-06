#' run CV
#' @description
#' Run Cross-Verified for GeM-LR with number of components being ncmp.
#'
#' @param kkk the number of folds of cross-validation
#' @param ncmp the number of clusters
#' @param nseeds the number of random seeds
#' @param rangeSeed the largest number among random seeds
#' @param vargmm the number of variables that are used in gmm
#' @param Y1 the response variable
#' @param X1 all non-dummy variables
#' @param Indi dummy variable
#' @param MLMoption all necessary variables
#'
#' @return all AUC informations and the final classification result
#' @export
#'
runCV <- function(kkk, ncmp, nseeds, rangeSeed, vargmm, Y1, X1, Indi, MLMoption) {
  library(caret)
  library(pROC)

  lcmp <- length(ncmp)
  dimgmm <- length(vargmm)
  labels <- vector("list", lcmp * kkk)
  guess <- vector("list", lcmp * kkk)
  cvAUCfinal <- matrix(0, nrow = kkk, ncol = lcmp)
  bestseed <- matrix(0, nrow = kkk, ncol = lcmp)

  set.seed(9)
  tuningK2 <- createFolds(Y1, k = kkk, list = TRUE)
  rseeds <- sample(1:rangeSeed, nseeds, replace = FALSE)
  dim <- ncol(X1)

  for (ifold in 1:kkk) {
    test_index <- tuningK2[[ifold]]
    training_index <- setdiff(seq_along(Y1), test_index)

    Xtraining <- scale(X1[training_index, ])
    Ctest1 <- attr(Xtraining, "scaled:center")
    Stest1 <- attr(Xtraining, "scaled:scale")
    Ytraining <- Y1[training_index]

    Xtt <- scale(X1[test_index, ], center = Ctest1, scale = Stest1)
    Ytt <- Y1[test_index]

    if (is.null(Indi)) {
      Xtrain_indi <- Xtraining
      Xtt_indi <- Xtt
    } else {
      Xtrain_indi <- cbind(Xtraining, Indi[training_index])
      Xtt_indi <- cbind(Xtt, Indi[test_index])
    }

    for (jj in 1:lcmp) {
      MLMoption$numcmp <- ncmp[jj]
      MLMoption$lambdaLasso <- rep(MLMoption$lambdaLasso[1], ncmp[jj])

      # cat('the runCV iterations is:',jj,'\n')
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

      roc_obj <- roc(Ytt, pyi,levels = c(0,1))
      cvAUCfinal[ifold, jj] <- roc_obj$auc
    }
  }

  return(list(cvAUCfinal = cvAUCfinal, labels = labels, guess = guess))
}
