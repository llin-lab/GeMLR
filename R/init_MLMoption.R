
#' @title init MLMoption parameters
#' @description Specify the GeM-LR model via MLMoption
#'
#'
#' @return a list of necessary parameters in a list named 'MLMoption'
#' @export
init_MLMoption <- function(alphaLasso = 0.5, vlasso = 1, numcmp = 2, stopratio = 1.0e-5,
                           verbose = 1, minloop = 3, maxloop = 50, constrain = 'DIAS',
                          diagshrink = 0.9, kmseed = 0, algorithm = 1, kappa = -1,
                          AUC = 1, DISTR = 'binomial', NOEM = 0,Yalpha = 1.0) {
  MLMoption <- list()
  MLMoption$AlphaLasso <- alphaLasso # if set to 1, this is Lasso, value smaller than 1 is elastic net
  MLMoption$lambdaLasso <- vlasso * rep(1, times = numcmp)
  MLMoption$stopratio <- stopratio # Threshold controling the number of EM iterations in em_MLM
  MLMoption$kappa <- kappa # Weights on instances can be part of the optimization
                           # if the option is evoked.
                           # For the paper, we set negative value that disables the option of weighted instances.
  MLMoption$verbose <- verbose # T
  if (minloop < 2){
    stop("minloop must be at least 2.")
  }
  MLMoption$minloop <- minloop # must be at least 2, otherwise the program automatically uses 2

  if (maxloop<minloop){
    stop("maxloop must be greater than minloop.")
  }
  MLMoption$maxloop <- maxloop # maximum number of iterations in EM

  MLMoption$constrain <- constrain
  MLMoption$diagshrink <- diagshrink  # larger value indicates more shrinkage towards diagonal, only used if the constraint is 'DIAS'
  MLMoption$kmseed <- kmseed
  MLMoption$algorithm <- algorithm  # 1 for Lasso, 0 for Logistic without variable selection
  MLMoption$numcmp <- numcmp
  MLMoption$AUC <- AUC # If set AUC=1, use AUC to pick the best seed in estimateBestSD, otherwise, use accuracy
  MLMoption$DISTR <- DISTR # 'binomial' for classification, 'normal' for regression
  MLMoption$NOEM <- NOEM # if equal 1, then only run initialization, NO EM update
  MLMoption$InitCluster <- data.frame() # If MLMoption.InitCluster is set as a vector of length=data size, it is
                                        # assumed to be the initial cluster labels used by initemMLM.
  MLMoption$Yalpha <- Yalpha

  return(MLMoption)
}
