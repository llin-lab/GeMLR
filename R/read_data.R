#' @title read dataset
#' @description read the dataset and relevant parameters.
#' @param dat_road the road of your .txt dataset file
#' @param num_gmm the number of variables that are used in the gmm model, default=5
#' @param ncmp a vector that store the different number of clusters
#' @param alphaLasso the alpha number that are used in the elastic net regression
#' @param nseeds the number of random number seeds
#' @param rangeSeed the largest number in seeds' list
#' @param kkk the number of folds of cross-validation
#'
#' @return a list of dataset and necessary variables
#' @export
read_data <- function(dat_road,num_gmm = 5,ncmp=c(2,3,4),alphaLasso=0.8,nseeds=20,rangeSeed=30,kkk=5){
  #library(reticulate)
  library(dplyr)
  library(glmnet)

  rawdat = read.table(dat_road)
  if (!all(rawdat[[1]] %in% c(0, 1))) {
    stop("The first column in your data must be 0 or 1(Indi)")
  }

  if (!all(rawdat[[ncol(rawdat)]] %in% c(0, 1))) {
    stop("The final column in your data must be 0 or 1(Response)")
  }

  numdata = nrow(rawdat)
  dim = ncol(rawdat)
  Y1 = rawdat[,dim]
  X1 = rawdat[,2:(dim-1)]
  Indi = rawdat[,1]
  dim = ncol(X1)
  X1s <- X1 %>%
    mutate_all(~ scale(.)[, 1])
  X1_var <- sapply(X1, var)
  vargmm <- order(X1_var, decreasing = TRUE)[1:num_gmm]

  X <- cbind(X1s, Indi)
  Y1 <- as.vector(Y1)
  cv_fit <- cv.glmnet(data.matrix(X), Y1, alpha = alphaLasso, family = "binomial", nfolds = 5)
  B <- coef(cv_fit, s = "lambda.min")
  vlasso <- cv_fit$lambda.min
  numcmp = ncmp[1]

  MLMoption <- list()
  MLMoption$AlphaLasso=alphaLasso;
  MLMoption$lambdaLasso=vlasso*rep(1, times = numcmp);
  MLMoption$stopratio=1.0e-5;
  MLMoption$kappa=-1.0;
  MLMoption$verbose=1;
  MLMoption$minloop=3;
  MLMoption$maxloop=50;
  MLMoption$constrain='DIAS';
  MLMoption$diagshrink=0.9;
  MLMoption$kmseed=0;
  MLMoption$algorithm=1;
  MLMoption$numcmp=numcmp;
  MLMoption$AUC=1;
  MLMoption$DISTR='binomial';
  MLMoption$NOEM=0;
  MLMoption$InitCluster = data.frame();
  MLMoption$Yalpha = 1.0;

  nseeds=nseeds;
  rangeSeed=rangeSeed;

  kkk=kkk;

  return(list(dim=dim,MLMoption=MLMoption,ncmp=ncmp,nseeds=nseeds,numcmp=numcmp,numdata=numdata,rangeSeed=rangeSeed,rawdat=rawdat,vargmm=vargmm,vlasso=vlasso,X1=X1,X1s=X1s,Y1=Y1,kkk=kkk,Indi=Indi))

}
