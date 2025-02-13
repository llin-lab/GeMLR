#' @title read dataset
#' @description read the dataset and relevant parameters.
#' @param dat_road the road of your .txt dataset file
#' @param num_gmm the number of variables that are used in the gmm model, default=5
#' @param alphaLasso the alpha number that are used in the elastic net regression
#'
#' @return a list of dataset and necessary variables
#' @export
read_data <- function(dat_road,num_gmm=NULL,alphaLasso=0.8, ycol=NULL, Indi_col=1){
  #library(reticulate)
  #dat_road =  "D:\\暑研\\DukeCHSI-GeM-LR-4874c53\\Data\\simulation_data.csv"
  library(dplyr)
  library(glmnet)

  rawdat = read.table(dat_road,sep = ',')
  if (!all(rawdat[[Indi_col]] %in% c(0, 1))) {
    stop("The Indi column in your data must be 0 or 1(Indi)")
  }


  if (is.null(ycol)) {
    # 如果用户没有输入ycol，使用ncol(df)
    ycol <- ncol(rawdat)
    cat("Using the last column of the dataframe as ycol.\n")
  } else {
    # 如果用户输入了ycol，使用用户指定的值
    cat("Using user-specified ycol value.\n")
  }

  if (!all(rawdat[[ycol]] %in% c(0, 1))) {
    stop("The Y column in your data must be 0 or 1(Response)")
  }

  numdata = nrow(rawdat)
  dim = ncol(rawdat)
  Y = rawdat[,ycol]
  X = rawdat[,-c(ycol,Indi_col)]
  Indi = rawdat[,Indi_col]
  dim = ncol(X) # dim是X的维度
  # 对X的每一列进行处理
  Xs <- X
  for (col in names(X)) {
    # 检查当前列是否是二元变量
    if (!all(X[[col]] %in% c(0, 1)) || length(unique(X[[col]])) != 2) {
      # 如果不是二元变量，则进行标准化
      Xs[[col]] <- scale(X[[col]])[, 1]
    }
    # 如果是二元变量，保持原值不变
  }
  X_var <- sapply(X, var)
  # 检查num_gmm是否被指定
  if (is.null(num_gmm)) {
    # 如果num_gmm未指定，vargmm赋值为X_var的所有列号
    vargmm <- order(X_var, decreasing = TRUE)
  } else {
    # 如果num_gmm已指定，vargmm赋值为1:num_gmm的列号
    vargmm <- order(X_var, decreasing = TRUE)[1:num_gmm]
  }


  X <- cbind(X, Indi)
  Y <- as.vector(Y)
  cv_fit <- suppressWarnings(cv.glmnet(data.matrix(X), Y, alpha = 1, family = "binomial", nfolds = 5))
  B <- coef(cv_fit, s = "lambda.min")
  vlasso <- cv_fit$lambda.min

  return(list(dim=dim,numdata=numdata,rawdat=rawdat,vargmm=vargmm,vlasso=vlasso,X=X,Xs=Xs,Y=Y,Indi=Indi))

}
