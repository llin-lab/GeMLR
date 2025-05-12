#' @title read dataset
#' @description read the dataset and relevant parameters.
#' @param dat_road the road of your .txt dataset file
#' @param num_gmm the number of variables that are used in the gmm model, default=5
#' @param alphaLasso the alpha number that are used in the elastic net regression
#'
#' @return a list of dataset and necessary variables
#' @export
read_data <- function(dat_road,sep_mark=' ', num_gmm=NULL,alphaLasso=0.8, ycol=NULL, Indi_col=1, gmm_var=NULL){
  library(dplyr)
  library(glmnet)

  rawdat = read.table(dat_road,sep = sep_mark)
  if (!all(rawdat[[Indi_col]] %in% c(0, 1))) {
    stop("The Indi column in your data must be 0 or 1(Indi)")
  }


  if (is.null(ycol)) {
    ycol <- ncol(rawdat)
    cat("Using the last column of the dataframe as ycol.\n")
  } else {
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
    if (!all(X[[col]] %in% c(0, 1)) || length(unique(X[[col]])) != 2) {
      Xs[[col]] <- scale(X[[col]])[, 1]
    }
  }
  X_var <- sapply(X, var)

  vargmm = numeric(0)

  if (is.null(num_gmm)) {
    # default: use all
    vargmm <- order(X_var, decreasing = TRUE)
  } else {
    # user-defined: Loop through each item in the user input vector
    if (num_gmm==0){
      if (is.null(gmm_var)){
        "Please provide the column names or column indexes of the variables you want to use for the GMM model!"
      } else {
        for (item in gmm_var) {
          if (item %in% colnames(rawdat)) {
            vargmm <- c(vargmm, which(colnames(rawdat) == item))
          } else if (is.numeric(as.numeric(item)) && item %in% 1:dim) {
            vargmm <- c(vargmm, item)
          } else {
            warning(paste("Invalid input:", item))
          }
        }
      }
    } else if (is.numeric(num_gmm) && num_gmm > 0 && floor(num_gmm) == num_gmm){ # use variance, but with user-defined number
      if (num_gmm<=dim){
        vargmm <- order(X_var, decreasing = TRUE)[1:num_gmm]
      } else {
        vargmm <- order(X_var, decreasing = TRUE)[1:dim]
        print("The number of columns you input is too large. By default, all variables are selected to participate in the GMM model.")
      }
    } else {
      print('Invalid num_gmm!')
    }
  }


  #X <- cbind(X, Indi)
  Y <- as.vector(Y)
  cv_fit <- suppressWarnings(cv.glmnet(data.matrix(X), Y, alpha = 1, family = "binomial", nfolds = 5))
  B <- coef(cv_fit, s = "lambda.min")
  vlasso <- cv_fit$lambda.min

  return(list(dim=dim,numdata=numdata,rawdat=rawdat,vargmm=vargmm,vlasso=vlasso,X=X,Xs=Xs,Y=Y,Indi=Indi))

}
