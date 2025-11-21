#' Compute vargmm and vlasso for GeMLR Models
#'
#' @description
#' Replicates the parameter computation logic from the original read_data() function.
#'
#' @param X Raw (unstandardized) predictor matrix (n x p).
#' @param Y Binary response vector (0/1).
#' @param Indi Indicator variables. Should match the Indi used in modeling.
#' @param num_gmm Number of top-variance features (default 5).
#'                If NULL, uses all features.
#'                If 0, must provide gmm_var.
#' @param gmm_var Manual feature specification (only when num_gmm = 0).
#'
#' @return A list with vargmm (indices) and vlasso (numeric).
#' @export

compute_gemlr_params <- function(X, Y, Indi = NULL, num_gmm = 5, gmm_var = NULL) {
  library(glmnet)
  
  dim <- ncol(X)
  X_var <- apply(X, 2, var)
  vargmm <- numeric(0)
  
  if (is.null(num_gmm)) {
    vargmm <- order(X_var, decreasing = TRUE)
  } else {
    if (num_gmm == 0) {
      if (is.null(gmm_var)) {
        stop("Please provide the column names or column indexes of the variables you want to use for the GMM model!")
      } else {
        for (item in gmm_var) {
          if (item %in% colnames(X)) {
            vargmm <- c(vargmm, which(colnames(X) == item))
          } else if (is.numeric(as.numeric(item)) && item %in% 1:dim) {
            vargmm <- c(vargmm, item)
          } else {
            warning(paste("Invalid input:", item))
          }
        }
      }
    } else if (is.numeric(num_gmm) && num_gmm > 0 && floor(num_gmm) == num_gmm) {
      if (num_gmm <= dim) {
        vargmm <- order(X_var, decreasing = TRUE)[1:num_gmm]
      } else {
        vargmm <- order(X_var, decreasing = TRUE)[1:dim]
        print("The number of columns you input is too large. By default, all variables are selected to participate in the GMM model.")
      }
    } else {
      print('Invalid num_gmm!')
    }
  }
  
  Y <- as.vector(Y)
  set.seed(1)
  cv_fit <- suppressWarnings(
    cv.glmnet(data.matrix(cbind(X, Indi)), Y, alpha = 1, family = "binomial", nfolds = 5)
  )
  vlasso <- cv_fit$lambda.min
  
  return(list(vargmm = vargmm, vlasso = vlasso))
}