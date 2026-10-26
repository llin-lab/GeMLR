#' VS
#' @description
#' variable selection for Gaussian mixtured model.
#'
#' @param comp use vs on which components
#' @param a2 the weight of each cluster
#' @param mu2 the mean of each cluster
#' @param sigma2 the covariance matrix of each cluster
#' @param dimgmm dimensions of gmm
#'
#' @return a list of information about the explainable variables
#' @export


VS <- function(comp, a2, mu2, sigma2, dimgmm = nrow(mu2), vargmm = NULL, varnames = NULL) {
  if (dimgmm == 1) {
    # only one GMM variable: nothing to select (the MATLAB code returns an empty set here)
    varsel <- 1
    wm <- accuracy(a2, mu2, sigma2, matrix(1, 1, 1))$Acc[comp, 1]
    remain <- integer(0)
  } else {
    index <- creatindex(dimgmm, dimgmm - 1)
    Acc <- accuracy(a2, mu2, sigma2, index)$Acc

    Accfull <- Acc[, 1]
    Acc <- Acc[, -1, drop = FALSE]
    wm <- max(Acc[comp, ])
    vm <- which.max(Acc[comp, ])
    varsel <- index[vm + 1, 1]
    remain <- setdiff(1:dimgmm, varsel)
    Accref <- 0.01

    while ((wm - Accref) / Accref > 0.1 && abs(wm - Accfull[comp]) > 0.01) {
      Accref <- wm
      index <- cbind(matrix(varsel, nrow = length(remain), ncol = length(varsel), byrow = TRUE), remain)
      Acc <- accuracy(a2, mu2, sigma2, index)$Acc
      wm <- max(Acc[comp, ])
      vm <- which.max(Acc[comp, ])
      varsel <- index[vm, ]
      remain <- setdiff(1:dimgmm, varsel)
    }
  }

  pos <- as.integer(varsel)
  col <- if (is.null(vargmm)) NULL else vargmm[pos]
  if (is.numeric(col)) col <- as.integer(col)
  nm <- if (is.character(col)) {
    col
  } else if (!is.null(varnames) && !is.null(col)) {
    varnames[col]
  } else if (!is.null(rownames(mu2))) {
    rownames(mu2)[pos]
  } else {
    NULL
  }

  return(list(varsel = pos, varsel_col = col, varsel_name = nm, wm = wm))
}
