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


VS <- function(comp, a2, mu2, sigma2, dimgmm) {
  index <- creatindex(dimgmm, dimgmm - 1)             # requires the item-3 fix
  Acc <- accuracy(a2, mu2, sigma2, index)$Acc         # A: take the matrix
  Accfull <- Acc[, 1]
  Acc <- Acc[, -1, drop = FALSE]                      # keep k x d even if d = 1
  vm <- which.max(Acc[comp, ]); wm <- Acc[comp, vm]
  varsel <- index[vm + 1, 1]
  remain <- setdiff(1:dimgmm, varsel)
  Accref <- 0.01
  while ((wm - Accref) / Accref > 0.1 && abs(wm - Accfull[comp]) > 0.01) {
    Accref <- wm
    index <- cbind(matrix(varsel, nrow = length(remain),          # B: rows [varsel, r]
                          ncol = length(varsel), byrow = TRUE), remain)
    Acc <- accuracy(a2, mu2, sigma2, index)$Acc
    vm <- which.max(Acc[comp, ]); wm <- Acc[comp, vm]
    varsel <- index[vm, ]
    remain <- setdiff(1:dimgmm, varsel)
  }
  list(varsel = varsel, wm = wm)
}
