#' createindex
#' @description
#' generate an index for a specific number of tag subsets.
#'
#' @return the index result

creatindex <- function(dimofx, m) {
  index <- NULL
  if (m == 0) {
    for (mm in 1:(dimofx - 1)) {
      index1 <- t(combn(1:dimofx, dimofx - mm))
      index <- if (is.null(index)) index1 else combine(index, index1)
    }
  } else {
    index <- t(combn(1:dimofx, dimofx - m))
  }
  index <- combine(matrix(1:dimofx, nrow = 1), index)
  index
}
