#' createindex
#' @description
#' generate an index for a specific number of tag subsets.
#'
#' @return the index result

creatindex <- function(dimofx, m) {
  index <- NULL
  if (m == 0) {
    for (mm in 1:(dimofx - 1)) {
      index1 <- t(combn(1:dimofx, (dimofx - mm)))
      index <- rbind(index, index1)
    }
  } else {
    index <- t(combn(1:dimofx, (dimofx - m)))
  }
  index <- rbind(1:dimofx, index)
  return(index)
}

