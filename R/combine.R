#' combine
#'
#'
#' @return the combined result
#'

combine <- function(...) {
  args <- list(...)
  numInputs <- length(args)

  numRows <- sapply(args, nrow)
  numCols <- sapply(args, ncol)

  maxRows <- sum(numRows, na.rm = TRUE)
  maxCols <- max(numCols, na.rm = TRUE)

  buffer <- matrix(0, nrow = maxRows, ncol = maxCols)

  currRow <- 1
  for (count in seq_along(args)) {
    if (is.character(args[[count]])) {
      args[[count]] <- utf8ToInt(args[[count]])
    }
    rowIndex <- currRow:(currRow + numRows[count] - 1)
    colIndex <- 1:numCols[count]
    buffer[rowIndex, colIndex] <- args[[count]]
    currRow <- currRow + numRows[count]
  }

  return(buffer)
}
