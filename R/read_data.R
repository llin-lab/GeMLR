#' Read dataset and prepare X/Xs/Y/Indi (with automatic column name fixing)
#'
#' @param dat_path Path to data file.
#' @param ycol     Optional Y column (index or name). If NULL, use last column.
#' @param Indi_col Optional indicator column (single index or single name).
#' @param encoding Optional file encoding for text files.
#' @return list(dim, numdata, rawdat, X, Xs, Y, Indi)
#' @export
read_data <- function(dat_path, ycol = NULL, Indi_col = NULL, encoding = NULL) {
  stopifnot(length(dat_path) == 1, file.exists(dat_path))
  ext <- tolower(tools::file_ext(dat_path))
  
  guess_sep <- function(path, enc = encoding) {
    con <- if (is.null(enc)) file(path, "r") else file(path, "r", encoding = enc)
    on.exit(close(con), add = TRUE)
    ln <- readLines(con, n = 1)
    if (length(ln) == 0) return("")
    if (grepl(",", ln)) return(",")
    if (grepl("\t", ln)) return("\t")
    ""
  }
  
  # 1) read
  if (ext %in% c("csv")) {
    rawdat <- if (is.null(encoding)) {
      utils::read.csv(dat_path, header = TRUE, stringsAsFactors = FALSE)
    } else {
      utils::read.csv(dat_path, header = TRUE, stringsAsFactors = FALSE, fileEncoding = encoding)
    }
  } else if (ext %in% c("tsv", "tab")) {
    rawdat <- if (is.null(encoding)) {
      utils::read.delim(dat_path, header = TRUE, stringsAsFactors = FALSE)
    } else {
      utils::read.delim(dat_path, header = TRUE, stringsAsFactors = FALSE, fileEncoding = encoding)
    }
  } else if (ext %in% c("txt", "dat", "")) {
    sep_guess <- guess_sep(dat_path, encoding)
    rawdat <- if (is.null(encoding)) {
      utils::read.table(dat_path, header = TRUE, sep = sep_guess, stringsAsFactors = FALSE)
    } else {
      utils::read.table(dat_path, header = TRUE, sep = sep_guess, stringsAsFactors = FALSE, fileEncoding = encoding)
    }
  } else if (ext %in% c("rds")) {
    rawdat <- as.data.frame(readRDS(dat_path), stringsAsFactors = FALSE)
  } else if (ext %in% c("xlsx", "xls")) {
    if (!requireNamespace("readxl", quietly = TRUE)) stop("readxl is required for Excel files.")
    rawdat <- as.data.frame(readxl::read_excel(dat_path), stringsAsFactors = FALSE)
  } else if (ext %in% c("sav", "sas7bdat", "dta")) {
    if (!requireNamespace("haven", quietly = TRUE)) stop("haven is required for SPSS/SAS/Stata files.")
    if (ext == "sav")      rawdat <- as.data.frame(haven::read_sav(dat_path), stringsAsFactors = FALSE)
    if (ext == "sas7bdat") rawdat <- as.data.frame(haven::read_sas(dat_path), stringsAsFactors = FALSE)
    if (ext == "dta")      rawdat <- as.data.frame(haven::read_dta(dat_path), stringsAsFactors = FALSE)
  } else {
    rawdat <- if (is.null(encoding)) {
      utils::read.table(dat_path, header = TRUE, stringsAsFactors = FALSE)
    } else {
      utils::read.table(dat_path, header = TRUE, stringsAsFactors = FALSE, fileEncoding = encoding)
    }
  }
  if (!is.data.frame(rawdat) || ncol(rawdat) < 2L) stop("Invalid data frame.")
  
  # 2) Y (binary {0,1})
  if (is.null(ycol)) {
    ycol <- ncol(rawdat)
  } else if (is.character(ycol)) {
    stopifnot(length(ycol) == 1, ycol %in% colnames(rawdat))
    ycol <- match(ycol, colnames(rawdat))
  } else {
    stopifnot(is.numeric(ycol), length(ycol) == 1, ycol >= 1, ycol <= ncol(rawdat))
  }
  Y <- rawdat[[ycol]]
  if (!all(is.finite(Y)) || !all(Y %in% c(0, 1))) stop("Y must be strict {0,1}.")
  Y <- as.numeric(Y)
  
  # 3) Indi (single column if provided; else auto-detect all strict 0/1 excluding Y)
  n_all <- ncol(rawdat)
  other_cols <- setdiff(seq_len(n_all), ycol)
  is_bin <- function(v) is.numeric(v) && all(is.finite(v)) && all(v %in% c(0, 1))
  
  if (!is.null(Indi_col)) {
    if (is.character(Indi_col)) {
      stopifnot(length(Indi_col) == 1, Indi_col %in% colnames(rawdat))
      indi_col <- match(Indi_col, colnames(rawdat))
    } else {
      stopifnot(is.numeric(Indi_col), length(Indi_col) == 1, Indi_col >= 1, Indi_col <= n_all)
      indi_col <- as.integer(Indi_col)
    }
    if (indi_col == ycol) stop("Indi cannot be the same as Y.")
    if (!is_bin(rawdat[[indi_col]])) stop("Specified Indi column is not strict {0,1}.")
    Indi <- as.matrix(rawdat[, indi_col, drop = FALSE])
    colnames(Indi) <- colnames(rawdat)[indi_col]
  } else {
    cand <- other_cols[vapply(rawdat[other_cols], is_bin, logical(1))]
    Indi <- if (length(cand) == 0L) NULL else as.matrix(rawdat[, cand, drop = FALSE])
  }
  
  # 4) X / Xs
  X_cols <- setdiff(seq_len(n_all), c(ycol, if (is.null(Indi)) integer(0) else match(colnames(Indi), colnames(rawdat))))
  if (length(X_cols) == 0L) stop("X is empty after removing Y and Indi.")
  X <- as.data.frame(rawdat[, X_cols, drop = FALSE])
  
  is_binary_vec <- function(v) is.numeric(v) && all(is.finite(v)) && all(sort(unique(v)) %in% c(0, 1))
  Xs <- X
  for (nm in colnames(X)) {
    v <- X[[nm]]
    if (!is_binary_vec(v)) Xs[[nm]] <- scale(v)[, 1]
  }
  
  dim <- ncol(X)
  numdata <- nrow(rawdat)
  
  return(list(dim = dim,
              numdata = numdata,
              rawdat = rawdat,
              X = X,
              Xs = Xs,
              Y = Y,
              Indi = Indi))
}