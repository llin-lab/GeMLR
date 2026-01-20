#' Read dataset and prepare X/Xs/Y/Indi (with auto header detection and RData support)
#'
#' @param dat_path Path to data file.
#' @param ycol     Optional Y column (index or name). If NULL, use last column.
#' @param Indi_col Optional indicator column (single index or single name).
#' @param encoding Optional file encoding for text files.
#' @return list(dim, numdata, rawdat, X, Xs, Y, Indi)
#' @export
read_data <- function(dat_path, ycol = NULL, Indi_col = NULL, encoding = NULL) {
  if (is.data.frame(dat_path) || is.matrix(dat_path)) {
    rawdat <- as.data.frame(dat_path, stringsAsFactors = FALSE)
    has_header <- TRUE
  } else {
    stopifnot(length(dat_path) == 1, file.exists(dat_path))
    ext <- tolower(tools::file_ext(dat_path))
    
    guess_sep <- function(path, enc = encoding) {
      con <- if (is.null(enc)) file(path, "r") else file(path, "r", encoding = enc)
      on.exit(close(con), add = TRUE)
      ln <- readLines(con, n = 1)
      if (length(ln) == 0) return("")
      if (grepl(",", ln)) return(",")
      if (grepl("\t", ln)) return("\t")
      " "  
    }
    
    # Auto-detect if file has header
    detect_header <- function(path, sep_val, enc = encoding) {
      con <- if (is.null(enc)) file(path, "r") else file(path, "r", encoding = enc)
      on.exit(close(con), add = TRUE)
      
      # Read first two lines
      lines <- readLines(con, n = 2)
      if (length(lines) < 2) return(TRUE)  
      
      # Split by separator
      line1 <- strsplit(lines[1], sep_val)[[1]]
      line2 <- strsplit(lines[2], sep_val)[[1]]
      
      # Remove empty strings and whitespace
      line1 <- trimws(line1[nzchar(trimws(line1))])
      line2 <- trimws(line2[nzchar(trimws(line2))])
      
      if (length(line1) == 0 || length(line2) == 0) return(TRUE)
      
      # Count numeric values in each line
      numeric_count_1 <- sum(suppressWarnings(!is.na(as.numeric(line1))))
      numeric_count_2 <- sum(suppressWarnings(!is.na(as.numeric(line2))))
      
      # If first line has significantly fewer numbers than second line, likely header
      # If both lines are mostly numeric, likely no header
      ratio_1 <- numeric_count_1 / length(line1)
      ratio_2 <- numeric_count_2 / length(line2)
      
      # Logic:
      # - If line1 < 50% numeric and line2 > 80% numeric → has header
      # - If line1 > 80% numeric and line2 > 80% numeric → no header
      if (ratio_1 < 0.5 && ratio_2 > 0.8) {
        return(TRUE)   # Has header
      } else if (ratio_1 > 0.8 && ratio_2 > 0.8) {
        return(FALSE)  # No header
      } else {
        return(TRUE)   # Default: assume header
      }
    }
    
    # 1) read with auto-detected header
    has_header <- TRUE  # Default
    sep_val <- " "      # Default
    
    if (ext %in% c("csv")) {
      sep_val <- ","
      has_header <- detect_header(dat_path, sep_val, encoding)
      rawdat <- if (is.null(encoding)) {
        utils::read.csv(dat_path, header = has_header, stringsAsFactors = FALSE)
      } else {
        utils::read.csv(dat_path, header = has_header, stringsAsFactors = FALSE, fileEncoding = encoding)
      }
    } else if (ext %in% c("tsv", "tab")) {
      sep_val <- "\t"
      has_header <- detect_header(dat_path, sep_val, encoding)
      rawdat <- if (is.null(encoding)) {
        utils::read.delim(dat_path, header = has_header, stringsAsFactors = FALSE)
      } else {
        utils::read.delim(dat_path, header = has_header, stringsAsFactors = FALSE, fileEncoding = encoding)
      }
    } else if (ext %in% c("txt", "dat", "")) {
      sep_val <- guess_sep(dat_path, encoding)
      if (sep_val == "") sep_val <- " "
      has_header <- detect_header(dat_path, sep_val, encoding)
      rawdat <- if (is.null(encoding)) {
        utils::read.table(dat_path, header = has_header, sep = sep_val, stringsAsFactors = FALSE)
      } else {
        utils::read.table(dat_path, header = has_header, sep = sep_val, stringsAsFactors = FALSE, fileEncoding = encoding)
      }
    } else if (ext %in% c("rds")) {
      rawdat <- as.data.frame(readRDS(dat_path), stringsAsFactors = FALSE)
      has_header <- TRUE
    } else if (ext %in% c("rdata", "rda")) {
      # Load RData file
      cat("Loading RData file...\n")
      env <- new.env()
      loaded_objects <- load(dat_path, envir = env)
      
      cat("  Objects in RData:", paste(loaded_objects, collapse = ", "), "\n")
      
      # Check what was loaded
      if (length(loaded_objects) == 0) {
        stop("RData file is empty")
      }
      
      # If multiple objects, try common names first
      if (length(loaded_objects) > 1) {
        # Look for common data frame names
        common_names <- c("data", "rawdat", "df", "dataset", "dat")
        found <- intersect(common_names, loaded_objects)
        if (length(found) > 0) {
          obj_name <- found[1]
          rawdat <- as.data.frame(get(obj_name, envir = env), stringsAsFactors = FALSE)
          cat("  Using object:", obj_name, "\n")
        } else {
          # Use the first data frame found
          rawdat_found <- FALSE
          for (obj_name in loaded_objects) {
            obj <- get(obj_name, envir = env)
            if (is.data.frame(obj) || is.matrix(obj)) {
              rawdat <- as.data.frame(obj, stringsAsFactors = FALSE)
              cat("  Using object:", obj_name, "\n")
              rawdat_found <- TRUE
              break
            }
          }
          if (!rawdat_found) {
            stop("No data frame found in RData file. Available objects: ", paste(loaded_objects, collapse = ", "))
          }
        }
      } else {
        # Only one object
        obj_name <- loaded_objects[1]
        obj <- get(obj_name, envir = env)
        if (is.data.frame(obj) || is.matrix(obj)) {
          rawdat <- as.data.frame(obj, stringsAsFactors = FALSE)
          cat("  Using object:", obj_name, "\n")
        } else {
          stop("Object '", obj_name, "' is not a data frame or matrix. Type: ", class(obj)[1])
        }
      }
      has_header <- TRUE
    } else if (ext %in% c("xlsx", "xls")) {
      if (!requireNamespace("readxl", quietly = TRUE)) stop("readxl is required for Excel files.")
      rawdat <- as.data.frame(readxl::read_excel(dat_path), stringsAsFactors = FALSE)
      has_header <- TRUE
    } else if (ext %in% c("sav", "sas7bdat", "dta")) {
      if (!requireNamespace("haven", quietly = TRUE)) stop("haven is required for SPSS/SAS/Stata files.")
      if (ext == "sav")      rawdat <- as.data.frame(haven::read_sav(dat_path), stringsAsFactors = FALSE)
      if (ext == "sas7bdat") rawdat <- as.data.frame(haven::read_sas(dat_path), stringsAsFactors = FALSE)
      if (ext == "dta")      rawdat <- as.data.frame(haven::read_dta(dat_path), stringsAsFactors = FALSE)
      has_header <- TRUE
    } else {
      sep_val <- guess_sep(dat_path, encoding)
      if (sep_val == "") sep_val <- " "
      has_header <- detect_header(dat_path, sep_val, encoding)
      rawdat <- if (is.null(encoding)) {
        utils::read.table(dat_path, header = has_header, sep = sep_val, stringsAsFactors = FALSE)
      } else {
        utils::read.table(dat_path, header = has_header, sep = sep_val, stringsAsFactors = FALSE, fileEncoding = encoding)
      }
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
