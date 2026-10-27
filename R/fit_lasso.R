#' fit_lasso
#' @description Elastic-net logistic fit at a single lambda, with a convergence check.
#'   glmnet can fail to converge at a single lambda (fit$jerr != 0); it then returns an
#'   empty model (intercept 0, all coefficients 0) with only a warning. In that case refit
#'   along a warm-start lambda path ending at \code{lambda}; if that also fails, return
#'   \code{fallback}.
#' @param x design matrix (n x p)
#' @param y binary response
#' @param w observation weights
#' @param alpha elastic-net mixing parameter
#' @param lambda penalty
#' @param family glmnet family
#' @param fallback coefficients to return if glmnet cannot converge
#' @return numeric vector c(intercept, coefficients)
#' @keywords internal
fit_lasso <- function(x, y, w, alpha, lambda, family = "binomial", fallback = NULL) {
  fit <- glmnet::glmnet(x, y, family = family, weights = w, alpha = alpha, lambda = lambda)
  if (fit$jerr == 0) return(c(as.numeric(fit$a0), as.numeric(fit$beta)))

  # single-lambda fit failed: approach lambda along a warm-start path
  path <- exp(seq(log(100 * lambda), log(lambda), length.out = 50))
  fit <- glmnet::glmnet(x, y, family = family, weights = w, alpha = alpha, lambda = path)
  k <- length(fit$lambda)
  if (fit$jerr == 0 && k == length(path)) {
    return(c(as.numeric(fit$a0[k]), as.numeric(fit$beta[, k])))
  }

  if (!is.null(fallback)) return(fallback)
  stop("fit_lasso: glmnet did not converge")
}
