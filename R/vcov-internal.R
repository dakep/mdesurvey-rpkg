#' Calculate Variance-Covariance Matrix for Minimum Divergence Estimators
#'
#' @param object a MDE fitted to survey data
#' @param type which type of covariance matrix to compute.
#' @param n either how to estimate the effective sample size, or (if numeric) the
#'   sample size.
#' @param enforce_type should the `type` be enforced, or use the model-based estimate
#'   as fallback with a warning.
#' @importFrom stats vcov coef
#' @importFrom rlang abort
#' @keywords internal
.vcov <- function (object, type, n, enforce_type = FALSE, finf_inv) {
  if (!missing(type) && !is.null(type)) {
    if (identical(attr(object$cov, "type"), type)) {
      if (is.character(object$cov)) {
        abort(paste0("Covariance estimate of type \"", type, "\" is ", object$cov, "."))
      }
      return(object$cov)
    } else if (enforce_type && !identical(type, "model")) {
      warn(paste0("Covariance estimate of type \"", type, "\" is unavalable. ",
                  "Returning model-based estimate instead."))
    }
  }

  if (missing(n) || is.null(n) || !is.numeric(n)) {
    n <- object$neff_kish
  }

  if (missing(finf_inv) || is.null(finf_inv)) {
    finf_inv <- solve(object$family$fisher_inf(coef(object)))
  }

  if (length(n) == 1L) {
    finf_inv / n
  } else if (length(n) == nrow(finf_inv)) {
    finf_inv / outer(sqrt(n), sqrt(n))
  } else {
    abort(sprintf("`n` must be either of length 1 or %d (= the number of parameters)",
                  nrow(finf_inv)))
  }
}
