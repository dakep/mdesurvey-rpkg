#' Extract Model Coefficients
#'
#' Extract coefficients from a fitted MDE.
#'
#' @export
coef.survey_mde <- function (object, ...) {
  object$estimates
}

#' Calculate Variance-Covariance Matrix for Minimum Divergence Estimators
#'
#' @param object a MDE fitted to survey data.
#' @param type which type of covariance matrix to compute.
#' @param n either how to estimate the effective sample size, or (if numeric) the
#'   sample size.
#' @param ... currently unused.
#' @name vcov
#' @export
vcov.survey_mde <- function (object, type = c("sandwich", "model"),
                              n = c("score", "kish"), ...) {
  type_missing <- missing(type)
  if (!is.numeric(n)) {
    n <- switch(match.arg(n),
                score = object$neff_score,
                kish  = object$neff_kish,
                object$nobs)
  }
  .vcov(object, type = match.arg(type), n = n, enforce_type = !type_missing)
}
