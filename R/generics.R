#' Extract Model Coefficients
#'
#' Extract coefficients from a fitted MDE.
#'
#' @export
coef.survey_mde <- function (object, ...) {
  object$estimates
}
