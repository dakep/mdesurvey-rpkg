#' Extract Information from Fitted Minimum Divergence Estimators
#'
#' Methods to extract coefficient estimates and variance/covariance components from
#' minimum divergence estimators fitted to complex survey data.
#'
#' @details
#' `coef()` extracts the parameter estimates from the fitted MDE and the coefficient estimates
#' from the regression MDE.
#'
#' @param object an MDE object fitted to survey data
#' @param ... currently unused.
#' @importFrom rlang check_dots_empty
#' @rdname generics
#' @export
coef.survey_mde <- function (object, ...) {
  check_dots_empty()
  object$estimates
}

#' @importFrom rlang check_dots_empty
#' @rdname generics
#' @export
coef.survey_reg_mde <- function (object, ...) {
  check_dots_empty()
  object$coefficients
}

#' @details
#' `sigma()` extracts the standard deviation from the fitted regression model.
#'
#' @importFrom stats sigma
#' @rdname generics
#' @export
sigma.survey_reg_mde <- function (object, ...) {
  check_dots_empty()
  object$nuisance
}

#' @details
#' `vcov()` calculates the variance-covariance matrix for fitted MDEs.
#'
#' @param type which type of covariance matrix to compute.
#' @param n either how to estimate the effective sample size, or (if numeric) the
#'   sample size.
#' @param ... currently unused.
#' @importFrom rlang check_dots_empty
#' @name vcov
#' @rdname generics
#' @export
vcov.survey_mde <- function (object, type = c("sandwich", "model"),
                             n = c("score", "kish"), ...) {
  check_dots_empty()
  type_missing <- missing(type)
  if (!is.numeric(n)) {
    n <- switch(match.arg(n),
                score = object$neff_score,
                kish  = object$neff_kish,
                object$nobs)
  }
  .vcov(object, type = match.arg(type), n = n, enforce_type = !type_missing)
}

#'
#' @param which which components to include in the covariance matrix: only the regression
#'   coefficients (`"regression"`), only the nuisance parameters (`"nuisance"`), or
#'   all (`"all"`).
#' @rdname generics
#' @export
vcov.survey_reg_mde <- function (object, type = c("sandwich", "model"),
                                 which = c('regression', 'nuisance', 'all'),
                                 n = c("score", "kish"), ...) {
  check_dots_empty()
  covmat_ind <- switch(match.arg(which),
                       regression = seq_along(object$coefficients),
                       nuisance   = length(object$coefficients) + seq_along(object$nuisance),
                       seq_len(length(object$coefficients) + length(object$nuisance)))

  type_missing <- missing(type)
  if (!is.numeric(n)) {
    n <- switch(match.arg(n),
                score = object$neff_score,
                kish  = object$neff_kish,
                object$nobs)
  }

  if (length(n) == length(covmat_ind)) {
    n <- c(n, rep.int(1, nrow(object$finf_inv) - length(n)))
  }

  # finf_inv <- object$finf_inv[covmat_ind, covmat_ind, drop = FALSE]
  covmat <- .vcov(object, type = match.arg(type), n = n, enforce_type = !type_missing,
                  finf_inv = object$finf_inv)
  covmat[covmat_ind, covmat_ind, drop = FALSE]
}
