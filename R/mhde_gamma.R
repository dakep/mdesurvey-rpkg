#' Minimum Hellinger Distance Estimator for the Gamma Model
#'
#' Compute the MHDE for the shape \eqn{\alpha} and scale \eqn{\theta} of the Gamma model.
#'
#' @inheritParams survey_mhde
#' @param initial an initial estimate for the shape and scale parameters. If missing,
#'   use the MLE.
#' @param sandwich_cov if `TRUE`, compute the sandwich estimator for the covariance matrix.
#'   If `FALSE`, uses the inverse of the Fisher information under independence and the
#'   effective number of observations estimated from the scores.
#'   The only reason for setting this to `FALSE` is to increase computation speed when the
#'   covariance matrix is not needed.
#'   Otherwise, use [vcov()] to get the desired estimator of the covariance matrix.
#' @importFrom stats dgamma weighted.mean vcov
#' @importFrom rlang warn abort enquo
#' @importFrom survey svydesign svymean
#' @family Minimum Hellinger Distance Estimator
#' @seealso [stats::dgamma()] for the parametrization by shape and scale.
#' @export
mhd_gamma <- function (x, design, initial,
                       na.rm = FALSE, sandwich_cov = TRUE,
                       integration_subdivisions = 256, bw,
                       optim_method = "Nelder-Mead", optim_control = list()) {
  x <- enquo(x)
  svy <- .extract_survey_values(!!x, design, na.rm = na.rm)
  wgts <- svy$wgts

  if (dim(svy$x)[[2]] > 1L) {
    abort("Only univariate values are supported")
  }
  svy$x <- svy$x[,1]
  nobs <- length(svy$x)

  if (missing(initial)) {
    initial <- mle_gamma(svy$x, wgts)$estimates
  }
  if (!is.null(names(initial))) {
    initial <- initial[c("shape", "scale")]
  } else {
    names(initial) <- c("shape", "scale")
  }

  raw_scores <- function (x, est) {
    cbind(log(x) - log(est[['scale']]) - digamma(est[['shape']]),
          (x / est[['scale']] - est[['shape']]) / est[['scale']])
  }

  fisher_inf <- function (estimates) {
    matrix(c(trigamma(estimates[['shape']]),
             1 / estimates[['scale']],
             1 / estimates[['scale']],
             estimates[['shape']] / estimates[['scale']]^2),
           ncol = 2)
  }

  covfun <- if (isFALSE(sandwich_cov)) {
    function (estimates, mhde_integral) {
      score <- raw_scores(svy$x, estimates)
      design$variables <- data.frame(score = score)
      design_vcov <- vcov(svymean(~ score.1 + score.2, design = design))
      finf <- fisher_inf(estimates)
      neff <- diag(finf) / diag(design_vcov)

      covest <- solve(finf) / outer(sqrt(neff), sqrt(neff))

      dimnames(covest) <- list(names(estimates), names(estimates))
      attr(covest, "neff") <- neff
      attr(covest, "type") <- 'model'
      covest
    }
  } else {
    function (estimates, mhde_integral) {
      # Sandwich estimator of the covariance matrix
      A_est_els <- mhde_integral(log = FALSE, non_negative_integrand = FALSE, \(xint) {
        score <- raw_scores(xint, estimates)
        hess <- cbind(0.5 * score[, 1]^2 -
                        trigamma(estimates[['shape']]),
                      0.5 * score[, 1] * score[, 2] -
                        1 / estimates[['scale']],
                      0.5 * score[, 2]^2 +
                        estimates[['shape']] / estimates[['scale']]^2 -
                        2 * xint / estimates[['scale']]^3)

        f_theta <- dgamma(xint, shape = estimates[['shape']], scale = estimates[['scale']])
        sqrt(f_theta) * hess
      })

      A_est <- matrix(0.5 * A_est_els[c(1, 2, 2, 3)], ncol = 2)
      score <- raw_scores(svy$x, estimates)
      design$variables <- data.frame(score = score)
      design_vcov <- vcov(svymean(~ score.1 + score.2, design = design))

      covest <- tryCatch({
        solve(A_est, design_vcov) %*% solve(A_est) / 16
      }, error = \(cnd) {
        warning("Sandwich estimator is singular.")
        'singular'
      })

      dimnames(covest) <- list(names(estimates), names(estimates))
      attr(covest, "neff") <- diag(fisher_inf(estimates)) / diag(design_vcov)
      attr(covest, "type") <- 'sandwich'
      covest
    }
  }

  obj <- survey_mhde(!!x,
                     design,
                     initial,
                     model_dfun = \(x, parameters, log) {
                       dgamma(x, shape = parameters[['shape']], scale = parameters[['scale']],
                              log = log)
                     },
                     cov_fun                  = covfun,
                     model_domain             = c(0, Inf),
                     parameter_transform      = log,
                     parameter_transform_inv  = exp,
                     integration_subdivisions = integration_subdivisions,
                     bw                       = bw,
                     optim_method             = optim_method,
                     optim_control            = optim_control)

  obj$fisher_inf <- fisher_inf(obj$estimates)
  obj$neff_scores <- attr(obj$cov, "neff")
  names(obj$neff_scores) <- names(obj$estimates)
  attr(obj$cov, "neff") <- NULL
  class(obj) <- c("survey_mde_gamma", class(obj))

  obj
}

#' @rdname vcov
#' @export
vcov.survey_mde_gamma <- function (object, type = c("sandwich", "model"),
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

#' Maximum Likelihood Estimator for the Gamma Model
#'
#' @param x univariate observations from the finite population
#' @param wgts sampling weights
#' @return a list with `estimates`, the covariance `cov` and the density function.
#' @importFrom stats dgamma weighted.mean uniroot
#' @keywords internal
mle_gamma <- function (x, wgts = NULL) {
  avg <- if (is.null(wgts)) {
    mean
  } else {
    \(x) weighted.mean(x, w = wgts)
  }
  wmx <- avg(x)
  loglik_root <- log(wmx) - avg(log(x))
  shape <- uniroot(\(a) log(a) - digamma(a) - loglik_root,
                   interval = c(1e-3, 1e3))

  estimates <- c(shape = shape$root, scale = wmx / shape$root)

  covest <- matrix(c(trigamma(estimates[[1]]),
                     1 / estimates[[2]],
                     1 / estimates[[2]],
                     estimates[[1]] / estimates[[2]]^2),
                   ncol = 2) |>
    solve()

  list(estimates = estimates,
       cov       = covest,
       dfun      = dgamma)
}
