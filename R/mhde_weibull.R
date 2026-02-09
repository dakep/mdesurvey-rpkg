#' Minimum Hellinger Distance Estimator for the Weibull Model
#'
#' Compute the MHDE for the shape \eqn{k} and scale \eqn{\lambda} of the Weibull model.
#'
#' @inheritParams survey_mhde
#' @param initial an initial estimate for the shape and scale parameters. If missing,
#'   use the method of moments.
#' @param sandwich_cov if `TRUE`, compute the sandwich estimator for the covariance matrix.
#'   If `FALSE`, uses the inverse of the Fisher information under independence and the
#'   effective number of observations estimated from the scores.
#'   The only reason for setting this to `FALSE` is to increase computation speed when the
#'   covariance matrix is not needed.
#'   Otherwise, use [vcov()] to get the desired estimator of the covariance matrix.
#' @importFrom stats dweibull
#' @importFrom rlang warn enquo
#' @family Minimum Hellinger Distance Estimator
#' @export
mhd_weibull <- function (x, design, initial, na.rm = FALSE,
                         sandwich_cov = TRUE,
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

  euler_gamma <- 0.57721566490153286060651209008240
  pisq_6 <- 1.6449340668482264364724151666460 # pi^2/6


  if (missing(initial)) {
    initial <- mom_weibull(svy$x, wgts)$estimates
  }

  if (!is.null(names(initial))) {
    initial <- initial[c("shape", "scale")]
  } else {
    names(initial) <- c("shape", "scale")
  }

  raw_scores <- function (x, est, z) {
    if (missing(z)) {
      z <- x / est[['scale']]
    }
    cbind(1 / est[['shape']] + log(z) * (1 - z^est[['shape']]),
          est[['shape']] * (z^est[['shape']] - 1) / est[['scale']])
  }

  fisher_inf <- function (est, ...) {
    matrix(c((1 - 2 * euler_gamma + euler_gamma^2 + pisq_6) / est[['shape']],
             (euler_gamma - 1) / est[['scale']],
             (euler_gamma - 1) / est[['scale']],
             est[['shape']]^2 / est[['scale']]^2),
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
        z <- xint / estimates[['scale']]
        score <- raw_scores(z = z, est = estimates)
        hess <- cbind(0.5 * score[, 1]^2 -
                        1 / estimates[['shape']]^2 - log(z)^2 * z^estimates[['shape']],
                      0.5 * score[, 1] * score[, 2] +
                        (z^estimates[['shape']] - 1 +
                           estimates[['shape']] * z^estimates[['shape']] * log(z)) / estimates[['scale']],
                      0.5 * score[, 2]^2 +
                        (estimates[['shape']] / estimates[['scale']]^2) *
                        (1 - (1 + estimates[['shape']]) * z^estimates[['shape']]))

        f_theta <- dweibull(xint, shape = estimates[['shape']], scale = estimates[['scale']])
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

  obj <- survey_mhde(x, design,
                     initial,
                     model_dfun = \(x, parameters, log) {
                       dweibull(x, shape = parameters[['shape']], scale = parameters[['scale']],
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
  class(obj) <- c("survey_mde_weibull", class(obj))

  obj
}

#' @rdname vcov
#' @export
vcov.survey_mde_weibull <- function (object, type = c("sandwich", "model"),
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


#' Method of Moments Estimator for the Weibull Model
#'
#' @param x univariate observations from the finite population
#' @param wgts sampling weights
#' @return a list with `estimates`, the covariance `cov` and the density function.
#' @importFrom stats dweibull weighted.mean uniroot
#' @keywords internal
mom_weibull <- function (x, wgts, shape_range = c(0.02, 100), cov = FALSE) {
  mu <- weighted.mean(x, w = wgts)
  sig2 <- weighted.mean((x - mu)^2, w = wgts)
  cv2_hat <- sig2 / mu^2

  shape_opt <- tryCatch({
    uniroot(f = \(k) {
      g1 <- gamma(1 + 1/k)
      g2 <- gamma(1 + 2/k)
      theo_cv2 <- (g2 / g1^2) - 1

      theo_cv2 - cv2_hat
    }, interval = c(0.02, 100))
  }, error = \(cnd) list(root = NA_real_))

  k_est <- shape_opt$root

  if (!is.na(k_est)) {
    lambda_est <- mu / gamma(1 + 1/k_est)
  } else {
    lambda_est <- NA_real_
  }

  estimates <- c(shape = k_est, scale = lambda_est)

  covest <- if (isTRUE(cov)) {
    euler_gamma <- 0.57721566490153286060651209008240
    pisq_6 <- 1.6449340668482264364724151666460 # pi^2/6
    covest <- matrix(c((1 - 2 * euler_gamma + euler_gamma^2 + pisq_6) / estimates[['shape']],
                       (euler_gamma - 1) / estimates[['scale']],
                       (euler_gamma - 1) / estimates[['scale']],
                       estimates[['shape']]^2 / estimates[['scale']]^2),
                     ncol = 2) |>
      solve()
  } else {
    NULL
  }

  list(estimates = estimates,
       cov       = covest)
}
