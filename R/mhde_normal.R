#' Minimum Hellinger Distance Estimator for the Normal Model
#'
#' Compute the MHDE for the parameters \eqn{\mu,\sigma} of the Normal model.
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
#' @importFrom stats dnorm
#' @importFrom rlang warn enquo
#' @family Minimum Hellinger Distance Estimator
#' @export
mhd_norm <- function (x, design, initial,
                      sandwich_cov = TRUE,
                      na.rm = FALSE,
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
    initial <- c(mean = svymean(svy$x, design = design),
                 sd   = sqrt(svyvar(svy$x, design = design)))
  }

  if (!is.null(names(initial))) {
    initial <- initial[c("mean", "sd")]
  } else {
    names(initial) <- c("mean", "sd")
  }

  raw_scores <- function (x, est, z) {
    if (missing(z)) {
      z <- (x - est[['mean']]) / est[['sd']]
    }
    cbind(z / est[['sd']],
          ((z^2 - 1) / est[['sd']]))
  }

  fisher_inf <- function (estimates) {
    matrix(c(1 / estimates[['sd']]^2,
             0,
             0,
             2 / estimates[['sd']]^2),
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
        z <- (xint - estimates[['mean']]) / estimates[['sd']]
        score <- raw_scores(z = z, est = estimates)
        hess <- cbind(0.5 * score[, 1]^2 -
                        1 / estimates[['sd']]^2,
                      0.5 * score[, 1] * score[, 2] -
                        2 * z / estimates[['sd']]^2,
                      0.5 * score[, 2]^2 +
                        (1 - 3 * z^2) / estimates[['sd']]^2)

        f_theta <- dnorm(xint, mean = estimates[[1]], sd = estimates[[2]], log = FALSE)
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

  obj <- survey_mhde(svy$x, design,
                     initial,
                     model_dfun = \ (x, parameters, log) {
                       dnorm(x, mean = parameters[[1]], sd = parameters[[2]], log = log)
                     },
                     cov_fun                  = covfun,
                     model_domain             = c(-Inf, Inf),
                     parameter_transform      = \(x) { c(x[[1]], log(x[[2]])) },
                     parameter_transform_inv  = \(x) { c(x[[1]], exp(x[[2]])) },
                     integration_subdivisions = integration_subdivisions,
                     bw                       = bw,
                     optim_method             = optim_method,
                     optim_control            = optim_control)

  obj$fisher_inf <- fisher_inf(obj$estimates)
  obj$neff_scores <- attr(obj$cov, "neff")
  names(obj$neff_scores) <- names(obj$estimates)
  attr(obj$cov, "neff") <- NULL
  class(obj) <- c("survey_mde_norm", class(obj))

  obj
}

#' @rdname vcov
#' @export
vcov.survey_mde_norm <- function (object, type = c("sandwich", "model"),
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
