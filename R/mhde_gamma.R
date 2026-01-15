#' Minimum Hellinger Distance Estimator for the Gamma Model
#'
#' Compute the MHDE for the shape \eqn{\alpha} and scale \eqn{\theta} of the Gamma model.
#'
#' @param x univariate observations from the finite population
#' @param wgts sampling weights
#' @param initial an initial estimate for the shape and scale parameters. If missing,
#'   use the MLE.
#' @param cov_type whether to use the sandwich estimator for the covariance matrix
#'   the inverse Fisher information under the model, or not estimate the covariance matrix.
#' @param integration_subdivisions number of partitions to divide the domain of \eqn{\hat f()}
#'   for Gauss-Kronrod quadrature.
#' @param optim_method,optim_control method and control options passed on to [stats::optim()].
#' @importFrom stats dgamma weighted.mean uniroot optim
#' @importFrom rlang warn
#' @family Minimum Hellinger Distance Estimator
#' @seealso [stats::dgamma()] for the parametrization by shape and scale.
#' @export
mhd_gamma <- function (x, wgts = NULL, initial, cov_type = c("sandwich", "model", "none"),
                       integration_subdivisions = 256,
                       optim_method = "Nelder-Mead", optim_control = list()) {
  cov_type <- match.arg(cov_type)
  if (is.null(wgts)) {
    wgts <- rep.int(1 / length(x), length(x))
  } else if (length(wgts) == 1L) {
    wgts <- rep.int(wgts, length(x))
  }
  if (missing(initial)) {
    initial <- gamma_mle(x, wgts)$estimates
  }

  bandwidth <- 1.06 * sd(x) * length(x)^(-1/5)

  mhde_integral <- hd_gauss_quadrature(x, wgts, bandwidth, range = c(0, Inf),
                                       n_subdivisions = integration_subdivisions)

  mhd_est <- optim(log(initial), \(params) {
    params <- exp(params)
    mhde_integral(\(xint) {
      dgamma(xint, shape = params[[1]], scale = params[[2]], log = TRUE)
    })
  },
  method = optim_method,
  control = optim_control)

  if (mhd_est$convergence != 0) {
    warn(sprintf("Optimizer did not converge (code %d): %s",
                 mhd_est$convergence, paste0('', mhd_est$message)))
  }

  estimates <- exp(mhd_est$par)
  names(estimates) <- c("shape", "scale")

  if (identical(cov_type, "sandwich")) {
    # Sandwich estimator of the covariance matrix
    A_est_els <- mhde_integral(log = FALSE, \(xint) {
      score <- cbind(log(xint) - log(estimates[['scale']]) - digamma(estimates[['shape']]),
                     (xint/estimates[['scale']] - estimates[['shape']]) / estimates[['scale']])

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

    # Combine the two steps
    # (1) A_est_els = -0.5 * A_est_els + 1  # because mhde_integral returns 2 - 2 * A_est_els
    # (2) A_est = -0.5 * A_est_els
    # into one
    A_est <- matrix(A_est_els[c(1, 2, 2, 3)], ncol = 2) / 4 - 0.5

    fhat <- vapply(x, FUN.VALUE = numeric(1), FUN = \(gp) {
      u <- (x - gp) / bandwidth
      sum(wgts * epanechnikov_kernel(u))
    }) / (bandwidth * sum(wgts))

    score <- cbind(log(x) - log(estimates[['scale']]) - digamma(estimates[['shape']]),
                   (x / estimates[['scale']] - estimates[['shape']]) / estimates[['scale']])

    sigma_hat <- cov(0.25 * score *
                       sqrt(dgamma(x, shape = estimates[['shape']], scale = estimates[['scale']])) /
                       sqrt(fhat))

    covest <- tryCatch({
      solve(A_est, sigma_hat) %*% solve(A_est)
    }, error = \(cnd) {
      warning("Sandwich estimator is singular. Returning model-based covariance estimate.")
      matrix(c(trigamma(estimates[[1]]),
               1 / estimates[[2]],
               1 / estimates[[2]],
               estimates[[1]] / estimates[[2]]^2),
             ncol = 2) |>
        solve()
    })

    estimated_bias <- c(
      3 / (length(x) * 2 * estimates[['shape']]^2 * psigamma(estimates[['shape']], deriv = 2)),
      -estimates[['scale']] * (
        3 * estimates[['shape']] * psigamma(estimates[['shape']], deriv = 2) +
          psigamma(estimates[['shape']], deriv = 1) +
          estimates[['shape']] * psigamma(estimates[['shape']], deriv = 1)^2) /
        (length(x) * estimates[['shape']] *
           (estimates[['shape']] * psigamma(estimates[['shape']], deriv = 2) - 1)))

    ## Bias correction ##
    # bias_vec <- mhde_integral(log = FALSE, \(xint) {
    #   score <- cbind(log(xint) - log(estimates[['scale']]) - digamma(estimates[['shape']]),
    #                  (xint/estimates[['scale']] - estimates[['shape']]) / estimates[['scale']])
    #
    #   # hessian elements [1, 1], [1,2], and [2, 2]
    #   hess <- cbind(0.5 * score[, 1]^2 -
    #                   trigamma(estimates[['shape']]),
    #                 0.5 * score[, 1] * score[, 2] -
    #                   1 / estimates[['scale']],
    #                 0.5 * score[, 2]^2 +
    #                   estimates[['shape']] / estimates[['scale']]^2 -
    #                   2 * xint / estimates[['scale']]^3)
    #
    #   hess_col <- function (j, k) {
    #     col <- c(1, 2, 2, 3)[[(j - 1) * 2 + k]]
    #     hess[, col]
    #   }
    #
    #   third_score <- array(dim = c(2, 2, 2, length(xint)), data = 0)
    #   third_score[1, 1, 1, ] <- -psigamma(estimates[['shape']], deriv = 2)
    #   third_score[2, 2, 2, ] <- 6 * xint / estimates[['scale']]^4 - 2 * estimates[['shape']] / estimates[['scale']]^3
    #   third_score[1, 2, 2, ] <- 1 / estimates[['scale']]^2
    #
    #   third_deriv <- matrix(0, nrow = length(xint), ncol = 2^3)
    #
    #   for (j in 1:2) {
    #     for (k in 1:2) {
    #       for (l in 1:2) {
    #         col <- (j - 1) * (2^2) + (k - 1) * 2 + l
    #         third_deriv[, col] <- third_score[j, k, l, ] +
    #           0.5 * (score[, j] * hess_col(k, l) + score[, k] * hess_col(j, l) + score[, l] * hess_col(j, k)) +
    #           0.25 * score[, j] * score[, k] * score[, l]
    #       }
    #     }
    #   }
    #
    #   f_theta <- dgamma(xint, shape = estimates[['shape']], scale = estimates[['scale']])
    #   0.5 * sqrt(f_theta) * third_deriv
    # })
    #
    # bias_mat <- array(0.5 * bias_vec - 1, dim = rep.int(2, 3))
    # ds <- apply(bias_mat, 3, \(b) sum(colSums(covest * b)))
    #
    # estimated_bias <- -solve(A_est, ds) / length(x)
  } else if (identical(cov_type, "model")) {
    covest <- matrix(c(trigamma(estimates[[1]]),
                       1 / estimates[[2]],
                       1 / estimates[[2]],
                       estimates[[1]] / estimates[[2]]^2),
                     ncol = 2) |>
      solve()

    estimated_bias <- c(
      3 / (length(x) * 2 * estimates[['shape']]^2 * psigamma(estimates[['shape']], deriv = 2)),
      -estimates[['scale']] * (
        3 * estimates[['shape']] * psigamma(estimates[['shape']], deriv = 2) +
          psigamma(estimates[['shape']], deriv = 1) +
          estimates[['shape']] * psigamma(estimates[['shape']], deriv = 1)^2) /
        (length(x) * estimates[['shape']] *
           (estimates[['shape']] * psigamma(estimates[['shape']], deriv = 2) - 1)))
  } else {
    covest <- NULL
    estimated_bias <- numeric(2)
  }

  structure(
    list(estimates      = estimates,
         bias           = estimated_bias,
         cov            = covest,
         mhd            = mhd_est$value,
         initial        = initial,
         dfun           = dgamma,
         model          = "gamma",
         optimizer_code = mhd_est$convergence,
         optimizer_msg  = mhd_est$message),
    class = 'survey_mde')
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
