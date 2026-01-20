#' Minimum Hellinger Distance Estimator for the Weibull Model
#'
#' Compute the MHDE for the shape \eqn{k} and scale \eqn{\lambda} of the Weibull model.
#'
#' @param x univariate observations from the finite population
#' @param wgts sampling weights
#' @param initial an initial estimate for the shape and scale parameters. If missing,
#'   use the method of moments.
#' @param integration_subdivisions number of partitions to divide the domain of \eqn{\hat f()}
#'   for Gauss-Kronrod quadrature.
#' @param bw bandwidth for the HT-adjusted KDE. By default,
#'   use Scott's rule (see [stats::bw.nrd()]).
#' @param optim_method,optim_control method and control options passed on to [stats::optim()].
#' @importFrom stats dweibull uniroot optim
#' @importFrom rlang warn
#' @family Minimum Hellinger Distance Estimator
#' @export
mhd_weibull <- function (x, wgts = NULL, initial, integration_subdivisions = 256,
                         bw, optim_method = "Nelder-Mead", optim_control = list()) {
  if (is.null(wgts)) {
    wgts <- rep.int(1 / length(x), length(x))
  }
  if (missing(initial)) {
    initial <- mom_weibull(x, wgts)
  }

  if (missing(bw)) {
    bw <- 1.06 * sd(x) * length(x)^(-1/5)
  }
  mhde_integral <- hd_gauss_quadrature(x, wgts,
                                       bandwidth = bw,
                                       range = c(0, Inf),
                                       n_subdivisions = integration_subdivisions)

  mhd_est <- optim(log(initial), \(params) {
    params <- exp(params)
    mhde_integral(\(x) {
      dweibull(x, shape = params[[1]], scale = params[[2]], log = TRUE)
    })
  },
  method = optim_method,
  control = optim_control)

  names(mhd_est$par) <- c("shape", "scale")

  if (mhd_est$convergence != 0) {
    warn(sprintf("Optimizer did not converge (code %d): %s",
                 mhd_est$convergence, paste0('', mhd_est$message)))
  }

  euler_gamma <- 0.57721566490153286060651209008240
  pisq_6 <- 1.6449340668482264364724151666460 # pi^2/6
  covest <- matrix(c((1 - 2 * euler_gamma + euler_gamma^2 + pisq_6) / mhd_est$par[['shape']],
                     (euler_gamma - 1) / mhd_est$par[['scale']],
                     (euler_gamma - 1) / mhd_est$par[['scale']],
                     mhd_est$par[['shape']]^2 / mhd_est$par[['scale']]^2),
                   ncol = 2) |>
    solve()

  list(estimates      = exp(mhd_est$par),
       mhd            = mhd_est$value,
       cov            = covest,
       dfun           = dweibull,
       initial        = initial,
       optimizer_code = mhd_est$convergence,
       optimizer_msg  = mhd_est$message)
}


#' Method of Moments Estimator for the Weibull Model
#'
#' @param x univariate observations from the finite population
#' @param wgts sampling weights
#' @return a list with `estimates`, the covariance `cov` and the density function.
#' @importFrom stats dweibull weighted.mean uniroot
#' @keywords internal
mom_weibull <- function (x, wgts, shape_range = c(0.02, 100)) {
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

  c(shape = k_est, scale = lambda_est)
}
