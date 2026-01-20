#' Minimum Hellinger Distance Estimator for the Lognormal Model
#'
#' Compute the MHDE for the parameters \eqn{\mu,\sigma} of the Lognormal model.
#'
#' @param x univariate observations from the finite population
#' @param wgts sampling weights
#' @param initial an initial estimate for the shape and scale parameters. If missing,
#'   use the MLE.
#' @param bw bandwidth for the HT-adjusted KDE. By default,
#'   use Scott's rule (see [stats::bw.nrd()]).
#' @param log_transform whether to log-transform the data before estimation or not.
#' @param integration_subdivisions number of partitions to divide the domain of \eqn{\hat f()}
#'   for Gauss-Kronrod quadrature.
#' @param optim_method,optim_control method and control options passed on to [stats::optim()].
#' @importFrom stats dnorm dlnorm weighted.mean uniroot optim
#' @importFrom rlang warn
#' @family Minimum Hellinger Distance Estimator
#' @seealso [stats::dlnorm()] for the parametrization by `meanlog` \eqn{\mu} and
#'   `sdlog` \eqn{\sigma}.
#' @export
mhd_lognorm <- function (x, wgts = NULL, initial, log_transform = FALSE,
                         bw,
                         integration_subdivisions = 256,
                         optim_method = "Nelder-Mead", optim_control = list()) {
  if (isTRUE(log_transform) && any(x < .Machine$double.eps)) {
    warn("Cannot fit the log-normal distribution on the log scale if 0's are present")
    log_transform <- FALSE
  }
  logx <- log(x[x > 0])
  if (is.null(wgts)) {
    wgts <- rep.int(1 / length(x), length(x))
  }
  if (missing(initial)) {
    initial <- c(weighted.mean(logx, w = wgts), sd(logx))
  }

  if (isTRUE(log_transform)) {
    x <- log(x)
    rg <- c(-Inf, Inf)
    dfun_factory <- \(params) {
      \(x) dnorm(x, mean = params[[1]], sd = params[[2]], log = TRUE)
    }
  } else {
    rg <- c(0, Inf)
    dfun_factory <- \(params) {
      \(x) dlnorm(x, meanlog = params[[1]], sdlog = params[[2]], log = TRUE)
    }
  }

  environment(dfun_factory) <- environment(mhd_lognorm)

  if (missing(bw)) {
    bw <- 1.06 * mad(x) * length(x)^(-1/5)
  }

  mhde_integral <- hd_gauss_quadrature(x, wgts,
                                       bandwidth = bw,
                                       range = rg,
                                       n_subdivisions = integration_subdivisions)

  initial[[2]] <- log(initial[[2]])
  mhd_est <- optim(initial, \(params) {
    params[[2]] <- exp(params[[2]])
    mhde_integral(dfun_factory(params))
  },
  method = optim_method,
  control = optim_control)

  if (mhd_est$convergence != 0) {
    warn(sprintf("Optimizer did not converge (code %d): %s",
                        mhd_est$convergence, paste0('', mhd_est$message)))
  }

  estimates <- c("mean" = mhd_est$par[[1]],
                 "sd"   = exp(mhd_est$par[[2]]))
  covest <- matrix(c(estimates[['sd']]^2, 0, 0, 0.5 * estimates[['sd']]^2), ncol = 2)

  structure(
    list(estimates      = estimates,
         bias           = numeric(2),
         cov            = covest,
         mhd            = mhd_est$value,
         initial        = c(initial[[1]], exp(initial[[2]])),
         dfun           = dlnorm,
         model          = "gamma",
         optimizer_code = mhd_est$convergence,
         optimizer_msg  = mhd_est$message),
    class = 'survey_mde')
}
