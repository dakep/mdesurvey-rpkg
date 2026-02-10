#' Minimum Hellinger Distance Estimator for the General Continuous Models with Complex
#' Survey Designs
#'
#' Compute the survey MHDE for the parameters of a user-specified continuous model.
#'
#' # Integration function
#'
#' The second argument to the `cov_fun` parameter is a function with signature
#' `function(g, log = FALSE, non_negative_integrand = TRUE)` that computes
#' \deqn{\int \sqrt{ \hat f(x) g(x) }\;\mathrm{d}x,}
#' or (if `log=TRUE`)
#' \deqn{\int \exp\{\frac{1}{2} (\log(\hat f(x)) + \log(g(x)) )\}\;\mathrm{d}x,}
#' for arbitrary functions \eqn{g()}.
#' The `integration` function is set up such that \eqn{\hat f} is evaluated a fixed number
#' of times in each subdivision, and does not need to be re-evaluated for different \eqn{g()}.
#'
#' @param x a formula, symbol, expression, vector, or matrix specifying the observations.
#'   The objects are first looked up in the provided `design`, then in the calling environment.
#' @param design the survey design created by [survey::svydesign()] and friends.
#' @param initial an initial estimate for the parameters.
#'   If the vector has names, those will be used for the returned estimates.
#' @param model_dfun the model density function with signature `function(x, parameters, log) {...}`,
#'   where `x` are the evaluation points, `parameters` are the parameter values and `log` is
#'   a logical value if the density should be evaluated on the log scale or not.
#' @param model_domain a matrix with 2 columns giving the lower and upper endpoints of the domain
#'   for each dimension. Can be infinite.
#' @param parameter_transform,parameter_transform_inv transformation function and its inverse
#'   to transform the parameters to the entire real space for easier optimization.
#'   These functions must take one argument (the parameters/transformed parameters) and return
#'   the parameters on the transformed/the original scale.
#' @param cov_fun function to compute the covariance matrix for the estimator.
#'   The signature of the function must be `function(estimates, integration)`, where
#'   `estimates` are the estimated parameters and `integration` is a function that integrates
#'   any given function times the KDE over the domain (see Details below).
#' @param na.rm a logical evaluating to `TRUE` or `FALSE` indicating whether `NA` values
#'   should be omitted.
#' @param integration_subdivisions number of partitions to divide the domain of \eqn{\hat f()}
#'   for Gauss-Kronrod quadrature.
#' @param bw bandwidth for the HT-adjusted KDE. Uses [stats::bw.nrd()] as default.
#' @param kernel Name of the kernel used in the HT-adjusted KDE.
#'
#' @param optim_method,optim_control method and control options passed on to [stats::optim()].
#' @importFrom stats uniroot optim bw.nrd
#' @importFrom rlang warn abort enquo
#' @family Minimum Hellinger Distance Estimator
#' @export
survey_mhde <- function (x, design,
                         initial,
                         na.rm = FALSE,
                         model_dfun,
                         model_domain,
                         parameter_transform,
                         parameter_transform_inv,
                         cov_fun,
                         integration_subdivisions = 256,
                         bw,
                         kernel = c("epanechnikov", "triangular", "rectangular", "biweight"),
                         optim_method = "Nelder-Mead", optim_control = list()) {
  x <- enquo(x)
  svy <- .extract_survey_values(!!x, design)
  wgts <- svy$wgts

  if (dim(svy$x)[[2]] > 1L) {
    abort("Only univariate values are supported.")
  }
  svy$x <- svy$x[,1]
  nobs <- length(svy$x)
  nobs_eff <- sum(wgts)^2 / sum(wgts^2)

  kernel <- match.arg(kernel)

  param_names <- names(initial)

  if (missing(bw) || is.null(bw)) {
    bw <- bw.nrd(svy$x)
  }

  mhde_integral <- hd_gauss_quadrature(svy$x, wgts,
                                       bandwidth = bw,
                                       range = model_domain,
                                       kernel = kernel,
                                       n_subdivisions = integration_subdivisions)

  mhd_est <- optim(parameter_transform(initial), \(params) {
    params <- parameter_transform_inv(params)
    hellinger_affinity <- mhde_integral(\(xint) {
      model_dfun(xint, parameters = params, log = TRUE)
    })
    2 - 2 * hellinger_affinity
  },
  method = optim_method,
  control = optim_control)

  if (mhd_est$convergence != 0) {
    warn(sprintf("Optimizer did not converge (code %d): %s",
                 mhd_est$convergence, paste0('', mhd_est$message)))
  }

  estimates <- parameter_transform_inv(mhd_est$par)
  names(estimates) <- param_names
  covest <- cov_fun(estimates, mhde_integral)
  dimnames(covest) <- list(param_names, param_names)

  structure(
    list(estimates      = estimates,
         cov            = covest,
         mhd            = mhd_est$value,
         initial        = initial,
         nobs           = nobs,
         neff_kish      = nobs_eff,
         optimizer_code = mhd_est$convergence,
         optimizer_msg  = mhd_est$message),
    class = 'survey_mde')
}
