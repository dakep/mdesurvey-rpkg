#' @importFrom stats dgamma pgamma
#' @importFrom survey svymean
#' @include model_family.R
Gamma <- ModelFamily$new(
  name            = "Gamma",
  parameter_names = c("shape", "scale"),
  trans           = log,
  inv_trans       = exp,
  range = matrix(c(0, Inf), ncol = 2),
  dfun  = \(x, params, log) {
    dgamma(x, shape = params[['shape']], scale = params[['scale']], log = log)
  },
  pfun = \(x, params) {
    pgamma(x, shape = params[['shape']], scale = params[['scale']])
  },
  raw_scores = \(x, params) {
    cbind(log(x) - log(params[['scale']]) - digamma(params[['shape']]),
          (x / params[['scale']] - params[['shape']]) / params[['scale']])
  },
  hessian = \(x, params) {
    cbind(
      -trigamma(params[['shape']]), # shape^2
      params[['shape']] / params[['scale']]^2 - 2 * x / params[['scale']]^3, # scale^2
      -1 / params[['scale']] # shape,scale
    )
  },
  fisher_inf = \(params) {
    matrix(c(trigamma(params[['shape']]),
             1 / params[['scale']],
             1 / params[['scale']],
             params[['shape']] / params[['scale']]^2),
           ncol = 2)
  },
  initial = \(x, design) {
    wmx <- svymean(x, design = design)
    loglik_root <- log(wmx) - svymean(log(x), design = design)
    shape <- uniroot(\(a) log(a) - digamma(a) - loglik_root,
                     interval = c(1e-6, max(x)))

    c(shape = shape$root, scale = wmx / shape$root)
  },
  # For the Gamma model, the nuisance parameter is the scale (standard deviation)
  nuisance_names = 'sd',
  default_link   = link('log', sd = 'log'),
  parameters_from_mean_par = \(mean, nuisance) {
    if (!isTRUE(mean > .Machine$double.eps) || !isTRUE(nuisance[[1]] > .Machine$double.eps)) {
      c(shape = NA_real_, scale = NA_real_)
    } else {
      c(shape = mean^2 / nuisance[[1]]^2, scale = nuisance[[1]]^2 / mean)
    }
  },
  mean_par = \ (params) {
    c(mean = prod(params),
      sd   = sqrt(params[['shape']]) * params[['scale']])
  },
  jacobian_mean_par_mapping = \ (mean, nuisance) {
    matrix(c(2 * mean / nuisance^2, -(nuisance/mean)^2,
             -2 * mean^2 / nuisance^3, 2 * nuisance / mean), ncol = 2)
  }
)

.model_family_register$gamma <- Gamma
