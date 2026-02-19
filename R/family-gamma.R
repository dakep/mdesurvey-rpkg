#' @importFrom stats dgamma pgamma
#' @importFrom survey svymean
#' @include model_family.R
Gamma <- ModelFamily$new(
  .register       = TRUE,
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

    c(shape$root, wmx / shape$root)
  }
)
