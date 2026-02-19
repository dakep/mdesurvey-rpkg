#' @importFrom stats dgamma pgamma
#' @importFrom survey svymean svyvar
#' @include model_family.R
Normal <- ModelFamily$new(
  .register       = TRUE,
  name            = "Normal",
  parameter_names = c("mean", "sd"),
  trans           = \(params) { params[[2]] <- log(params[[2]]); params; },
  inv_trans       = \(u_params) { u_params[[2]] <- exp(u_params[[2]]); u_params; },
  range = matrix(c(-Inf, Inf), ncol = 2),
  dfun  = \(x, params, log) {
    dnorm(x, mean = params[['mean']], sd = params[['sd']], log = log)
  },
  pfun = \(x, params) {
    pnorm(x, mean = params[['mean']], sd = params[['sd']])
  },
  raw_scores = \(x, params) {
    x[] <- (x - params[['mean']]) / params[['sd']]
    cbind(x / params[['sd']],
          ((x^2 - 1) / params[['sd']]))
  },
  hessian = \(x, params) {
    x[] <- (x - params[['mean']]) / params[['sd']]
    cbind(-1 / params[['sd']]^2,            # mean x mean
          (1 - 3 * x^2) / params[['sd']]^2, # sd x sd,
          -2 * x / params[['sd']]^2)        # mean x sd
  },
  fisher_inf = \(params) {
    matrix(c(1 / params[['sd']]^2, 0, 0, 2 / params[['sd']]^2),
           ncol = 2)
  },
  initial = \(x, design) {
    c(svymean(x, design = design),
      sqrt(svyvar(x, design = design)))
  }
)

