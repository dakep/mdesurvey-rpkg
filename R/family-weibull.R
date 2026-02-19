.k_euler_gamma <- 0.57721566490153286060651209008240
.k_pisq_6 <- 1.6449340668482264364724151666460 # pi^2/6

#' @importFrom stats dweibull pweibull
#' @importFrom survey svymean
#' @include model_family.R
Weibull <- ModelFamily$new(
  .register       = TRUE,
  name            = "Weibull",
  parameter_names = c("shape", "scale"),
  trans           = log,
  inv_trans       = exp,
  range = matrix(c(0, Inf), ncol = 2),
  dfun = \(x, params, log) {
    dweibull(x, shape = params[['shape']], scale = params[['scale']], log = log)
  },
  pfun = \(x, params) {
    pweibull(x, shape = params[['shape']], scale = params[['scale']])
  },
  raw_scores = \(x, params) {
    x[] <- x / params[['scale']]
    cbind(1 / params[['shape']] + log(x) * (1 - x^params[['shape']]),
          params[['shape']] * (x^params[['shape']] - 1) / params[['scale']])
  },
  hessian = \(x, params) {
    x[] <- x / params[['scale']]

    cbind(-1 / params[['shape']]^2 - log(x)^2 * x^params[['shape']], # shape^2
          (params[['shape']] / params[['scale']]^2) *
            (1 - (1 + params[['shape']]) * x^params[['shape']]),     # scale^2
          (x^params[['shape']] - 1 + params[['shape']] * x^params[['shape']] * log(x)) /
            params[['scale']] # shape,scale
    )
  },
  fisher_inf = \(params) {
    matrix(c((1 - 2 * .k_euler_gamma + .k_euler_gamma^2 + .k_pisq_6) / params[['shape']],
             (.k_euler_gamma - 1) / params[['scale']],
             (.k_euler_gamma - 1) / params[['scale']],
             params[['shape']]^2 / params[['scale']]^2),
           ncol = 2)
  },
  initial = \(x, design) {
    # A method-of-moments initial estimate
    mu <- svymean(x, design = design)
    sig2 <- svymean((x - mu)^2, design = design)
    cv2_hat <- sig2 / mu^2

    shape_opt <- tryCatch({
      uniroot(f = \(k) {
        g1 <- gamma(1 + 1/k)
        g2 <- gamma(1 + 2/k)
        theo_cv2 <- (g2 / g1^2) - 1

        theo_cv2 - cv2_hat
      }, interval = c(1e-8, max(x)))
    }, error = \(cnd) {
      warn("Cannot compute initial estimate for the Weibull family", parent = cnd)
      list(root = NA_real_)
    })

    k_est <- shape_opt$root

    if (!is.na(k_est)) {
      lambda_est <- mu / gamma(1 + 1/k_est)
    } else {
      lambda_est <- NA_real_
    }

    c(k_est, lambda_est)
  }
)
