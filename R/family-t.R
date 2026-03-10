#' @importFrom stats dt pt
#' @importFrom survey svymean svyvar
#' @include model_family.R link-functions.R
StudentT <- ModelFamily$new(
  name            = "Student-t",
  parameter_names = c("location", "scale", "df"),
  trans           = \(params) { params[2:3] <- log(params[2:3]); params },
  inv_trans       = \(u_params) { u_params[2:3] <- exp(u_params[2:3]); u_params },
  range = matrix(c(-Inf, Inf), ncol = 2),
  dfun  = \(x, params, log) {
    dens <- dt((x - params[['location']]) / params[['scale']], df = params[['df']], log = log)
    if (isTRUE(log)) {
      dens - log(params[['scale']])
    } else {
      dens / params[['scale']]
    }
  },
  pfun = \(x, params) {
    pt((x - params[['location']]) / params[['scale']], df = params[['df']])
  },

  raw_scores = \(x, params) {
    x[] <- (x - params[['location']]) / params[['scale']]
    df_p_x2 <- params[['df']] + x^2
    matrix(c(
      (params[['df']] + 1) * x / (df_p_x2 * params[['scale']]),
      params[['df']] * (x^2 - 1) / (df_p_x2 * params[['scale']]),
      0.5 * (
        digamma(0.5 * (params[['df']] + 1))
        - digamma(0.5 * params[['df']])
        - 1 / params[['df']]
        - log(1 + x^2 / params[['df']])
        + x^2 * (params[['df']] + 1) / (params[['df']] * df_p_x2))
    ),
    ncol = 3L)
  },

  hessian = \(x, params) {
    x[] <- (x - params[['location']]) / params[['scale']]
    x2 <- x^2
    denom <- params[['df']] + x2

    H <- matrix(0, nrow = length(x), ncol = 6L)

    # location x location
    H[, 1L] <- ((params[['df']] + 1) / params[['scale']]^2) * ((x2 - params[['df']]) / denom^2)

    # scale x scale
    H[, 2L] <- (params[['df']] * (1 - x2)) / (params[['scale']]^2 * denom) -
      (2 * params[['df']] * (params[['df']] + 1) * x2) / (params[['scale']]^2 * denom^2)

    # df x df
    H[, 3L] <- 0.25 * (trigamma(0.5 * (params[['df']]+1)) - trigamma(0.5 * params[['df']])) +
      1/(2 * params[['df']]^2) +
      x2 / (2 * params[['df']] * denom) -
      (x2 * (params[['df']]^2 + 2 * params[['df']] + x2)) / (2 * params[['df']]^2 * denom^2)

    # location x scale
    H[, 4L] <- - (2 * params[['df']] * (params[['df']] + 1) * x) / (params[['scale']]^2 * denom^2)

    # location x df
    H[, 5L] <- (x * (x2 - 1)) / (params[['scale']] * denom^2)

    # scale x df
    H[, 6L] <- (x2 * (x2 - 1)) / (params[['scale']] * denom^2)

    H
  },
  fisher_inf = \(params) {
    finf <- matrix(0, ncol = 3, nrow = 3)

    # location
    finf[[1, 1]] <- ((params[['df']] + 1) / (params[['df']] + 3)) / params[['scale']]^2
    # scale
    finf[[2, 2]] <- (2 * params[['df']] / (params[['df']] + 3)) / params[['scale']]^2
    # df
    finf[[3, 3]] <- 0.25 * trigamma(params[['df']]/2) - 0.25 * trigamma((params[['df']]+1)/2) -
      (params[['df']] + 5) / (2 * params[['df']] * (params[['df']] + 1) * (params[['df']] + 3))

    # scale x df
    finf[[2, 3]] <- finf[[3, 2]] <-
      -2 / (params[['scale']] * (params[['df']] + 1) * (params[['df']] + 3))

    dimnames(finf) <- list(c("location", "scale", "df"),
                           c("location", "scale", "df"))
    finf
  },
  #' @importFrom stats median mad
  initial = \(x, design) {
    med <- median(x)
    c(location = med,
      scale    = mad(x, center = med, constant = 1),
      df       = 4)
  },
  # For the Normal model, the nuisance parameter is the scale (standard deviation)
  nuisance_names = c('scale', 'df'),
  default_link   = link('identity', scale = 'log', df = 'log'),
  parameters_from_mean_par = \(mean, nuisance) {
    nuisance[!is.finite(nuisance)] <- NA_real_
    c(location = mean, scale = nuisance[[1]], df = nuisance[[2]])
  },
  mean_par = \(params) {
    names(params)[[1]] <- 'mean'
    params
  },
  jacobian_mean_par_mapping = \(mean, nuisance) {
    diag(nrow = 3, ncol = 3)
  }
)

.model_family_register$student_t <-
  .model_family_register$studentt <-
  .model_family_register$t <-
  .model_family_register$tdist <- StudentT
