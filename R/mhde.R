#' Minimum Hellinger Distance Estimator for the General Continuous Models with Complex
#' Survey Designs
#'
#' Compute the survey MHDE for the parameters in a given model family.
#'
#'
#' @param x a formula, symbol, expression, vector, or matrix specifying the observations.
#'   The objects are first looked up in the provided `design`, then in the calling environment.
#' @param design the survey design created by [survey::svydesign()] and friends.
#' @param initial an initial estimate for the parameters.
#'   If the vector has names, those will be used for the returned estimates.
#' @param family either the name of a known model family, or a
#'   model family as returned by [model_family()].
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
#' @importFrom R6 is.R6
#' @family Minimum Hellinger Distance Estimator
#' @export
survey_mhde <- function (x, design,
                         family,
                         initial,
                         na.rm = FALSE,
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

  if (is.character(family)) {
    family <- model_family(family)
  }

  if (!is.R6(family) || !is(family, "ModelFamily")) {
    abort("`family` must be a valid model family from model_family()")
  }

  if (missing(initial)) {
    initial <- family$initial(svy$x, design = design)
  }
  if (is.null(names(initial))) {
    names(initial) <- family$parameter_names
  }
  initial <- initial[family$parameter_names]

  if (missing(bw) || is.null(bw)) {
    bw <- bw.nrd(svy$x)
  }

  mhde_integral <- hd_gauss_quadrature(svy$x,
                                       wgts = wgts,
                                       bandwidth = bw,
                                       range = family$range[1, , drop = TRUE],
                                       kernel = kernel,
                                       n_subdivisions = integration_subdivisions)

  mhd_est <- optim(family$trans(initial), \(params) {
    params <- family$inv_trans(params)
    hellinger_affinity <- mhde_integral$hellinger_affinity(\(xint) {
      family$dfun(xint, params = params, log = TRUE)
    })
    2 - 2 * hellinger_affinity
  },
  method = optim_method,
  control = optim_control)

  if (mhd_est$convergence != 0) {
    warn(sprintf("Optimizer did not converge (code %d): %s",
                 mhd_est$convergence, paste0('', mhd_est$message)))
  }

  estimates <- family$inv_trans(mhd_est$par)
  names(estimates) <- family$parameter_names
  covest <- mhde_covest(svy$x, estimates, mhde_integral$int_sqrt_f, family, design)

  structure(
    list(estimates      = estimates,
         cov            = covest$cov,
         mhd            = mhd_est$value,
         family         = family,
         initial        = initial,
         nobs           = nobs,
         neff_scores    = covest$neff,
         neff_kish      = nobs_eff,
         optimizer_code = mhd_est$convergence,
         optimizer_msg  = mhd_est$message),
    class = c('survey_mhde', 'survey_mde'))
}

#' @importFrom survey svymean
#' @importFrom stats vcov
#' @importFrom rlang warn
mhde_covest <- function (x, estimates, int_fun, family, design) {
  # Sandwich estimator of the covariance matrix
  cov_entries <- matrix(rep(seq_along(estimates), each = 2L), nrow = 2L) |>
    cbind(combn(length(estimates), 2L))

  A_est_els <- int_fun(\(xint) {
    score <- family$raw_scores(xint, estimates)
    hess <- family$hessian(xint, estimates)

    A_int <- seq_len(ncol(cov_entries)) |>
      vapply(FUN.VALUE = score[, 1], FUN = \(ceind) {
        ij <- cov_entries[, ceind]
        0.5 * score[, ij[[1]]] * score[, ij[[2]]] + hess[, ceind]
      })

    f_theta <- family$dfun(xint, estimates, log = FALSE)
    0.5 * sqrt(f_theta) * A_int
  })

  A_est <- matrix(data = NA_real_, nrow = length(estimates), ncol = length(estimates))
  for (ceind in seq_len(ncol(cov_entries))) {
    A_est[[cov_entries[1, ceind], cov_entries[2, ceind]]] <-
      A_est[[cov_entries[2, ceind], cov_entries[1, ceind]]] <-
      A_est_els[[ceind]]
  }
  score <- family$raw_scores(x, estimates)
  design$variables <- data.frame(score = score)
  fmla <- as.formula(paste0("~`", paste0(colnames(design$variables), collapse = "`+`"), "`"))
  design_vcov <- vcov(svymean(fmla, design = design))
  finf <- family$fisher_inf(estimates)
  neff <- diag(finf) / diag(design_vcov)

  covest <- tryCatch({
    (solve(A_est, design_vcov) %*% solve(A_est) / 16) |>
      structure(type = "sandwich")
  }, error = \(cnd) {
    # Compute the model-based estimator instead.
    warn("Sandwich covariance estimate is singular. Computing model-based estimate instead.",
         parent = cnd)

    (solve(finf) / outer(sqrt(neff), sqrt(neff))) |>
      structure(type = "model")
  })

  dimnames(covest) <- list(family$parameter_names, family$parameter_names)
  names(neff) <- family$parameter_names
  list(cov = covest,
       neff = neff)
}

#' Functor to Create Fixed-grid Gaussian Quadrature For Hellinger Affinity
#'
#' This creates a function to compute
#' \deqn{\int \sqrt{ \hat f(x) g(x) }\;\mathrm{d}x,}
#' or (if `log=TRUE`)
#' \deqn{\int \exp\{\frac{1}{2} (\log(\hat f(x)) + \log(g(x)) )\}\;\mathrm{d}x,}
#' for arbitrary non-negative functions \eqn{g()}.
#' The function is set up such that \eqn{\hat f} is evaluated exactly 21 times in each
#' subdivision, and does not need to be re-evaluated for different \eqn{g()}.
#'
#'
#' @return a function to evaluate the integral
#'   \deqn{\int \sqrt{ \hat f(x) g(x) }\;\mathrm{d}x}
#'   \deqn{\int \sqrt{ \hat f(x) g(x) }\;\mathrm{d}x}
#'
#' @keywords internal
#' @importFrom SparseGrid createSparseGrid
#' @rdname hd_gauss_quadrature
hd_gauss_quadrature <- function (x, wgts, bandwidth, n_subdivisions = 256, poly_order = 20,
                                 range, kernel) {
  kernel <- substr(kernel[[1]], 1, 1)

  # Evaluate the KDE at the Gaussian evaluation points
  # We don't need to go beyond the range ±bandwidth because the KDE
  # will be 0 outside.
  range_x <- range(x) + c(-bandwidth, bandwidth)
  if (!missing(range)) {
    range_x[[1]] <- max(range_x[[1]], range[[1]])
    range_x[[2]] <- min(range_x[[2]], range[[2]])
  }

  subintervals <- partition_domain(x, subdivisions = n_subdivisions, bw = bandwidth,
                                   range = range_x)

  int_pts <- createSparseGrid('GQU', 1, poly_order)

  gkpts <- lapply(subintervals, \(subint) {
    lwr <- subint[[1]]
    upr <- subint[[2]]
    matrix(c((upr - lwr) * int_pts$nodes + lwr,
             (upr - lwr) * int_pts$weights),
           ncol = 2)
  }) |>
    do.call(what = rbind)

  colnames(gkpts) <- c("int_x", "int_wgt")
  gkpts <- cbind(gkpts,
                 log_f_hat = .kde(x       = x,
                                  wgts    = wgts,
                                  evalpts = gkpts[, "int_x"],
                                  bw      = bandwidth,
                                  kernel  = kernel) |>
                   log())

  list(
    # Compute int exp{0.5 [log(f(x)) + h(x)]} dx
    # for non-negative h(x)
    hellinger_affinity = \(h) {
      h <- h(gkpts[, "int_x"])
      if (is.null(dim(h))) {
        h <- matrix(h, ncol = 1)
      }

      vapply(seq_len(ncol(h)), FUN.VALUE = numeric(1), FUN = \(i) {
        crossprod(gkpts[, "int_wgt"], exp(0.5 * (gkpts[, "log_f_hat"] + h[, i])))[[1, 1]]
      }) |>
        pmax(0)
    },

    # Compute int sqrt(f(x)) h(x) dx
    # for arbitrary h(x)
    int_sqrt_f = \(h) {
      h <- h(gkpts[, "int_x"])
      if (is.null(dim(h))) {
        h <- matrix(h, ncol = 1)
      }

      vapply(seq_len(ncol(h)), FUN.VALUE = numeric(1), FUN = \(i) {
        crossprod(gkpts[, "int_wgt"], exp(0.5 * gkpts[, "log_f_hat"]) * h[, i])[[1, 1]]
      })
    }
  )
}
