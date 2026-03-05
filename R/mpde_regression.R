#' Minimum Phi Distance Estimator for the Linear Regression Models with Complex
#' Survey Designs
#'
#' Compute the survey MPDE for the mean regression parameters under a user-specified continuous
#' model for the conditional distribution \eqn{Y \mid X = x} for **discrete** \eqn{X}.
#'
#' @param x a formula or matrix specifying the linear model.
#'   The objects are first looked up in the provided `design`, then in the calling environment.
#' @param design the survey design created by [survey::svydesign()] and friends.
#' @param initial a list with two components:
#'   `coef` with the initial estimates for the regression coefficients and
#'   `nuisance` with the initial values for the nuisance components.
#'   The family's `nuisance_name` has information about the nuisance parameters for the
#'   respective model family.
#' @param family either the name of a known model family, or a
#'   model family as returned by [model_family()].
#' @param divergence either the name of a known phi divergence, or a
#'   phi divergence as returned by [phi_divergence()].
#' @param na.rm a logical evaluating to `TRUE` or `FALSE` indicating whether `NA` values
#'   should be omitted.
#' @param integration_subdivisions number of partitions to divide the domain of \eqn{\hat f()}
#'   for Gaussian quadrature.
#' @param bw bandwidth for the HT-adjusted KDE. Either a scalar bandwidth (not advised)
#'   used for all groups, or (the name of) a function which will be used to cmpute the
#'   bandwidth separately in each group. Uses [stats::bw.nrd()] as default.
#' @param kernel Name of the kernel used in the HT-adjusted KDE.
#'
#' @param optim_method,optim_control method and control options passed on to [stats::optim()].
#' @importFrom stats uniroot optim bw.nrd model.frame model.matrix model.response lm.fit
#' @importFrom stats na.omit na.fail
#' @importFrom rlang warn abort enquo
#' @importFrom survey svydesign
#' @importFrom methods is
#' @family Minimum Phi-Divergence Estimator
#' @export
survey_regression_mpde <- function (x,
                                    design,
                                    initial,
                                    initial_var,
                                    family,
                                    divergence = 'HD',
                                    na.rm = FALSE,
                                    integration_subdivisions = 256,
                                    bw = bw.nrd,
                                    kernel = c("epanechnikov", "triangular", "rectangular", "biweight"),
                                    optim_method = "Nelder-Mead", optim_control = list()) {
  if (is.character(divergence)) {
    divergence <- phi_divergence(divergence)
  }

  if (!is.R6(divergence) || !is(divergence, "PhiDivergence")) {
    abort("`divergence` must be a valid phi divergence from phi_divergence()")
  }

  kernel <- match.arg(kernel)

  if (is.character(family)) {
    family <- model_family(family)
  }

  if (!is.R6(family) || !is(family, "ModelFamily")) {
    abort("`family` must be a valid model family from model_family()")
  }

  # Define groups
  mf <- model.frame(x, design$variables,
                    na.action = if (isTRUE(na.rm)) na.omit else na.fail)
  tms <- terms(mf)
  y <- model.response(mf) |>
    unname()
  x_terms <- attr(tms, "variables")[-c(1, (1 + attr(tms, "response")))] |>
    as.character()

  xmat <- model.matrix(mf, design$variables)

  if (!all(match(x_terms, names(attr(xmat, "contrasts")), nomatch = 0L) > 0L)) {
    abort("All predictors must be of type `factor` or `character`.")
  }

  if (any(colSums(xmat) < 1)) {
    abort("Design matrix is not full rank")
  }
  grp <- apply(xmat, 1, paste, collapse = '|') |>
    unname()

  tot_wgts <- sum(weights(design))
  grouped <- split(seq_along(y), grp) |>
    lapply(\(i) {
      list(y    = y[i],
           design = subset_svydesign(design, i),
           prop = sum(weights(design)[i]) / tot_wgts,
           x    = xmat[i[[1]], ])
    }) |>
    unname()


  bwf <- if (missing(bw) || is.null(bw)) {
    bw.nrd
  } else if (!is.numeric(bw)) {
    match.fun(bw)
  } else {
    \(...) {
      as.numeric(bw[[1]])
    }
  }

  grouped <- lapply(grouped, \(gr) {
    gr$bw <- bwf(gr$y)
    c(gr,
      phidiv_gauss_quadrature(gr$y,
                              wgts = weights(gr$design),
                              divergence = divergence,
                              family = family,
                              bandwidth = gr$bw,
                              kernel = kernel,
                              n_subdivisions = integration_subdivisions))
  })

  initial_params <- if (missing(initial) || is.null(initial)) {
    grp_init <- vapply(grouped, FUN.VALUE = numeric(2), FUN = \(gr) {
      subset_init <- family$initial(gr$y, gr$design)
      family$mean_par(subset_init)
    })

    # Use the geometric mean of the group variances
    c(rowMeans(log(grp_init[-1, , drop = FALSE])),
      lm.fit(x = do.call(rbind, lapply(grouped, `[[`, 'x')),
             y = grp_init['mean', ])$coefficients)
  } else if (length(initial) != 2) {
    if (is.null(initial$coef)) {
      abort("`initial` does not contain information about the coefficients")
    }
    if (length(family$nuisance_names) > 0L && is.null(initial$nuisance)) {
      abort(paste("`initial` does not contain information about the nuisance parameters",
                  paste(family$nuisance_names, collapse = ", ")))
    }
  } else if (!is.null(names(initial$coefs))) {
    c(as.numeric(initial$nuisance), as.numeric(initial$coefs[colnames(xmat)]))
  } else {
    c(as.numeric(initial$nuisance), as.numeric(initial$coefs))
  }
  names(initial_params) <- NULL
  nuisance_components <- seq_along(family$nuisance_names)
  mean_components <- length(family$nuisance_names) + seq_len(ncol(xmat))
  mpd_est <- optim(initial_params,
                   fn = \(regpars) {
                     vapply(grouped, FUN.VALUE = numeric(1), FUN = \(gr) {
                       mean <- drop(gr$x %*% regpars[mean_components])
                       # We need to transform here because the `divergence_int()` back-transforms!
                       gr_params <- family$trans(
                         family$parameters_from_mean_par(mean, exp(regpars[nuisance_components])))
                       if (all(is.finite(gr_params))) {
                         gr$prop * gr$divergence_int(gr_params)
                       } else {
                         Inf
                       }
                     }) |>
                       sum()
                   },
                   method = optim_method,
                   control = optim_control)

  if (mpd_est$convergence != 0) {
    warn(sprintf("Optimizer did not converge (code %d): %s",
                 mpd_est$convergence, paste0('', mpd_est$message)))
  }

  initial_params[nuisance_components] <- exp(initial_params[nuisance_components])
  names(initial_params) <- c(family$nuisance_names, colnames(xmat))
  reg_coefs <- mpd_est$par[mean_components]
  names(reg_coefs) <- colnames(xmat)
  nuisance_est <- exp(mpd_est$par[nuisance_components])
  names(nuisance_est) <- family$nuisance_names

  # Covariance matrix estimation
  vapply(grouped, FUN.VALUE = numeric(1), FUN = \(gr) {
    mean <- drop(gr$x %*% reg_coefs)
    estimates <- family$parameters_from_mean_par(mean, nuisance_est) |>
      setNames(family$parameter_names)
    A_est <- gr$sandwich_cov_A(estimates)
    score <- family$raw_scores(gr$y, estimates)
    jac <- family$jacobian_mean_par_mapping(mean, nuisance_est)
    score_jac <- score %*% jac
    colnames(score_jac) <- c(".mean", family$nuisance_names)

    score_reg <- outer(score_jac[, 1L], gr$x) |>
      cbind(score_jac[, -1, drop = FALSE])
  })

  structure(
    list(coefficients   = reg_coefs,
         nuisance       = nuisance_est,
         family         = family,
         divergence     = divergence,
         mpd            = mpd_est$value,
         initial        = initial,
         optimizer_code = mpd_est$convergence,
         optimizer_msg  = mpd_est$message),
    class = 'survey_reg_mde')
}

subset_svydesign <- function (x, i) {
  d <- x[i, ]
  d$fpc$sampsize <- d$fpc$sampsize[i, ]
  d$dcheck <- lapply(d$dcheck, \(dc) {
    dc$id <- dc$id[i]
    dc$dcheck <- dc$dcheck[i, i]
    dc
  })
  d
}

#' @importFrom survey svymean
#' @importFrom stats vcov as.formula
#' @importFrom rlang warn
reg_mpde_covest <- function (x, estimates, integrator, family, design, divergence) {
  A_est <- integrator$sandwich_cov_A(estimates)
  score <- family$raw_scores(x, estimates)
  design$variables <- data.frame(score = score)
  fmla <- as.formula(paste0("~`", paste0(colnames(design$variables), collapse = "`+`"), "`"))
  design_vcov <- vcov(svymean(fmla, design = design))
  finf <- family$fisher_inf(estimates)
  neff <- diag(finf) / diag(design_vcov)

  covest <- tryCatch({
    (solve(A_est, design_vcov) %*% solve(A_est) * divergence$phi_2nd_deriv_at_1^2) |>
      structure(type = "sandwich")
  }, error = \(cnd) {
    # Compute the model-based estimator instead.
    warn("Sandwich covariance estimate is singular. Computing model-based estimate instead.",
         parent = cnd)

    finf_inv <- tryCatch(solve(finf),
                         error = \(cnd) {
                           # Try a Moore-Penrose pseudo-inverse
                           tryCatch(pseudo_inverse(finf),
                                    error = \(cnd) {
                                      warn("Fisher Information matrix is not positive semi-definite.",
                                           parent = cnd)
                                      NULL
                                    })
                         })

    (finf_inv / outer(sqrt(neff), sqrt(neff))) |>
      structure(type = "model")
  })

  dimnames(covest) <- list(family$parameter_names, family$parameter_names)
  names(neff) <- family$parameter_names
  list(cov = covest,
       neff = neff)
}
