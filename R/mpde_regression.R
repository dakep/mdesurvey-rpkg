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
#'   The initial values must be on the scale of the link function!
#' @param family either the name of a known model family, or a
#'   model family as returned by [model_family()].
#' @param link the link functions obtained with [link()] for the mean and nuisance components
#'   of the linear regression model.
#'   If missing, uses the default link for the family `family$default_link`.
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
                                    family,
                                    link,
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

  link <- if (missing(link)) {
    family$default_link
  } else {
    .check_link(link, family$nuisance_names)
  }

  if (identical(design$variance, "YG")) {
    warn(paste("The Yates-Grundy (YG) variance estimator is not suppored for MPD regression estimators.",
               "Using Horvitz-Thompson (HT) instead."))
    design$variance <- "HT"
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

  p <- ncol(xmat)
  d <- length(family$nuisance_names)

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
      sub_design <- subset_svydesign(design, i)
      list(y    = y[i],
           design = sub_design,
           prop = sum(weights(sub_design)) / tot_wgts,
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
    grp_init <- vapply(grouped, FUN.VALUE = numeric(1L + d), FUN = \(gr) {
      subset_init <- family$initial(gr$y, gr$design)
      family$mean_par(subset_init)
    })

    init_ns <- vapply(rownames(grp_init)[-1], FUN.VALUE = numeric(1), FUN = \(ns) {
      mean(link$nuisance[[ns]]$fun(grp_init[ns, ]))
    })

    c(init_ns,
      lm.fit(x = do.call(rbind, lapply(grouped, `[[`, 'x')),
             y = link$mean$fun(grp_init['mean', ]))$coefficients)
  } else if (length(initial) != 2) {
    if (is.null(initial$coef)) {
      abort("`initial` does not contain information about the coefficients")
    }
    if (length(family$nuisance_names) > 0L && is.null(initial$nuisance)) {
      abort(paste("`initial` does not contain information about the nuisance parameters",
                  paste(family$nuisance_names, collapse = ", ")))
    }
  } else if (!is.null(names(initial$coefs))) {
    c(initial$nuisance, initial$coefs[colnames(xmat)])
  } else {
    c(initial$nuisance, initial$coefs)
  }

  names(initial_params) <- NULL
  nuisance_components <- seq_along(family$nuisance_names)
  mean_components <- length(family$nuisance_names) + seq_len(ncol(xmat))

  nuisance_linkinv_apply <- function (params) {
    nuisance <- numeric(length(nuisance_components))
    for (ns in nuisance_components) {
      nuisance[[ns]] <- link$nuisance[[ns]]$inv(params[[ns]])
    }
    nuisance
  }

  if (is.null(optim_control$parscale)) {
    optim_control$parscale <- initial_params
  }

  mpd_est <- optim(initial_params,
                   fn = \(regpars) {
                     nuisance <- nuisance_linkinv_apply(regpars)
                     vapply(grouped, FUN.VALUE = numeric(1), FUN = \(gr) {
                       mean <- drop(gr$x %*% regpars[mean_components]) |>
                         link$mean$inv()
                       # We need to transform here because the `divergence_int()` back-transforms!
                       gr_params <- family$trans(
                         family$parameters_from_mean_par(mean, nuisance))
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

  initial_params[nuisance_components] <- nuisance_linkinv_apply(initial_params)
  names(initial_params) <- c(family$nuisance_names, colnames(xmat))
  reg_coefs <- mpd_est$par[mean_components]
  names(reg_coefs) <- colnames(xmat)
  nuisance_est <- nuisance_linkinv_apply(mpd_est$par)
  names(nuisance_est) <- family$nuisance_names

  # Covariance matrix estimation
  group_cov_info <- lapply(grouped, FUN = \(gr) {
    mean <- drop(gr$x %*% reg_coefs)
    mean_linked <- link$mean$inv(mean)
    link_inv_d1 <- link$mean$inv_d1(mean)
    link_inv_d2 <- link$mean$inv_d2(mean)

    estimates <- family$parameters_from_mean_par(mean_linked, nuisance_est) |>
      setNames(family$parameter_names)
    jac <- family$jacobian_mean_par_mapping(mean_linked, nuisance_est)

    A_est <- gr$sandwich_cov_A(estimates)
    G_est <- gr$sandwich_effective_score(estimates, jac)

    score <- family$raw_scores(gr$y, estimates)
    score_jac <- score %*% jac
    colnames(score_jac) <- c(".mean", family$nuisance_names)

    score_reg <- (link_inv_d1 * outer(score_jac[, 1L], gr$x)) |>
      cbind(score_jac[, -1, drop = FALSE])
    A_est[] <- crossprod(jac, A_est) %*% jac

    Fmat <- Amat <- matrix(NA_real_, nrow = p + d, ncol = p + d)

    finf_jac <- crossprod(jac, family$fisher_inf(estimates)) %*% jac
    Fmat[1:p, 1:p] <- link_inv_d1^2 * finf_jac[[1, 1]] * tcrossprod(gr$x)
    Fmat[(p + 1):(p + d), 1:p] <- Fmat[1:p, (p + 1):(p + d)] <-
      link_inv_d1 * drop(finf_jac[-1L, 1L, drop = FALSE] %*% gr$x)
    Fmat[(p + 1):(p + d), (p + 1):(p + d)] <- finf_jac[-1L, -1L, drop = FALSE]

    Amat[1:p, 1:p] <- (link_inv_d1^2 * A_est[[1, 1]] + link_inv_d2 * G_est) * tcrossprod(gr$x)
    Amat[(p + 1):(p + d), 1:p] <- Amat[1:p, (p + 1):(p + d)] <-
      link_inv_d1 * drop(A_est[-1L, 1L, drop = FALSE] %*% gr$x)
    Amat[(p + 1):(p + d), (p + 1):(p + d)] <- A_est[-1L, -1L, drop = FALSE]

    list(score = score_reg,
         A     = gr$prop * Amat,
         finf  = gr$prop * Fmat)
  })

  finf <- Reduce(`+`, lapply(group_cov_info, `[[`, 'finf'))
  A_est <- Reduce(`+`, lapply(group_cov_info, `[[`, 'A'))

  score_all <- do.call(rbind, lapply(group_cov_info, `[[`, 'score'))
  design_vcov <- vcov(svymean(score_all, design = design))

  neff <- diag(finf) / diag(design_vcov)

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

  covest <- tryCatch({
    (solve(A_est, design_vcov) %*% solve(A_est) * divergence$phi_2nd_deriv_at_1^2) |>
      structure(type = "sandwich")
  }, error = \(cnd) {
    # Compute the model-based estimator instead.
    warn("Sandwich covariance estimate is singular. Computing model-based estimate instead.",
         parent = cnd)

    (finf_inv / outer(sqrt(neff), sqrt(neff))) |>
      structure(type = "model")
  })

  dimnames(covest) <- dimnames(finf_inv) <-
    list(c(names(reg_coefs), paste0('.', names(nuisance_est))),
         c(names(reg_coefs), paste0('.', names(nuisance_est))))

  structure(
    list(coefficients   = reg_coefs,
         nuisance       = nuisance_est,
         cov            = covest,
         neff_scores    = neff,
         nobs           = nrow(xmat),
         neff_kish      = sum(weights(design))^2 / sum(weights(design)^2),
         family         = family,
         link           = link,
         divergence     = divergence,
         finf_inv       = finf_inv,
         mpd            = mpd_est$value,
         initial        = initial_params,
         optimizer_code = mpd_est$convergence,
         optimizer_msg  = mpd_est$message),
    class = 'survey_reg_mde')
}

subset_svydesign <- function (x, i) {
  d <- x[i, ]
  if (!is.null(d$fpc)) {
    d$fpc$sampsize <- d$fpc$sampsize[i, ]
  }
  if (!is.null(d$dcheck)) {
    d$dcheck <- lapply(d$dcheck, \(dc) {
      dc$id <- dc$id[i]
      dc$dcheck <- dc$dcheck[i, i]
      dc
    })
  }
  d
}
