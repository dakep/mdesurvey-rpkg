#' Minimum Phi Distance Estimator for the General Continuous Models with Complex
#' Survey Designs
#'
#' Compute the survey MPDE for the parameters of a user-specified continuous model.
#'
#' @param x a formula, symbol, expression, vector, or matrix specifying the observations.
#'   The objects are first looked up in the provided `design`, then in the calling environment.
#' @param design the survey design created by [survey::svydesign()] and friends.
#' @param initial an initial estimate for the parameters.
#'   If the vector has names, those will be used for the returned estimates.
#' @param family either the name of a known model family, or a
#'   model family as returned by [model_family()].
#' @param divergence either the name of a known phi divergence, or a
#'   phi divergence as returned by [phi_divergence()].
#' @param phi either the name of a built-in divergence measure or a convex function.
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
#' @family Minimum Phi-Divergence Estimator
#' @export
survey_mpde <- function (x, design,
                         initial,
                         family,
                         divergence = 'HD',
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

  if (is.character(divergence)) {
    divergence <- phi_divergence(divergence)
  }

  if (!is.R6(divergence) || !is(divergence, "PhiDivergence")) {
    abort("`divergence` must be a valid phi divergence from phi_divergence()")
  }

  if (missing(initial)) {
    initial <- family$initial(svy$x, design)
  }
  names(initial) <- family$parameter_names

  if (missing(bw) || is.null(bw)) {
    bw <- bw.nrd(svy$x)
  }

  mpde_integral <- phidiv_gauss_quadrature(svy$x,
                                           wgts = wgts,
                                           divergence = divergence,
                                           family = family,
                                           bandwidth = bw,
                                           kernel = kernel,
                                           n_subdivisions = integration_subdivisions)

  mpd_est <- optim(family$trans(initial),
                   fn = mpde_integral$divergence_int,
                   method = optim_method,
                   control = optim_control)

  if (mpd_est$convergence != 0) {
    warn(sprintf("Optimizer did not converge (code %d): %s",
                 mpd_est$convergence, paste0('', mpd_est$message)))
  }

  estimates <- family$inv_trans(mpd_est$par)
  names(estimates) <- family$parameter_names
  covest <- mpde_covest(svy$x, estimates, mpde_integral, family, design)

  structure(
    list(estimates      = estimates,
         family         = family,
         divergence     = divergence,
         cov            = covest$cov,
         mpd            = mpd_est$value,
         initial        = initial,
         nobs           = nobs,
         neff_kish      = nobs_eff,
         neff_scores    = covest$neff,
         optimizer_code = mpd_est$convergence,
         optimizer_msg  = mpd_est$message),
    class = c('survey_mpde', 'survey_mde'))
}

#' @importFrom survey svymean
#' @importFrom stats vcov
#' @importFrom rlang warn
mpde_covest <- function (x, estimates, integrator, family, design) {
  A_est <- integrator$sandwich_cov_A(estimates)
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

#' Functor to Create Fixed-grid Gaussian Quadrature For General Phi Divergences
#'
#' This creates two function.
#' The function `divergence_int` to computes integral (1):
#' \deqn{\int g_\theta(x) \phi(\hat{f}(x) / g_\theta(x)) \;\mathrm{d}x,}
#' for model density \eqn{g_\theta()} and convex function \eqn{\phi(x)}.
#' The function is set up such that \eqn{\hat f} is evaluated a fixed number of times in each
#' subdivision, and does not need to be re-evaluated for different \eqn{g_\theta()}.
#'
#' The function `sandwich_cov_A` to computes the integrals (2):
#' \deqn{\int g_\theta(x) \left[
#'    w_1(\hat{f}(x) / g_\theta(x)) H_\theta(x) +
#'    w_2(\hat{f}(x) / g_\theta(x)) s_\theta(x) s_\theta(x)^\intercal \right] \;\mathrm{d}x,
#'  }
#' for a model family with density \eqn{g_\theta()}, score function \eqn{s_\theta()}
#' and Hessian function \eqn{H_\theta()}.
#' The functions \eqn{w_j(x)} are associated with the phi divergence and defined as
#' \eqn{w_1(x) = \phi(x) - x \phi'(x)} and \eqn{w_2(x) = x^2 \phi''(x)}.
#' The function is set up such that \eqn{\hat f} is evaluated a fixed number of times in each
#' subdivision, and does not need to be re-evaluated for different members of the same model
#' family.
#'
#'
#' @returns a list of two functions:
#'  `divergence_int` computes the above integral (1), and
#'  `sandwich_cov_A` which computes the integral (2)
#'
#' @keywords internal
#' @importFrom SparseGrid createSparseGrid
#' @rdname phidiv_gauss_quadrature
phidiv_gauss_quadrature <- function (x, wgts, bandwidth, n_subdivisions = 256, poly_order = 20,
                                     divergence, family, kernel) {
  kernel <- substr(kernel[[1]], 1, 1)

  # Evaluate the KDE at the Gaussian evaluation points
  # We don't need to go beyond the range ±bandwidth because the KDE
  # will be 0 outside.
  range_family <- family$range[1, , drop = TRUE]
  range_x <- range(x) + c(-bandwidth, bandwidth)
  range_x[[1]] <- max(range_x[[1]], range_family[[1]])
  range_x[[2]] <- min(range_x[[2]], range_family[[2]])

  subintervals <- partition_domain(x, subdivisions = n_subdivisions, bw = bandwidth,
                                   range = range_x)

  # Determine the gaps in the partitioning.
  # In these gaps, the KDE is 0 and the contribution to the divergence
  # is given by phifun(0) times the probability mass from g() inside these gaps.
  gaps <- cbind(
    c(range_family[[1]], vapply(subintervals, FUN.VALUE = numeric(1), FUN = \(.x) .x[[2]])),
    c(vapply(subintervals, FUN.VALUE = numeric(1), FUN = \(.x) .x[[1]]), range_family[[2]]))
  gaps <- gaps[gaps[, 2] - gaps[, 1] > .Machine$double.eps, , drop = FALSE] |>
    t()

  gap_const <- divergence$phi(0)

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
                 f_hat = .kde(x       = x,
                              wgts    = wgts,
                              evalpts = gkpts[, "int_x"],
                              bw      = bandwidth,
                              kernel  = kernel))

  list(
    # divergence_int = \(dfun, pfun) {
    #   f_theta <- dfun(gkpts[, "int_x"])
    #   if (is.null(dim(f_theta))) {
    #     f_theta <- matrix(f_theta, ncol = 1)
    #   }
    #
    #   integrand <- phifun(gkpts[, "f_hat"] / f_theta) * f_theta
    #   integrand[f_theta < .Machine$double.eps] <- 0
    #   int_val <- crossprod(gkpts[, "int_wgt"], integrand) |>
    #     drop()
    #
    #   gap_contrib <- vapply(seq_len(ncol(gaps)), FUN.VALUE = numeric(1), FUN = \(i) {
    #     gap_const * diff(pfun(gaps[, i]))
    #   }) |>
    #     sum()
    #
    #   int_val + gap_contrib
    # },

    divergence_int = \(params) {
      params[] <- family$inv_trans(params)
      f_theta <- family$dfun(gkpts[, "int_x"], params = params, log = FALSE)
      integrand <- divergence$phi(gkpts[, "f_hat"] / f_theta) * f_theta
      integrand[f_theta < .Machine$double.eps] <- 0
      int_val <- crossprod(gkpts[, "int_wgt"], integrand) |>
        drop()

      gap_contrib <- vapply(seq_len(ncol(gaps)), FUN.VALUE = numeric(1), FUN = \(i) {
        gap_const * diff(family$pfun(gaps[, i], params = params))
      }) |>
        sum()

      int_val + gap_contrib
    },

    sandwich_cov_A = \(params, rel_tol = 1e-12) {
      cov_entries <- matrix(rep(seq_along(params), each = 2L), nrow = 2L) |>
        cbind(combn(length(params), 2L))

      # First we'll evaluate the integral over the regions where both
      # the KDE and the model density are positive.
      scores <- family$raw_scores(gkpts[, "int_x"], params = params)
      hessian <- family$hessian(gkpts[, "int_x"], params = params)
      f_theta <- family$dfun(gkpts[, "int_x"], params = params, log = FALSE)
      ratio <- gkpts[, "f_hat"] / f_theta
      ratio[f_theta < .Machine$double.eps] <- 0

      A_int_pos_kde <- seq_len(ncol(cov_entries)) |>
        vapply(FUN.VALUE = scores[, 1], FUN = \(ceind) {
          ij <- cov_entries[, ceind]
          f_theta * (divergence$w1(ratio) * hessian[, ceind] +
                       divergence$w2(ratio) * scores[, ij[[1]]] * scores[, ij[[2]]])
        })

      A_int_pos_kde <- crossprod(gkpts[, "int_wgt"], A_int_pos_kde) |>
        drop()

      # Then we'll evaluate the integral over the regions where the KDE is zero
      # but the model density is positive.
      w10 <- divergence$w1(0)
      w20 <- divergence$w2(0)
      A_int_gap <- seq_len(ncol(gaps)) |>
        vapply(FUN.VALUE = numeric(ncol(cov_entries)), FUN = \(gapi) {
          seq_len(ncol(cov_entries)) |>
            vapply(FUN.VALUE = numeric(1), FUN = \(ceind) {
              ij <- cov_entries[, ceind]
              if (isTRUE(abs(w10) > 0) && isTRUE(abs(w20) > 0)) {
                integrate(lower = gaps[[1, gapi]],
                          upper = gaps[[2, gapi]],
                          subdivisions = n_subdivisions,
                          f = \(xint) {
                            scores <- family$raw_scores(xint, params = params)
                            hessian <- family$hessian(xint, params = params)[, ceind]
                            f_theta <- family$dfun(xint, params = params, log = FALSE)
                            f_theta * (w10 * hessian +  w20 * scores[, ij[[1]]] * scores[, ij[[2]]])
                          },
                          rel.tol = rel_tol)$value
              } else if (isTRUE(abs(w10) > 0)) {
                integrate(lower = gaps[[1, gapi]],
                          upper = gaps[[2, gapi]],
                          subdivisions = n_subdivisions,
                          f = \(xint) {
                            hessian <- family$hessian(xint, params = params)[, ceind]
                            f_theta <- family$dfun(xint, params = params, log = FALSE)
                            f_theta * w10 * hessian
                          },
                          stop.on.error = FALSE,
                          rel.tol = rel_tol)$value
              } else if (isTRUE(abs(w20) > 0)) {
                integrate(lower = gaps[[1, gapi]],
                          upper = gaps[[2, gapi]],
                          subdivisions = n_subdivisions,
                          f = \(xint) {
                            scores <- family$raw_scores(xint, params = params)
                            f_theta <- family$dfun(xint, params = params, log = FALSE)
                            f_theta * w20 * scores[, ij[[1]]] * scores[, ij[[2]]]
                          },
                          abs.tol = rel_tol)$value
              } else {
                0
              }
            })
        }) |>
        rowSums()


      A_est <- matrix(data = NA_real_, nrow = length(params), ncol = length(params))
      for (ceind in seq_len(ncol(cov_entries))) {
        A_est[[cov_entries[1, ceind], cov_entries[2, ceind]]] <-
          A_est[[cov_entries[2, ceind], cov_entries[1, ceind]]] <-
          A_int_pos_kde[[ceind]] + A_int_gap[[ceind]]
      }

      A_est
    }
  )
}
