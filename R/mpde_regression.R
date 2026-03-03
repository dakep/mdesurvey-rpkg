#' Minimum Phi Distance Estimator for the Linear Regression Models with Complex
#' Survey Designs
#'
#' Compute the survey MPDE for the mean regression parameters under a user-specified continuous
#' model for the conditional distribution \eqn{Y \mid X = x}.
#'
#' @param x a formula, symbol, expression, vector, or matrix specifying the linear model.
#'   The objects are first looked up in the provided `design`, then in the calling environment.
#' @param design the survey design created by [survey::svydesign()] and friends.
#' @param initial an initial estimate for the parameters.
#'   If the vector has names, those will be used to match the elements of `x`.
#' @param family either the name of a known model family, or a
#'   model family as returned by [model_family()].
#' @param divergence either the name of a known phi divergence, or a
#'   phi divergence as returned by [phi_divergence()].
#' @param na.rm a logical evaluating to `TRUE` or `FALSE` indicating whether `NA` values
#'   should be omitted.
#' @param integration_subdivisions number of partitions to divide the domain of \eqn{\hat f()}
#'   for Gauss-Kronrod quadrature.
#' @param bw bandwidth for the HT-adjusted KDE. Uses [stats::bw.nrd()] as default.
#' @param kernel Name of the kernel used in the HT-adjusted KDE.
#'
#' @param optim_method,optim_control method and control options passed on to [stats::optim()].
#' @importFrom stats uniroot optim bw.nrd model.frame model.matrix model.response
#' @importFrom rlang warn abort enquo
#' @importFrom survey svydesign
#' @importFrom methods is
#' @family Minimum Phi-Divergence Estimator
#' @export
survey_regression_mpde <- function (x, design,
                                    initial,
                                    family,
                                    divergence = 'HD',
                                    na.rm = FALSE,
                                    integration_subdivisions = 256,
                                    bw,
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
  mf <- model.frame(x, design$variables)
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

  grouped <- split(seq_along(y), grp) |>
    lapply(\(i) {
      list(y    = y[i],
           wgts = weights(design)[i],
           x    = xmat[i[[1]], ])
    })


  bwf <- if (missing(bw) || is.null(bw)) {
    bw.nrd
  } else if (!is.numeric(bw)) {
    match.fun(bw)
  } else {
    \(...) {
      as.numeric(bw[[1]])
    }
  }

  mpd_ints <- lapply(grouped, \(subset) {
    bw <- bwf(subset$y)
    c(phidiv_gauss_quadrature(subset$y,
                              wgts = subset$wgts,
                              divergence = divergence,
                              family = family,
                              bandwidth = bw,
                              kernel = kernel,
                              n_subdivisions = integration_subdivisions),
      list(x = subset$x))
  })

  if (missing(initial)) {
    grp_init <- vapply(grouped, FUN.VALUE = numeric(2), FUN = \(subset) {
      subset_init <- family$initial(subset$y, svydesign(~ 1, weights = ~ wgts,
                                                        data = data.frame(wgts = subset$wgts)))
      family$reparameterize_to_mv(subset_init)
    })

    initial <- c(lm.fit(x = do.call(rbind, lapply(grouped, `[[`, 'x')),
                        y = grp_init['mean', ])$coefficients,
                 .var = mean(grp_init['var', ]))
  }

  if (length(initial['.var']) != 1) {
    abort("Initial estimate must contain entry for the variance under `.var`")
  }

  initial <- c(.var = log(initial[['.var']]),
               initial[names(initial) != '.var'])

  opt_res <- optim(initial,
                   fn = \(regpars) {
                     vapply(mpd_ints, FUN.VALUE = numeric(1), FUN = \(gr) {
                       mean <- drop(gr$x %*% regpars[-1])
                       gr_params <- family$reparameterize_from_mv(mean, var = exp(regpars[[1]])) |>
                         family$trans()
                       browser()
                       if (all(is.finite(gr_params))) {
                         gr$divergence_int(gr_params)
                       } else {
                         Inf
                       }
                     }) |>
                       sum()
                   },
                   method = optim_method,
                   control = optim_control)

}
