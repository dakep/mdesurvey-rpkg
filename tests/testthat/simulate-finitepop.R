#' Simulate a finite population from a superpopulation model
#'
#' @param size the number of units in the finite population
#' @param cor the correlation between the charateristic $Y$ and the PPS variable $Z$
#' @param ranf a function to generate random deviates from the superpopulation model
#' @return a data frame with 3 columns:
#'   * `y` … the characteristic of interest
#'   * `z` … the variable to base the inclusion weights on
simulate_finitepop <- function (size, cor, ranf) {
  y <- ranf(size)
  if (isTRUE(cor > 0)) {
    sd_z <- sd(y) * sqrt(1 - cor^2) / cor
    xi <- rlnorm(size)
    z <- y + (xi / sd(xi)) * sd_z

    xi <- rlnorm(size)
  } else {
    z <- rlnorm(size)
  }

  data.frame(y = y,
             z = z)
}

#' Simulate a finite population from a superpopulation regression model
#'
#' @param size the number of units in the finite population
#' @param terms a list of terms of the form
#'   `list(intercept, var_name = c(level1 = coef1, level2 = coef2, ...), ...)`
#'   The expected value for each unit is determined by the sum of the coefficients.
#'   If first element in `terms` is a scalara, it treated as the intercept.
#' @param cor the correlation between the charateristic \eqn{Y} and the PPS variable \eqn{Z}
#' @param ranf a function to generate random deviates from the superpopulation model given
#'  the expected value.
#' @return a data frame with columns:
#'   * `y` … the characteristic of interest
#'   * `z` … the variable to base the inclusion weights on
#'   * ...   all the terms from `terms`
simulate_finitepop_lm <- function (size, terms, cor, ranf, link_inv = identity) {
  if (length(terms[[1]]) == 1L) {
    intercept <- terms[[1]]
    terms <- terms[-1]
  }
  df <- lapply(terms, \(t) {
    sample(names(t), size = size, replace = TRUE) |>
      factor(levels = names(t))
  }) |>
    list2DF()

  dfext <- split(df, do.call(interaction, df), drop = TRUE) |>
    unname() |>
    lapply(\(sdf) {
      rownames(sdf) <- NULL
      mean <- vapply(colnames(df), FUN.VALUE = numeric(1), FUN = \(n) {
        unname(terms[[n]][sdf[[n]][[1]]])
      }) |>
        sum()
      mean[[1]] <- link_inv(mean[[1]] + intercept)
      fp <- simulate_finitepop(size = nrow(sdf), cor = cor,
                               ranf = \(n) ranf(n, mean))

      cbind(fp, sdf)
    }) |>
    do.call(what = rbind)
}

stratified_sampling <- function (n, strata, finite_pop, sampling = c('pps', 'srs', 'srswr'),
                                 pps_aux) {
  requireNamespace('survey')
  sampling <- match.arg(sampling)
  pps_aux_vals <- if (!missing(pps_aux)) {
    model.matrix(pps_aux, finite_pop)[, 2]
  } else {
    NULL
  }
  cf <- model.frame(strata, finite_pop)
  strat <- do.call(interaction, cf)

  n_per_stratum <- rep.int(n %/% nlevels(strat), nlevels(strat))
  rem <- n %% nlevels(strat)
  if (isTRUE(rem > 0)) {
    n_per_stratum[seq_len(rem)] <- n_per_stratum[seq_len(rem)] + 1L
    n_per_stratum <- sample(n_per_stratum)
  }

  names(n_per_stratum) <- levels(strat)

  ind <- seq_len(nrow(finite_pop))
  strat_ind <- split(ind, strat)

  sample <- lapply(levels(strat), \(sn) {
    sind <- ind[strat == sn]
    strat_sample_size <- n_per_stratum[[sn]]

    if (identical(sampling, 'pps')) {
      # Sample with Poisson-PPS
      pps_incl_prob <- strat_sample_size * pps_aux_vals[sind] / sum(pps_aux_vals[sind])
      sind[runif(length(sind)) <= pps_incl_prob]
    } else if (identical(sampling, 'srs')) {
      # SRS-WOR
      sample(sind, size = strat_sample_size, replace = FALSE)
    } else if (identical(sampling, 'srswr')) {
      # SRS-WR
      sample(sind, size = strat_sample_size, replace = TRUE)
    }
  }) |>
    unlist(recursive = FALSE, use.names = FALSE)

  design <- if (identical(sampling, 'pps')) {
    # Sample with Poisson-PPS
    pps_incl_prob <- n * pps_aux_vals[sample] / sum(pps_aux_vals)
    if (length(sample) <= sqrt(2^30 / 8)) {
      # If the PPS matrix takes less than 1GB in memory
      pps_info <- outer(pps_incl_prob, pps_incl_prob)
      diag(pps_info) <- pps_incl_prob
      pps_info <- survey::ppsmat(pps_info)
      variance <- 'HT'
    } else {
      # Otherwise, use the Brewer approximation
      pps_info <- 'brewer'
      variance <- 'HT'
    }

    design <- survey::svydesign(ids          = ~ 1,
                                probs        = pps_incl_prob,
                                pps          = pps_info,
                                variance     = variance,
                                strata       = strata,
                                check.strata = FALSE,
                                data         = finite_pop[sample, ])
    # With PPS, the survey package expects a vector for the strata, not a data frame
    design$strata <- design$strata[, 1L, drop = TRUE]
    design
  } else if (identical(sampling, 'srs')) {
    # SRS-WOR
    survey::svydesign(ids    = ~ 1,
                      probs  = NULL,
                      strata = strata,
                      data   = finite_pop[sample, ])
  } else if (identical(sampling, 'srswr')) {
    # SRS-WR
    probs <- 1 - (1 - 1 / as.numeric(table(strat)))^n_per_stratum
    survey::svydesign(ids    = ~ 1,
                      probs  = probs[strat[sample]],
                      strata = strata,
                      check.strata = FALSE,
                      data   = finite_pop[sample, ])
  }
}
