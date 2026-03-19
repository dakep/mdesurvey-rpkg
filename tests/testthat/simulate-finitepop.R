#' Simulate a finite population from a superpopulation model
#'
#' @param size the number of units in the finite population
#' @param cor the correlation between the charateristic $Y$ and the PPS variable $Z$
#' @param clusters the number of clusters
#' @param ranf a function to generate random deviates from the superpopulation model
#' @return a data frame with 3 columns:
#'   * `y` … the characteristic of interest
#'   * `z` … the variable to base the inclusion weights on
#'   * `x` … a calibration variable
#'   * `cluster` … the cluster assignment
simulate_finitepop <- function (size, cor, clusters = 5, ranf) {
  y <- ranf(size)
  if (isTRUE(cor > 0)) {
    sd_z <- sd(y) * sqrt(1 - cor^2) / cor
    xi <- rlnorm(size)
    z <- y + (xi / sd(xi)) * sd_z
  } else {
    z <- rlnorm(size)
  }

  data.frame(y = y,
             z = z)
}

#' Contaminate a sample from the finite population
#'
#' @param design the sampling frame
#' @param y a character vector with the names of columns in `design$variables` that should
#'   be contaminated.
#' @param type the type of contamination to introduce. Either point contamination at `contparam`
#'   or contamination by a Student t-distribution with `contparam` degrees of freedom.
#' @param qfun the model's quantile function
#' @param corr if `'none'`, no systematic correlation between the auxiliary size variable and
#'   the contamination indicator.
#'   If `corr='pos'` they are positively correlated (i.e., large auxiliary size variable leads
#'   to higher probability of contamination).
#'   If `corr='neg'`, they are anti-correlated.
#' @param prop proportion of contaminated observations
#' @param contparam  contamination parameter.
#'   For `type = "point"`, this is the location of the point contamination in terms
#'   of the quantile from `qfun`.
#'   For `type = "abst"`, this is the degrees of freedom of the (absolute) Student t-distribution.
#'   For `type = "t"`, this is a two element vector with first element the degrees of
#'     freedom of the (absolute) Student t-distribution and the second element the non-centrality
#'     parameter.
#' @return a list with contaminated observations `y` and the indices of the contaminated observations
contaminate_sample <- function (design, y = 'y', type = c("point", "abst", "t"), qfun,
                                contparam, prop, corr = c("none", "pos", "neg")) {
  sample_size <- length(design$prob)
  type <- match.arg(type)
  corr <- match.arg(corr)
  # Determine contamination probability
  cont_prob <- switch(corr,
                      pos  = design$prob^10,
                      neg  = (max(design$prob) - design$prob)^50,
                      none = NULL)

  cont_ind <- sample.int(n = sample_size,
                         size = as.integer(prop * sample_size),
                         prob = cont_prob)

  if (identical(type, "point")) {
    # Point contamination
    # Ensure `contparam` is a probability
    if (length(contparam) != 1 || !isTRUE(contparam[[1]] > 0) || !isTRUE(contparam[[1]] < 1)) {
      stop("`contparam` for `type=\"point\"` must be between 0 and 1")
    }
    cont_mu <- drop(qfun(contparam, cont_ind))
    eps <- 10^floor(log10(1 - contparam))
    cont_sd <- 0.3 * pmin(abs(cont_mu - qfun(pmin(1 - .Machine$double.eps,
                                                  pmax(.Machine$double.eps,
                                                       contparam + c(-eps, eps))),
                                             cont_ind)))

    design$variables[cont_ind, y] <- pmax(.Machine$double.eps,
                                          rnorm(length(cont_ind), mean = cont_mu, sd = cont_sd))
  } else if (identical(type, "abst")) {
    if (length(contparam) != 1 || !isTRUE(contparam[[1]] > 0)) {
      stop("`contparam` for `type=\"abst\"` must be positive")
    }

    model_iqr <- apply(qfun(c(0.25, 0.75), cont_ind), 2, diff)
    t_iqr <- diff(qt(c(0.25, 0.75), df = contparam))
    scaling <- model_iqr / t_iqr

    design$variables[cont_ind, y] <- design$variables[cont_ind, y] +
      abs(scaling * rt(length(cont_ind), df = contparam))
  } else if (identical(type, "t")) {
    if (length(contparam) != 2 || !isTRUE(contparam[[1]] > 0)) {
      stop("`contparam` for `type=\"t\"` must be two numerics, with the first element positive")
    }

    model_iqr <- apply(qfun(c(0.25, 0.75), , cont_ind), 2, diff)
    t_iqr <- diff(qt(c(0.25, 0.75), df = contparam))
    scaling <- model_iqr / t_iqr

    design$variables[cont_ind, y] <- contparam[[2]] +
      scaling * rt(length(cont_ind), df = contparam[[1]])
  }

  design
}

#' Simulate a finite population from a superpopulation regression model
#'
#' @param size the number of units in the finite population
#' @param terms a list of terms of the form
#'   `list(intercept, var_name = c(level1 = coef1, level2 = coef2, ...), ...)`
#'   The expected value for each unit is determined by the inverse link applied
#'   to the sum of the coefficients.
#'   If first element in `terms` is a scalar, it is treated as the intercept.
#'   \deqn{
#'   E[Y | X = x] = \eta^{-1}(x'\beta)
#'   }
#' @param cor the correlation between the characteristic \eqn{Y | X = x} and the
#'   PPS variable \eqn{Z | X = x}
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

determine_neff <- function (incl_prob) {
  list(neff  = length(incl_prob)^2 / sum(1 / incl_prob),
       nveff = length(incl_prob)^2 / sum((1 - incl_prob) / incl_prob))
}

get_design <- function (n, strata = NULL, finite_pop, sampling = c('pps', 'srswor', 'srswr'), pps_aux) {
  requireNamespace('survey')
  sampling <- match.arg(sampling)
  pps_aux_vals <- if (!missing(pps_aux)) {
    model.matrix(pps_aux, finite_pop)[, 2]
  } else {
    NULL
  }

  strat <- if (!is.null(strata)) {
    model.frame(strata, finite_pop) |>
      do.call(what = interaction)
  } else {
    factor(integer(nrow(finite_pop)), levels = 0L, labels = "1")
  }


  n_per_stratum <- rep.int(n %/% nlevels(strat), nlevels(strat))
  rem <- n %% nlevels(strat)
  if (isTRUE(rem > 0)) {
    n_per_stratum[seq_len(rem)] <- n_per_stratum[seq_len(rem)] + 1L
    n_per_stratum <- sample(n_per_stratum)
  }

  names(n_per_stratum) <- levels(strat)

  ind <- seq_len(nrow(finite_pop))
  strat_ind <- split(ind, strat)

  sample_info <- lapply(levels(strat), \(sn) {
    sind <- ind[strat == sn]
    strat_sample_size <- n_per_stratum[[sn]]

    if (identical(sampling, 'pps')) {
      # Sample with Poisson-PPS
      pps_incl_prob <- strat_sample_size * pps_aux_vals[sind] / sum(pps_aux_vals[sind])
      incl <- which(runif(length(sind)) <= pps_incl_prob)
      list(incl_prob = pps_incl_prob[incl],
           all_incl_prob = pps_incl_prob,
           sample    = sind[incl])
    } else if (identical(sampling, 'srswor')) {
      # SRS-WOR
      list(incl_prob = rep.int(strat_sample_size / length(sind), strat_sample_size),
           sample    = sample(sind, size = strat_sample_size, replace = FALSE))

    } else if (identical(sampling, 'srswr')) {
      # SRS-WR
      list(incl_prob = rep.int(1 - (1 - 1 / length(sind))^strat_sample_size, strat_sample_size),
           sample    = sample(sind, size = strat_sample_size, replace = TRUE))
    }
  })

  incl_prob <- lapply(sample_info, `[[`, 'incl_prob') |> unlist(FALSE, FALSE)
  sample <- lapply(sample_info, `[[`, 'sample') |> unlist(FALSE, FALSE)
  if (identical(sampling, 'pps')) {
    # Sample with Poisson-PPS
    if (length(sample) <= sqrt(2^30 / 8)) {
      # If the PPS matrix takes less than 1GB in memory
      pps_info <- outer(incl_prob, incl_prob)
      diag(pps_info) <- incl_prob
      pps_info <- survey::ppsmat(pps_info)
      variance <- 'HT'
    } else {
      # Otherwise, use the Brewer approximation
      pps_info <- 'brewer'
      variance <- 'HT'
    }

    design <- survey::svydesign(ids          = ~ 1,
                                probs        = incl_prob,
                                pps          = pps_info,
                                variance     = variance,
                                strata       = strata,
                                check.strata = FALSE,
                                data         = finite_pop[sample, ])
    design$neff <- lapply(sample_info, `[[`, 'all_incl_prob') |>
      unlist(FALSE, FALSE) |>
      determine_neff()
  } else if (identical(sampling, 'srswor')) {
    # SRS-WOR
    design <- survey::svydesign(ids    = ~ 1,
                                probs  = incl_prob,
                                strata = strata,
                                data   = finite_pop[sample, ])
    design$neff <- determine_neff(rep.int(n / nrow(finite_pop), nrow(finite_pop)))
  } else if (identical(sampling, 'srswr')) {
    # SRS-WR
    design <- survey::svydesign(ids    = ~ 1,
                                probs  = incl_prob,
                                strata = strata,
                                check.strata = FALSE,
                                data   = finite_pop[sample, ])
    probs <- 1 - (1 - 1 / as.numeric(table(strat)))^n_per_stratum
    design$neff <- determine_neff(probs[strat])
  }

  design
}
