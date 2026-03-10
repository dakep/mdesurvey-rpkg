library(testthat)
library(mdesurvey)
library(survey)
source(test_path("simulate-finitepop.R"))

test_that("Gamma Regression (identity link)", {
  N <- 1e6
  n <- 1e3
  gamma_sd <- sqrt(2) * 35000

  gamma_family <- model_family('gamma')
  set.seed(1)
  finite_pop <- simulate_finitepop_lm(
    size = N,
    terms = list(Intercept = 50000,
                 x1 = c(a = 0, b = 40000, c = 50000),
                 x2 = c(a = 0, b = 30000),
                 x3 = c(a = 0, b = 15000, c = 40000)),
    cor = 0.8,
    ranf = \(n, mean) {
      params <- gamma_family$parameters_from_mean_par(mean, nuisance = gamma_sd)
      rgamma(n, shape = params[['shape']], scale = params[['scale']])
    })

  design <- stratified_sampling(n = n, strata = ~ x1, finite_pop = finite_pop, sampling = 'pps',
                                pps_aux = ~ z)

  res <- expect_no_error(
    survey_regression_mpde(y ~ x1 + x2 + x3,
                           design     = design,
                           link       = link('identity', sd = 'log'),
                           divergence = 'ned',
                           family     = gamma_family,
                           optim_control = list(maxit = 5001)))
  expect_length(coef(res), 6)
  expect_length(sigma(res), 1)
  vcov_sandwich <- expect_no_error(vcov(res))
  expect_shape(vcov_sandwich, dim = c(6L, 6L))
  expect_all_true(sqrt(diag(vcov_sandwich)) > 1e3)
  expect_no_error(vcov(res, which = 'nuisance')) |>
    expect_shape(dim = c(1L, 1L))
  expect_no_error(vcov(res, which = 'all')) |>
    expect_shape(dim = c(7L, 7L))

  vcov_model <- expect_no_error(vcov(res, type = "model"))
  expect_shape(vcov_model, dim = c(6L, 6L))
  expect_all_true(sqrt(diag(vcov_model)) > 1e3)
  expect_no_error(vcov(res, which = 'nuisance', type = "model")) |>
    expect_shape(dim = c(1L, 1L))
  expect_no_error(vcov(res, which = 'all', type = "model")) |>
    expect_shape(dim = c(7L, 7L))
})

test_that("Gamma Regression (log link)", {
  N <- 1e6
  n <- 1e3
  gamma_sd <- sqrt(2) * 35000

  gamma_family <- model_family('gamma')
  set.seed(1)
  finite_pop <- simulate_finitepop_lm(
    size = N,
    terms = list(Intercept = log(50000),
                 x1 = c(a = 0, b = log(1.8), c = log(2)),
                 x2 = c(a = 0, b = log(1.5)),
                 x3 = c(a = 0, b = log(1.2), c = log(1.4))),
    cor = 0.8,
    link_inv = exp,
    ranf = \(n, mean) {
      params <- gamma_family$parameters_from_mean_par(mean, nuisance = gamma_sd)
      rgamma(n, shape = params[['shape']], scale = params[['scale']])
    })

  design <- stratified_sampling(n = n, strata = ~ x1, finite_pop = finite_pop, sampling = 'pps',
                                pps_aux = ~ z)

  res <- expect_no_error(
    survey_regression_mpde(y ~ x1 + x2 + x3,
                           design     = design,
                           link       = link('log', sd = 'log'),
                           divergence = 'ned',
                           family     = gamma_family,
                           optim_control = list(maxit = 5001)))
  expect_length(coef(res), 6)
  expect_length(sigma(res), 1)
  vcov_sandwich <- expect_no_error(vcov(res))
  expect_shape(vcov_sandwich, dim = c(6L, 6L))
  expect_all_true(sqrt(diag(vcov_sandwich)) < 0.1)
  expect_no_error(vcov(res, which = 'nuisance')) |>
    expect_shape(dim = c(1L, 1L))
  expect_no_error(vcov(res, which = 'all')) |>
    expect_shape(dim = c(7L, 7L))

  vcov_model <- expect_no_error(vcov(res, type = "model"))
  expect_shape(vcov_model, dim = c(6L, 6L))
  expect_all_true(sqrt(diag(vcov_model)) < 0.1)
  expect_no_error(vcov(res, which = 'nuisance', type = "model")) |>
    expect_shape(dim = c(1L, 1L))
  expect_no_error(vcov(res, which = 'all', type = "model")) |>
    expect_shape(dim = c(7L, 7L))

})

test_that("Normal Regression", {
  N <- 1e6
  n <- 1e3
  sd <- 1e4

  normal_family <- model_family('normal')
  set.seed(1)
  finite_pop <- simulate_finitepop_lm(
    size = N,
    terms = list(Intercept = 50000,
                 x1 = c(a = 0, b = 40000, c = 50000),
                 x2 = c(a = 0, b = 30000),
                 x3 = c(a = 0, b = 15000, c = 40000)),
    cor = 0.8,
    ranf = \(n, mean) {
      params <- normal_family$parameters_from_mean_par(mean, sd)
      rnorm(n, mean = params[['mean']], sd = params[['sd']])
    })

  design <- stratified_sampling(n = n, strata = ~ x1, finite_pop = finite_pop, sampling = 'pps',
                                pps_aux = ~ z)

  res <- expect_no_error(
    survey_regression_mpde(y ~ x1 + x2 + x3,
                           design     = design,
                           link       = link('identity', sd = 'log'),
                           divergence = 'ned',
                           family     = normal_family,
                           optim_control = list(maxit = 5001)))

  expect_length(coef(res), 6)
  expect_length(sigma(res), 1)
  vcov_sandwich <- expect_no_error(vcov(res))
  expect_shape(vcov_sandwich, dim = c(6L, 6L))
  expect_all_true(sqrt(diag(vcov_sandwich)) > 500)
  expect_no_error(vcov(res, which = 'nuisance')) |>
    expect_shape(dim = c(1L, 1L))
  expect_no_error(vcov(res, which = 'all')) |>
    expect_shape(dim = c(7L, 7L))

  vcov_model <- expect_no_error(vcov(res, type = "model"))
  expect_shape(vcov_model, dim = c(6L, 6L))
  expect_all_true(sqrt(diag(vcov_model)) > 500)
  expect_no_error(vcov(res, which = 'nuisance', type = "model")) |>
    expect_shape(dim = c(1L, 1L))
  expect_no_error(vcov(res, which = 'all', type = "model")) |>
    expect_shape(dim = c(7L, 7L))
})

test_that("Student-t Regression", {
  N <- 1e6
  n <- 1e3
  scale <- 5e3
  df <- 4

  t_family <- model_family('t')
  set.seed(1)
  finite_pop <- simulate_finitepop_lm(
    size = N,
    terms = list(Intercept = 50000,
                 x1 = c(a = 0, b = 40000, c = 50000),
                 x2 = c(a = 0, b = 30000),
                 x3 = c(a = 0, b = 15000, c = 40000)),
    cor = 0.8,
    ranf = \(n, mean) {
      params <- t_family$parameters_from_mean_par(mean, nuisance = c(scale = scale, df = df))
      x <- rt(n, df = params[['df']])
      params[['location']] + x * params[['scale']]
    })

  design <- stratified_sampling(n = n, strata = ~ x1, finite_pop = finite_pop, sampling = 'pps',
                                pps_aux = ~ z)

  res <- expect_no_error(
    survey_regression_mpde(y ~ x1 + x2 + x3,
                           design        = design,
                           divergence    = 'ned',
                           family        = t_family,
                           optim_method  = 'BFGS',
                           optim_control = list(maxit = 5001)))

  expect_length(coef(res), 6)
  expect_length(sigma(res), 2)
  vcov_sandwich <- expect_no_error(vcov(res))
  expect_shape(vcov_sandwich, dim = c(6L, 6L))
  expect_all_true(sqrt(diag(vcov_sandwich)) > 100)
  expect_no_error(vcov(res, which = 'nuisance')) |>
    expect_shape(dim = c(2L, 2L))
  expect_no_error(vcov(res, which = 'all')) |>
    expect_shape(dim = c(8L, 8L))

  vcov_model <- expect_no_error(vcov(res, type = "model"))
  expect_shape(vcov_model, dim = c(6L, 6L))
  expect_all_true(sqrt(diag(vcov_model)) > 100)
  expect_no_error(vcov(res, which = 'nuisance', type = "model")) |>
    expect_shape(dim = c(2L, 2L))
  expect_no_error(vcov(res, which = 'all', type = "model")) |>
    expect_shape(dim = c(8L, 8L))
})
