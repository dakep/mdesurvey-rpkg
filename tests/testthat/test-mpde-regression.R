library(testthat)
library(mdesurvey)
library(survey)
source(test_path("simulate-finitepop.R"))

test_that("Gamma Regression", {
  N <- 1e6
  n <- 1e3
  gamma_var <- 2 * 35000^2

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
      params <- gamma_family$reparameterize_from_mv(mean, var = gamma_var)
      rgamma(n, shape = params[['shape']], scale = params[['scale']])
    })

  design <- stratified_sampling(n = n, strata = ~ x1, finite_pop = finite_pop, sampling = 'pps',
                                pps_aux = ~ z)

  res <- survey_regression_mpde(y ~ x1 + x2 + x3,
                                design     = design,
                                divergence = 'ned',
                                family     = gamma_family,
                                optim_control = list(maxit = 5001))


})

test_that("Normal Regression", {
  N <- 1e6
  n <- 1e3
  var <- 1e8

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
      params <- normal_family$reparameterize_from_mv(mean, var = var)
      rnorm(n, mean = params[['mean']], sd = params[['sd']])
    })

  design <- stratified_sampling(n = n, strata = ~ x1, finite_pop = finite_pop, sampling = 'pps',
                                pps_aux = ~ z)

  res <- survey_regression_mpde(y ~ x1 + x2 + x3,
                                design     = design,
                                divergence = 'ned',
                                family     = normal_family,
                                optim_control = list(maxit = 5001))
})
