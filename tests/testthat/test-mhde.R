library(testthat)
library(mdesurvey)
library(survey)
source(test_path("simulate-finitepop.R"))

test_that("Gamma MHDE", {
  TRUE_SHAPE <- 2
  TRUE_SCALE <- 35000
  N <- 1e6
  n <- 1e3

  set.seed(1)
  finite_pop <- simulate_finitepop(size = N, cor = 0.8,
                                   ranf = \(n) rgamma(n, shape = TRUE_SHAPE, scale = TRUE_SCALE))

  finite_pop$incl_prob <- n * finite_pop$z / sum(finite_pop$z)

  set.seed(1)
  sample <- which(runif(N) <= finite_pop$incl_prob)

  des <- svydesign(id = ~ 1,
                   variables = ~ y,
                   probs = finite_pop$incl_prob[sample],
                   data = data.frame(y = finite_pop$y[sample]),
                   pps = HR(sum(finite_pop$incl_prob^2)))

  res <- survey_mhde(~ y,
                     design = des,
                     family = model_family('gamma'))

  expect_equal(coef(res), c(shape = 1.752, scale = 40063), tolerance = 1e-3)
  expect_equal(res$neff_scores, c(shape = 319, scale = 713), tolerance = 1e-3)
  expect_equal(res$family$fisher_inf(coef(res)),
               matrix(c(0.763077600534, 0.000024960726, 0.000024960726, 0.000000001091), ncol = 2),
               tolerance = 1e-3)
  expect_equal(vcov(res, "sandwich"),
               matrix(c(0.02229989, -383.4634, -383.4634, 7247553.0253), ncol = 2,
                      dimnames = list(c("shape", "scale"), c("shape", "scale"))) |>
                 structure(type = "sandwich"),
               tolerance = 1e-3)
  expect_equal(vcov(res, "model"),
               matrix(c(0.01629254, -249.2733, -249.2733, 5098261.0461), ncol = 2,
                      dimnames = list(c("shape", "scale"), c("shape", "scale"))),
               tolerance = 1e-3)
})

test_that("Normal MHDE", {
  TRUE_MEAN <- 2
  TRUE_SD <- 1.4
  N <- 1e6
  n <- 1e3

  set.seed(1)
  finite_pop <- simulate_finitepop(size = N, cor = 0.8,
                                   ranf = \(n) rnorm(n, mean = TRUE_MEAN, sd = TRUE_SD))

  finite_pop$incl_prob <- n * finite_pop$z / sum(finite_pop$z)

  set.seed(1)
  sample <- which(runif(N) <= finite_pop$incl_prob)

  des <- svydesign(id = ~ 1,
                   variables = ~ y,
                   probs = finite_pop$incl_prob[sample],
                   data = data.frame(y = finite_pop$y[sample]),
                   pps = HR(sum(finite_pop$incl_prob^2)))

  res <- survey_mhde(~ y,
                     design = des,
                     family = model_family('normal'))

  expect_equal(coef(res), c(mean = 2.162, sd = 1.255), tolerance = 1e-3)
  expect_equal(res$neff_scores, c(mean = 386, sd = 575), tolerance = 1e-3)
  expect_equal(vcov(res, "sandwich"),
               matrix(c(0.0038827681, -0.0007518841, -0.0007518841, 0.0010299188), ncol = 2,
                      dimnames = list(c("mean", "sd"), c("mean", "sd"))) |>
                 structure(type = "sandwich"),
               tolerance = 1e-3)
  expect_equal(vcov(res, "model"),
               matrix(c(0.004078841, 0, 0, 0.001370791), ncol = 2,
                      dimnames = list(c("mean", "sd"), c("mean", "sd"))),
               tolerance = 1e-3)
})
