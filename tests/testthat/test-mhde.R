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
                     family = model_family('gamma'),
                     optim_method = 'BFGS')

  expect_equal(coef(res), c(shape = 1.6, scale = 43380), tolerance = 1e-3)
  expect_equal(res$neff_scores, c(shape = 363.6, scale = 756.8), tolerance = 1e-3)
  expect_equal(res$family$fisher_inf(coef(res)),
               matrix(c(0.8690425, 0.000002305101, 0.000002305101, 8.422395e-10), ncol = 2),
               tolerance = 1e-3)
  expect_equal(vcov(res, "sandwich"),
               matrix(c(0.01590053, -314.3263, -314.3263, 6959807.9486), ncol = 2,
                      dimnames = list(c("shape", "scale"), c("shape", "scale"))) |>
                 structure(type = "sandwich"),
               tolerance = 1e-3)
  expect_equal(vcov(res, "model"),
               matrix(c(0.01154757, -219.0616, -219.0616, 5724513.7698), ncol = 2,
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

test_that("MHDE vs MPDE-HD in Gamma models", {
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

  mhde <- survey_mhde(~ y,
                      design = des,
                      family = 'gamma')
  mpde <- survey_mpde(~ y,
                      design = des,
                      divergence = 'hd',
                      family = 'gamma')

  expect_equal(coef(mhde), coef(mpde), tolerance = 1e-3)
  expect_equal(vcov(mhde, "sandwich"), vcov(mpde, "sandwich"), tolerance = 1e-3)
  expect_equal(vcov(mhde, "model"), vcov(mpde, "model"), tolerance = 1e-3)
})

test_that("MHDE vs MPDE-HD in Normal models", {
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

  mhde <- survey_mhde(~ y,
                      design = des,
                      family = 'normal')
  mpde <- survey_mpde(~ y,
                      design = des,
                      divergence = 'hd',
                      family = 'normal')

  expect_equal(coef(mhde), coef(mpde), tolerance = 1e-3)
  expect_equal(vcov(mhde, "sandwich"), vcov(mpde, "sandwich"), tolerance = 1e-3)
  expect_equal(vcov(mhde, "model"), vcov(mpde, "model"), tolerance = 1e-3)
})

