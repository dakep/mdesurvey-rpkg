library(testthat)
library(mdesurvey)
library(survey)
source(test_path("simulate-finitepop.R"))

test_that("Gamma MNEDE", {
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

  res <- survey_mpde(~ y,
                     design     = des,
                     divergence = 'ned',
                     family     = 'gamma')

  res <- survey_mpde(~ y,
                     design     = des,
                     divergence = 'HD',
                     family     = 'gamma')

  res_hd <- survey_mhde(~ y,
                        design     = des,
                        family     = 'gamma')

  coef(res)
  coef(res_hd)
  se <- res$cov |> diag() |> sqrt()
  se_hd <- res_hd$cov |> diag() |> sqrt()
  se * 2.6 / se_hd

  expect_equal(coef(res), c(shape = 1.733, scale = 40780), tolerance = 1e-3)
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

  res <- survey_mpde(~ y,
                     design     = des,
                     divergence = 'ned',
                     family     = 'normal')

  expect_equal(coef(res), c(mean = 2.117, sd = 1.3), tolerance = 1e-3)
})
