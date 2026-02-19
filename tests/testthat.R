if (require(testthat)) {
  library(mdesurvey)
  test_check("mdesurvey")
} else {
  warning("'mdesurvey' requires 'testthat' for tests.")
}
