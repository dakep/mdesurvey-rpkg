#' Evaluate the Univariate Weighted KDE
#'
#' None of the parameters are checked! They are assumed to be of the right type and size.
#'
#' @param x       observations (numeric vector)
#' @param wgts    weights (numeric vector of same length as `x`)
#' @param evalpts **ordered(!)** evaluations points (numeric vector)
#' @param bw      bandwidth (positive scalar)
#' @param kernel the first letter of the kernel function to use:
#'   e ... Epanechnikov,
#'   b ... Biweight
#'   t ... Triangular
#'   r ... Rectangular
#' @keywords internal
.kde <- function (x, wgts, evalpts, bw, kernel) {
  .Call(C_kde_univariate,
        as.numeric(x),
        as.numeric(wgts),
        as.numeric(evalpts),
        as.numeric(bw),
        as.character(kernel))
}

kde_boundary_correction <- function (evalpts, kde, bw, range, kernel) {
  kernel <- substr(kernel[[1]], 1, 1)
  corr_fun <- switch(kernel,
                     e = \(q) 0.25 * (2 + 3 * q - q^3), # epanechnikov
                     t = \(q) 0.5  * (1 + 2 * q - q^2), # triangular
                     r = \(q) 0.5 * (1 + q), # rectangular
                     b = \(q) 0.5 + 15 * q / 16 - 5 * q^3 / 8 + 3 * q^5 / 16) # biweight

  if (is.finite(range[[1]])) {
    bc <- which(evalpts < range[[1]] + bw)
    if (length(bc) > 0L) {
      kde[bc] <- kde[bc] / corr_fun((evalpts[bc] - range[[1]]) / bw)
    }
  }
  if (is.finite(range[[2]])) {
    bc <- which(evalpts > range[[2]] - bw)
    if (length(bc) > 0L) {
      kde[bc] <- kde[bc] / corr_fun((range[[2]] - evalpts[bc]) / bw)
    }
  }
  kde
}

