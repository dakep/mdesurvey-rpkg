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
