#' Epanechnikov Kernel for KDE
#'
#' Evaluate the Epanechnikov kernel
#' \deqn{
#'   k(x) = \frac{3}{4} (1 - x^2) 1_{[|x| \leq 1]}
#' }
#'
#' @param x numeric vector
#' @return the evaluated kernel
epanechnikov_kernel <- function (x) {
  0.75 * (1 - x^2) * (abs(x) <= 1)
}
