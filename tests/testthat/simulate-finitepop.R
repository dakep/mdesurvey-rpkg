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

    xi <- rlnorm(size)
    x <- y + (xi / sd(xi)) * sd_z
  } else {
    z <- rlnorm(size)
    x <- rlnorm(size)
  }

  cluster <- rep(factor(seq_len(clusters), levels = seq_len(clusters)),
                 length.out = size)

  cluster_seq <- rep(seq_len(ceiling(size / clusters)), each = clusters, length.out = size)

  data.frame(y = y,
             z = z,
             x = x,
             cluster = cluster,
             cluster_seq = cluster_seq)
}
