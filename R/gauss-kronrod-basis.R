# Matrix of Gauss-Kronrod evaluation points and weights from
# https://www.advanpix.com/2011/11/07/gauss-kronrod-quadrature-nodes-weights/
.kronrod_basis <- matrix(c(
  -9.956571630258080896069827758765314e-01, 1.169463886737187423292549937059448e-02, 0.000000000000000000000000000000000e+00,
  -9.739065285171717434309357486199588e-01, 3.255816230796472476871628032313311e-02, 6.667134430868813799175853773704148e-02,
  -9.301574913557082435744405302102678e-01, 5.475589657435199486545940317228087e-02, 0.000000000000000000000000000000000e+00,
  -8.650633666889845363456856830453034e-01, 7.503967481091995683772921665877220e-02, 1.494513491505805868886369580650353e-01,
  -7.808177265864169047659970601671375e-01, 9.312545458369760054129216086948873e-02, 0.000000000000000000000000000000000e+00,
  -6.794095682990244355892173189204186e-01, 1.093871588022976432119648393381794e-01, 2.190863625159820415877476307286997e-01,
  -5.627571346686046638296829769387841e-01, 1.234919762620658445495536170710693e-01, 0.000000000000000000000000000000000e+00,
  -4.333953941292472133994806426926516e-01, 1.347092173114733393290975982381497e-01, 2.692667193099963496294435572053771e-01,
  -2.943928627014602006362053998600459e-01, 1.427759385770600852882949993727379e-01, 0.000000000000000000000000000000000e+00,
  -1.488743389816312157059030596428784e-01, 1.477391049013384860533193432274857e-01, 2.955242247147528700246255084493896e-01,
  0.000000000000000000000000000000000e+00, 1.494455540029168971738471327626030e-01, 0.000000000000000000000000000000000e+00,
  1.488743389816312157059030596428784e-01, 1.477391049013384860533193432274857e-01, 2.955242247147528700246255084493896e-01,
  2.943928627014602006362053998600459e-01, 1.427759385770600852882949993727379e-01, 0.000000000000000000000000000000000e+00,
  4.333953941292472133994806426926516e-01, 1.347092173114733393290975982381497e-01, 2.692667193099963496294435572053771e-01,
  5.627571346686046638296829769387841e-01, 1.234919762620658445495536170710693e-01, 0.000000000000000000000000000000000e+00,
  6.794095682990244355892173189204186e-01, 1.093871588022976432119648393381794e-01, 2.190863625159820415877476307286997e-01,
  7.808177265864169047659970601671375e-01, 9.312545458369760054129216086948873e-02, 0.000000000000000000000000000000000e+00,
  8.650633666889845363456856830453034e-01, 7.503967481091995683772921665877220e-02, 1.494513491505805868886369580650353e-01,
  9.301574913557082435744405302102678e-01, 5.475589657435199486545940317228087e-02, 0.000000000000000000000000000000000e+00,
  9.739065285171717434309357486199588e-01, 3.255816230796472476871628032313311e-02, 6.667134430868813799175853773704148e-02,
  9.956571630258080896069827758765314e-01, 1.169463886737187423292549937059448e-02, 0.000000000000000000000000000000000e+00),
  ncol = 3,
  byrow = TRUE)

#' Partition the Domain Into Subdivisions
#'
#' This partitions the domain spanned by `x` into a fixed number of partitions such that the
#' KDE \eqn{\hat{f}(x) > 0} within each of those partitions and \eqn{\hat{f}(x) = 0} outside
#' those partitions.
#' The algorithm tries to have all partitions roughly the same length, but the
#' apportionment is greedy and hence likely not the best possible option.
#'
#' @param x numeric vector
#' @param bw bandwidth of the bounded kernel function
#' @param subdivisions desired number of subdivisions
#' @param range the allowable range of `x`
#' @importFrom rlang warn abort
#' @return a list of intervals
#' @keywords internal
partition_domain <- function(x, bw, subdivisions, range) {
  if (missing(range)) {
    range <- c(-Inf, Inf)
  }
  xdim <- dim(x)
  if (!is.null(xdim) && length(xdim) > 1L && xdim[[2]] > 1L) {
    abort("Domain partitioning only available for univariate data.")
  }
  if (length(x) < 1L) return(list())

  # --- Step 1: Union of intervals x[[i]] +/- bw ---
  x <- sort(x)

  # 2. Identify indices where a new disjoint interval begins.
  #    A break occurs if distance to previous point > 2 * bw.
  #    We prepend 1 because the first element always starts a group.
  group_starts_idx <- c(1L, which(diff(x) > 2 * bw) + 1L)

  # 3. Identify indices where an interval ends.
  #    These correspond to the points just before the breaks, plus the final point.
  group_ends_idx   <- c(group_starts_idx[-1] - 1L, length(x))

  # 4. Construct the disjoint intervals
  #    Start of interval i is min(group i) - bw
  #    End of interval i is max(group i) + bw
  merged_starts <- pmax(x[group_starts_idx] - bw, range[[1]])
  merged_ends   <- pmin(x[group_ends_idx]   + bw, range[[2]])

  # --- Step 2: Apportionment (Assign Counts to Intervals) ---

  num_intervals <- length(merged_starts)
  interval_lens <- merged_ends - merged_starts

  # Edge case: If user asks for fewer partitions than there are islands,
  # we must return at least one per island to cover the support.
  if (subdivisions < num_intervals) {
    warn(paste("The requested number of subdivisions is smaller than the number of disjoint intervals.",
               "Using as many subdivisions as intervals."))
    subdivisions <- num_intervals
  }

  # Initialize: Give every interval at least 1 partition
  counts <- rep.int(1, num_intervals)
  partitions_assigned <- sum(counts)

  # Greedy Loop: Assign remaining partitions to the interval
  # that currently has the LARGEST individual segment size.
  while (partitions_assigned < subdivisions) {
    # Current length of a single segment within each interval
    current_seg_lens <- interval_lens / counts

    # Find interval with the largest segments
    idx <- which.max(current_seg_lens)

    # Add a partition there to shrink its segment size
    counts[[idx]] <- counts[[idx]] + 1L
    partitions_assigned[[1]] <- partitions_assigned[[1]] + 1L
  }

  # --- Step 3: Construct the Final Partitions ---
  final_partitions <- vector("list", subdivisions)
  counter <- 1L

  for (i in seq_len(num_intervals)) {
    # Start/End of this consecutive stretch of non-zero KDE
    s <- merged_starts[[i]]
    e <- merged_ends[[i]]
    k <- counts[[i]] # Number of pieces for this stretch

    # Generate simple equal cuts for this stretch
    cuts <- seq(s, e, length.out = k + 1L)

    # Store each segment as a separate partition in the output list
    for (j in seq_len(k)) {
      final_partitions[[counter]] <- c(cuts[[j]], cuts[[j + 1]])
      counter[[1]] <- counter[[1]] + 1
    }
  }

  return(final_partitions)
}

#' Functor to Create Fixed-grid Gauss-Kronrod Quadrature For Hellinger Affinity
#'
#' This creates a function to compute
#' \deqn{\int \sqrt{ \hat f(x) g(x) }\;\mathrm{d}x}
#' or
#' \deqn{\int \exp\{\frac{1}{2} (\log(\hat f(x)) + \log(g(x)) )\}\;\mathrm{d}x}
#' for arbitrary non-negative functions \eqn{g()}.
#' The function is set up such that \eqn{\hat f} is evaluated exactly 21 times in each
#' subdivision, and does not need to be re-evaluated for different \eqn{g()}.
#'
#' The subdivisions are chosen to cover only regions where \eqn{f(x) > 0} using
#' a greedy partitioning algorithm.
#'
#' @return a function to evaluate the integral
#'   \deqn{2 - 2 \int \sqrt{ \hat f(x) g(x) }\;\mathrm{d}x}
#' @keywords internal
hd_gauss_quadrature <- function (x, wgts, bandwidth, n_subdivisions = 256,
                                 range, kernel = epanechnikov_kernel) {
  # Evaluate the KDE at the Gauss-Kronrod evaluation points
  # We don't need to go beyond the range ±bandwidth because the KDE
  # will be 0 outside.
  range_x <- range(x) + c(-bandwidth, bandwidth)
  if (!missing(range)) {
    range_x[[1]] <- max(range_x[[1]], range[[1]])
    range_x[[2]] <- min(range_x[[2]], range[[2]])
  }

  subintervals <- partition_domain(x, subdivisions = n_subdivisions, bw = bandwidth,
                                   range = range_x)

  log_sum_wgts <- log(sum(wgts))
  log_bw <- log(bandwidth)

  gkpts <- lapply(subintervals, \(subint) {
    lwr <- subint[[1]]
    upr <- subint[[2]]
    kronrod <- cbind(.kronrod_basis, numeric(nrow(.kronrod_basis)))

    gauss_scale <- 0.5 * (upr - lwr)

    kronrod[, 1L] <- gauss_scale * (kronrod[, 1L] + 1) + lwr
    kronrod[, 2L] <- gauss_scale * kronrod[, 2L]
    kronrod[, 3L] <- gauss_scale * kronrod[, 3L]

    fhat <- vapply(kronrod[, 1L] , FUN.VALUE = numeric(1), FUN = \(gp) {
      u <- (x - gp) / bandwidth
      log(sum(wgts * kernel(u)))
    })
    kronrod[, 4L] <- fhat - log_bw - log_sum_wgts

    kronrod
  }) |>
    do.call(what = rbind)

  colnames(gkpts) <- c("x", "wgt_kronrod", "wgt_gauss", "log_f_hat")

  function (h, only_value = TRUE, log = TRUE) {
    h <- h(gkpts[, "x"])
    if (is.null(dim(h))) {
      h <- matrix(h, ncol = 1)
    }
    kronrod <- if(isTRUE(log)) {
      apply(h, 2, \(hi) {
        crossprod(gkpts[, "wgt_kronrod"], exp(0.5 * (gkpts[, "log_f_hat"] + hi)))[[1, 1]]
      })
    } else {
      apply(h, 2, \(hi) {
        crossprod(gkpts[, "wgt_kronrod"], exp(0.5 * gkpts[, "log_f_hat"]) * hi)[[1, 1]]
      })
    }
    if (isTRUE(only_value)) {
      return(2 - 2 * kronrod)
    }

    gauss <- if(isTRUE(log)) {
      apply(h, 2, \(hi) {
        crossprod(gkpts[, "wgt_gauss"], exp(0.5 * (gkpts[, "log_f_hat"] + hi)))[[1, 1]]
      })
    } else {
      apply(h, 2, \(hi) {
        crossprod(gkpts[, "wgt_gauss"], exp(0.5 * gkpts[, "log_f_hat"]) * hi)[[1, 1]]
      })
    }
    list(value = 2 - 2 * kronrod,
         abs.error = 2 * abs(kronrod - gauss),
         subdivisions = length(subintervals))
  }
}
