#' Partition the Domain Into Subdivisions
#'
#' @details
#' ## Partitioning
#' The integration considers only the "relevant" domain spanned by `x`.
#' The span of `x` is partitioned into a fixed number of partitions such that the
#' KDE \eqn{\hat{f}(x) > 0} within each of those partitions and \eqn{\hat{f}(x) = 0} outside
#' those partitions.
#' The algorithm for `partition_domain()` tries to have all partitions roughly the same length, but the
#' apportionment is greedy and hence possibly not the globally best solution.
#'
#' @param x numeric vector
#' @param bw bandwidth of the bounded kernel function
#' @param subdivisions number of partitions
#' @param range the allowable range of `x`. If missing, the assumed allowable range is
#'   \eqn{[\min_i x_i - bw, \max_i x_i + bw]}.
#' @importFrom rlang warn abort
#' @keywords internal
#' @rdname hd_gauss_quadrature
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
