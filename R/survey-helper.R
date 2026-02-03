#' Extract the Numeric Values from a Survey Design
#' @importFrom stats na.pass weights
#' @importFrom rlang enquo eval_tidy
#' @keywords internal
.extract_survey_values <- function (x, design, na.rm = FALSE) {
  if (!inherits(design, "survey.design")) {
    abort("`design` must be a survey.design object from the survey package")
  }
  x <- eval_tidy(enquo(x), data = design$variables)

  if (inherits(x, "formula")) {
    mf <- model.frame(x, model.frame(design), na.action = na.pass)
    xx <- lapply(attr(terms(x), "variables")[-1],
                 \(term) model.matrix(eval(bquote(~0 + .(term))), mf))
    x <- do.call(cbind, xx)
  } else if (typeof(x) %in% c("expression", "symbol")) {
    x <- eval(x, design$variables)
  }

  x <- as.matrix(x)

  if (isTRUE(na.rm)) {
    nas <- rowSums(is.na(x))
    if (any(nas > 0)) {
      design <- design[nas == 0, ]
    }
    x[nas > 0, ] <-  0
  }

  list(x = x,
       wgts = weights(design))
}
