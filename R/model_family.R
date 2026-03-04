.model_family_register <- new.env(parent = emptyenv())

#' Access Model Families
#'
#' Get existing model families or define a new model family.
#'
#' @param name name of the model family.
#'   If missing, list all currently known model families.
#' @param ... passed on to initialize a new [ModelFamily].
#' @return an object of type [ModelFamily].
#' @importFrom stringr str_to_lower
#' @importFrom rlang dots_n abort
#' @export
model_family <- function (name, ...) {
  if (missing(name)) {
    return(objects(name = .model_family_register))
  } else if (is.character(name)) {
    if (exists(str_to_lower(name), where = .model_family_register)) {
      if (dots_n(...) > 0L) {
        abort(sprintf("Model family with name %s exists. In this case ... must be empty!",
                      name))
      }

      return(.model_family_register[[str_to_lower(name)]])
    } else {
      return(ModelFamily$new(name = name, ...))
    }
  }
  abort("`name` must be a character")
}


#' @title ModelFamily Class
#' @description A class to encapsulate statistical model families for survey analysis.
#' @importFrom R6 R6Class
ModelFamily <- R6Class(
  classname = "ModelFamily",
  public = list(
    #' @field name name/description of the model family
    name = NULL,

    #' @field parameter_names optional names of the parameters
    parameter_names = NULL,

    #' @field range the range of the model density (numeric matrix); (lower, upper) per dimension.
    range = NULL,

    #' @field dfun Density function `dfun(x, params, log) -> (log) density`
    dfun = NULL,

    #' @field pfun Distribution function `pfun(x, params) -> CDF`
    pfun = NULL,

    #' @field raw_scores function to compute the raw scores for the model. The function
    #'   should return a matrix with one column per parameter.
    #'   `raw_scores(x, params) -> numeric matrix of scores`
    raw_scores = NULL,

    #' @field hessian function to compute the entries of the Hessian for the model,
    #'   \eqn{\frac{\partial}{\partial \theta_{ij} f_\theta}}.
    #'   The function should return a matrix with one column per entry of the covariance matrix.
    #'   It is assumed that first the diagonal elements of the covariance matrix are returned,
    #'   and then the upper-triagonal elements in row-first fashion (as returned by `combn()`), e.g.,
    #'   `(1,1), (2,2), (3,3), (1,2), (1,3), (2,3)`, for a 3-dimensional parameter.
    #'   `hessian(x, params) -> numeric matrix of entries of the model Hessian`
    hessian = NULL,

    #' @field fisher_inf function to compute the Fisher Information matrix.
    #'   `fisher_inf(params) -> Fisher Information (numeric matrix)`
    fisher_inf = NULL,

    # -- Parameter Constraints --
    #' @field trans function to transform the natural parameters to unrestricted parameters (e.g. `log()`).
    #'   `function(params) -> vector of unrestricted parameters`
    trans = NULL,

    #' @field inv_trans function to transform the unrestricted parameters to the natural parameters (e.g. `exp()`).
    #'   `function(u_params) -> vector of natural parameters`
    inv_trans = NULL,

    #' @field initial function to get initial estimates.
    #'   `function(x, design) -> vector of initial estimates`
    initial = NULL,

    #' @field reparameterize_from_mv function to get the family's natural parameters from
    #'   a mean & variance parameterization.
    #'   May not be well defined, in which case the function should return a parameter vector
    #'   of NA's.
    reparameterize_from_mv = NULL,

    #' @field reparameterize_to_mv function to get the mean and variance from the
    #'   family's parameters.
    #'   May not be well defined, in which case the undefined moments should be NA.
    reparameterize_to_mv = NULL,

    #' @description
    #' Define a new model family.
    #'
    #' @param name name of the model family
    #' @param parameter_names names for the parameters.
    #' @param range the range of the model density.
    #'   The density is assumed to be 0 outside this range.
    #' @param dfun Density function of the model distribution.
    #'   `dfun(x, params, log) -> numeric vector of (log) density values`
    #' @param pfun Distribution functions of the model distribution.
    #'   `pfun(x, params) -> numeric vector of CDF values`
    #' @param trans,inv_trans transformation and inverse transformation function to
    #'   map to natural parameters to the unrestricted parameter space (`trans()`)
    #'   and the other direction (`inv_trans()`).
    #' @param raw_scores function to compute the raw scores for the model. The function
    #'   should return a matrix with one column per parameter.
    #'   `raw_scores(x, params) -> numeric matrix of scores`
    #' @param hessian function to compute the entries of the Hessian for the model,
    #'   \eqn{\frac{\partial}{\partial \theta_{ij} f_\theta}}.
    #'   The function should return a matrix with one column per entry of the covariance matrix.
    #'   It is assumed that first the diagonal elements of the covariance matrix are returned,
    #'   and then the upper-triagonal elements in row-first fashion (as returned by `combn()`), e.g.,
    #'   `(1,1), (2,2), (3,3), (1,2), (1,3), (2,3)`, for a 3-dimensional parameter.
    #'   `hessian(x, params) -> numeric matrix of entries of the model Hessian`
    #' @param fisher_inf function to compute the Fisher Information matrix.
    #'   `fisher_inf(params) -> Fisher Information (numeric matrix)`
    #' @param initial function to get initial estimates.
    #'   `function(x, design) -> vector of initial estimates`
    #' @param reparameterize_from_mv function to get the
    #'   family's natural parameters from a mean & variance parameterization.
    #'   May not be well defined, in which case the functions should return a parameter vector
    #'   of NA's.
    #' @param reparameterize_to_mv function to get the mean and variance from the
    #'   family's parameters.
    #'   May not be well defined, in which case the undefined moments should be NA.
    initialize = \ (name, parameter_names, range, dfun, pfun, raw_scores, hessian, fisher_inf,
                    trans, inv_trans, initial, reparameterize_from_mv, reparameterize_to_mv) {
      self$name <- name
      self$parameter_names <- parameter_names
      if (is.null(dim(range))) {
        range <- matrix(range, nrow = 1)
      }
      self$range <- range
      self$dfun <- match.fun(dfun)
      self$pfun <- match.fun(pfun)
      self$raw_scores <- match.fun(raw_scores)
      self$hessian <- match.fun(hessian)
      self$fisher_inf <- match.fun(fisher_inf)
      self$trans <- match.fun(trans)
      self$inv_trans <- match.fun(inv_trans)
      self$initial <- match.fun(initial)

      self$reparameterize_to_mv <- match.fun(reparameterize_to_mv)
      self$reparameterize_from_mv <- if (missing(reparameterize_from_mv)) {
        \(...) {
          params <- rep.int(NA_real_, length(self$reparameterize_from_mv))
          names(params) <- self$parameter_names
          params
        }
      } else {
        match.fun(reparameterize_from_mv)
      }

    },

    #' @description
    #' Print the name and parameters of a model family
    #'
    #' @param ... ignored
    #' @importFrom stringr str_to_title
    print = \(...) {
      cat(sprintf("%s model family with parameters %s",
                  str_to_title(self$name), paste(self$parameter_names, collapse = ", ")))
    }
  )
)

