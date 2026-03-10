.link_register <- new.env(parent = emptyenv())

#' Link Functions
#'
#' The `link()` function can be used to define/access link functions linking the location
#' and nuisance (variance) components to the predictors.
#'
#' @param mean either a character naming one of the built-in link functions (see Details)
#'   or an instance of the [Link] R6 class.
#' @param ... named link functions (either character or an instance of the [Link] class),
#'   one for each nuisance parameter..
#' @export
#' @importFrom rlang abort dots_list
#' @importFrom R6 is.R6
#' @importFrom methods is
link <- function (mean, ...) {
  if (is.character(mean)) {
    mean <- .link_register[[mean]]
  }

  if (!is.R6(mean) || !is(mean, "Link")) {
    abort("`mean` must be an instance of the Link R6 class or the name of a built-in link function")
  }

  nuisance <- dots_list(..., .ignore_empty = 'none', .homonyms = 'error') |>
    lapply(\(ns) {
      if (is.character(ns)) {
        ns <- .link_register[[ns]]
      }

      if (!is.R6(ns) || !is(ns, "Link")) {
        abort(paste("All elements of `nuisance` must be an instance of the Link R6 class,",
                    "or the name of a built-in link function."))
      }
      ns
    })

  if (length(nuisance) > 0L && (is.null(names(nuisance)) || !all(nzchar(names(nuisance))))) {
    abort("All links for the nuisance parameter(s) must be named")
  }

  list(mean = mean, nuisance = nuisance)
}


#' @importFrom rlang abort
#' @importFrom R6 is.R6
#' @importFrom methods is
.check_link <- function (link, nuisance_names) {
  if (!is.R6(link$mean) || !is(link$mean, "Link")) {
    abort("Cannot determine mean link.")
  }
  for (ns in nuisance_names) {
    if (!is.R6(link$nuisance[[ns]]) || !is(link$nuisance[[ns]], "Link")) {
      abort(sprintf("Cannot determin link for the %s nuisance parameter.", ns))
    }
  }
  link$nuisance <- link$nuisance[nuisance_names]
  link
}

#' @title Link Class
#' @description The `Link` class encapsulates the link function and it's inverse & derivative
#'   information.
#' @importFrom R6 R6Class
#' @export
#' @rdname link
Link <- R6Class(
  classname = "Link",
  public = list(
    #' @field name name/description of the link function
    name = NULL,

    #' @field fun link function mapping the mean (location) parameter of a model family
    #'   to the predictors
    fun = NULL,

    #' @field inv inverse link function mapping the mean (location) parameter of a model family
    #'   to the predictors
    inv = NULL,

    #' @field inv_d1 first derivative of the inverse link function
    inv_d1 = NULL,

    #' @field inv_d2 second derivative of the inverse link function
    inv_d2 = NULL,

    #' @description
    #' Define a new link function \eqn{\eta} and its inverse \eqn{\eta^{-1}}
    #' to map the location parameter of a model family to the predictors.
    #'
    #' @param name name of the link function
    #' @param fun the link function \eqn{\eta(\mu) = f(x)}
    #' @param inv the inverse link function \eqn{\mu = \eta^{-1}(f(x))}
    #' @param inv_d1 the first derivative of the inverse link function,
    #'   \eqn{\frac{\partial}{\partial u} \eta^{-1}(u)}
    #' @param inv_d2 the second derivative of the inverse link function,
    #'   \eqn{\frac{\partial^2}{\partial u^2} \eta^{-1}(u)}
    #'
    #' @importFrom rlang abort
    initialize = \ (name, fun, inv, inv_d1, inv_d2) {
      self$name <- name
      self$fun <- match.fun(fun)
      self$inv <- match.fun(inv)
      self$inv_d1 <- match.fun(inv_d1)
      self$inv_d2 <- match.fun(inv_d2)
    },

    #' @description
    #' Print the name and parameters of a model family
    #'
    #' @param ... ignored
    #' @importFrom stringr str_to_title
    print = \(...) {
      cat(paste(str_to_title(self$name), "link"))
    }
  )
)


.link_register$identity <- Link$new(
  name   = "identity",
  fun    = identity,
  inv    = identity,
  inv_d1 = \(...) 1,
  inv_d2 = \(...) 0)

.link_register$log <- Link$new(
  name   = "log",
  fun    = log,
  inv    = exp,
  inv_d1 = exp,
  inv_d2 = exp)

.link_register$inverse <- Link$new(
  name   = "inverse",
  fun    = \(mu) 1 / mu,
  inv    = \(x)  1 / x,
  inv_d1 = \(x) -1 / x^2,
  inv_d2 = \(x)  2 / x^3)


#' @importFrom stats binomial
.bin_logit <- binomial("logit")

.link_register$logit <- Link$new(
  name   = "logit",
  fun    = .bin_logit$linkfun,
  inv    = .bin_logit$linkinv,
  inv_d1 = \(x) {
    x[] <- .bin_logit$linkinv(x)
    x * (1 - x)
  },
  inv_d2 = \(x) {
    x[] <- .bin_logit$linkinv(x)
    x * (1 - x) * (1 - 2 * x)
  })

lockEnvironment(.link_register, bindings = TRUE)
