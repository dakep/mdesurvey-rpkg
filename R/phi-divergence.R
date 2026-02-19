.phi_divergence_register <- new.env(parent = emptyenv())

#' Get or Define a Phi Divergence
#'
#' A phi divergence is defined through its generator \eqn{\phi} as
#' \deqn{D_\phi(p \| q) := \int q(x) \phi(p(x) / q(x)) \mathrm{d} x.}
#'
#' @param name name of the phi divergence.
#'   If missing, list all currently known divergences.
#' @param phi,w1,w2 the generator of the phi divergence and function of its derivatives.
#'   Passed on to initialize a new [PhiDivergence].
#' @returns an object of type [PhiDivergence].
#' @importFrom rlang abort
#' @export
phi_divergence <- function (name, phi, w1, w2) {
  if (missing(name)) {
    return(objects(name = .phi_divergence_register))
  } else if (is.character(name)) {
    if (exists(name, where = .phi_divergence_register) ||
        exists(str_to_lower(name), where = .phi_divergence_register)) {
      if (!missing(phi) || !missing(w1) || !missing(w2)) {
        abort(sprintf("Divergence with name %s already exists. In this case `phi` must be empty!"))
      }

      if (exists(name, where = .phi_divergence_register)) {
        return(.phi_divergence_register[[name]])
      }

      return(.phi_divergence_register[[str_to_lower(name)]])
    } else {
      return(PhiDivergence$new(name = name, phi = phi, w1 = w1, w2 = w2))
    }
  }
  abort("`name` must be a character")
}

#' @description
#' The convenience function `alpha_divergence()` defines a specific
#' \eqn{\alpha} divergence with generator
#' \deqn{f(x) = \frac{x^\alpha - \alpha x - (1 - \alpha)}{\alpha (\alpha - 1)}.}
#'
#' @param alpha value for the \eqn{\alpha} parameter.
#' @export
#' @rdname phi_divergence
alpha_divergence <- function (alpha) {
  if (abs(alpha - 1) < .Machine$double.eps) {
    return(phi_divergence("KL"))
  } else if (abs(alpha) < .Machine$double.eps) {
    return(phi_divergence("Reverse-KL"))
  } else if (abs(alpha - 0.5) < .Machine$double.eps) {
    return(phi_divergence("Hellinger"))
  }

  PhiDivergence$new(name = sprintf("alpha divergence (alpha=%f)", alpha),
                    phi  = \(x) (x^alpha - alpha * x - (1 - alpha)) / (alpha * (alpha - 1)),
                    w1   = \(x) (1 - x^alpha) / alpha,
                    w2   = \(x) x^alpha)
}

#' @description
#' The convenience function `power_divergence()` defines a specific
#' power divergence with generator
#' \deqn{f(x) = \frac{1}{\lambda (\lambda + 1)} (x^{\lambda + 1} - 1).}
#'
#' @param lambda value for the \eqn{\lambda} power parameter.
#' @export
#' @rdname phi_divergence
power_divergence <- function (lambda) {
  if (abs(lambda + 1) < .Machine$double.eps) {
    return(phi_divergence("ReverseKL"))
  } else if (abs(lambda) < .Machine$double.eps) {
    return(phi_divergence("KL"))
  } else if (abs(lambda + 0.5) < .Machine$double.eps) {
    return(phi_divergence("Hellinger"))
  }

  PhiDivergence$new(name = sprintf("power divergence (lambda=%f)", lambda),
                    phi  = \(x) (x^(lambda + 1) - 1) / (lambda * (lambda + 1)),
                    w1   = \(x) -(1 + lambda * x^(lambda + 1)) / (lambda * (lambda + 1)),
                    w2   = \(x) x^(lambda + 1))
}

#' @title PhiDivergence Class
#' @description A class to encapsulate different phi divergences.
#' @importFrom R6 R6Class
PhiDivergence <- R6Class(
  classname = "PhiDivergence",

  public = list(
    #' @field name full name of the divergence measure
    name = NULL,

    #' @field phi divergence generator \eqn{\phi()}
    phi = NULL,

    #' @field w1 function to evaluate \eqn{\phi(x) - x \phi'(x)}
    w1 = NULL,

    #' @field w2 function to evaluate \eqn{x^2 \phi''(x)}
    w2 = NULL,

    #' @description
    #' Define a new phi divergence.
    #'
    #' @param name the name of phi divergence.
    #' @param phi the generator function.
    #'   The generator must be convex function \eqn{\phi\colon [0, \infty) \to (-\infty, \infty]}
    #'   and satisfy \eqn{\phi(1) = 1}, \eqn{\lim_{t \to 0^+} \phi(t) = \phi(0)}.
    #' @param w1 function to evaluate \eqn{\phi(x) - x \phi'(x)}
    #' @param w2 function to evaluate \eqn{x^2 \phi''(x)}
    #' @importFrom stringr str_to_lower
    initialize = \(name, phi, w1, w2) {
      self$name <- name

      self$phi <- match.fun(phi)
      self$w1 <- match.fun(w1)
      self$w2 <- match.fun(w2)

      tryCatch({
        f0 <- self$phi(0)
        if (is.na(f0) || is.nan(f0)) {
          abort("phi(0) is not a number")
        }
      }, error = \(cnd) {
        abort(parent = cnd,
              message = 'Cannot evaluate the phi function at 0.')
      })
      if (!isTRUE(abs(self$phi(1)) < .Machine$double.eps)) {
        abort("phi(1) must be exactly 0")
      }
    }
  )
)

.phi_divergence_register$hellinger <-
  .phi_divergence_register$hd <-
  .phi_divergence_register$hellingerdistance <- PhiDivergence$new(
    name        = "Hellinger",
    phi         = \(x) (sqrt(x) - 1)^2,
    w1          = \(x) (1 - sqrt(x)),
    w2          = \(x) 0.25 * sqrt(x))
    # phi         = \(x) 1 - sqrt(x),
    # w1          = \(x) 1 - 0.5 * sqrt(x),
    # w2          = \(x) 0.25 * sqrt(x))
    # f_deriv     = \(x) -0.5 / sqrt(x),
    # f_deriv_2nd = \(x) 0.25 * x^-1.5,

.phi_divergence_register$negativeexponential <-
  .phi_divergence_register$ned <- PhiDivergence$new(
  name        = "NegativeExponential",
  phi         = \(x) expm1(1 - x),
  w1          = \(x) expm1(1 - x) + x * exp(1 - x),
  w2          = \(x) x^2 * exp(1 - x))
  # f_deriv     = \(x) -exp(1 - x),
  # f_deriv_2nd = \(x) exp(1 - x),

.phi_divergence_register$kl <-
  .phi_divergence_register$kullbackleibler <- PhiDivergence$new(
    name = "KL",
    phi = \(x) {
      zeros <- which(x < .Machine$double.eps)
      x[] <- x * log(x)
      x[zeros] <- 0
      x
    },
    w1          = \(x) -x,
    w2          = \(x) x)
    # f_deriv     = \(x) 1 + log(x),
    # f_deriv_2nd = \(x) 1 / x,

.phi_divergence_register$reversekl <-
  .phi_divergence_register$revkl <-
  .phi_divergence_register$reversekullbackleibler <- PhiDivergence$new(
    name        = "Reverse-KL",
    phi         = \(x) -log(x),
    w1          = \(x) 1 - log(x),
    w2          = \(x) 1)
    # f_deriv     = \(x) -1/x,
    # f_deriv_2nd = \(x) 1/x^2,

lockEnvironment(.phi_divergence_register, bindings = TRUE)
