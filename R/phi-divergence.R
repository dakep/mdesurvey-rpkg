.phi_divergence_register <- new.env(parent = emptyenv())

#' Get or Define a Phi Divergence
#'
#' A phi divergence is defined through its generator \eqn{\phi} and the model density
#' exponent \eqn{\lambda} as
#' \deqn{D_{\phi,\lambda}(p \| q) := \int q^\lambda(x) \phi(p(x) / q(x)) \mathrm{d} x.}
#'
#' @param name name of the phi divergence.
#'   If missing, list all currently known divergences.
#' @param ... passed on to initialize a new [PhiDivergence].
#'   Note that if the exponent for (model density) \eqn{q} is different from 1
#'   care must be taken that \eqn{D_{\phi,\lambda}} is a valid divergence!
#' @returns an object of type [PhiDivergence].
#' @importFrom rlang abort dots_n
#' @export
phi_divergence <- function (name, ...) {
  if (missing(name)) {
    return(objects(name = .phi_divergence_register))
  } else if (is.character(name)) {
    if (exists(name, where = .phi_divergence_register) ||
        exists(str_to_lower(name), where = .phi_divergence_register)) {
      if (dots_n(...) > 0L) {
        abort(paste(sprintf("Divergence with name %s already exists.", name),
                    "In this case `...` must be empty!"))
      }
      if (exists(name, where = .phi_divergence_register)) {
        return(.phi_divergence_register[[name]])
      }

      return(.phi_divergence_register[[str_to_lower(name)]])
    } else {
      return(PhiDivergence$new(name = name, ...))
    }
  }
  abort("`name` must be a character")
}

#' @description
#' The convenience function `alpha_divergence()` defines a specific
#' \eqn{\alpha} divergence with \eqn{\lambda = 1} and  generator
#' \deqn{f(x) = \frac{x^\alpha - \alpha x - (1 - \alpha)}{\alpha (\alpha - 1)}.}
#'
#' @param alpha value for the \eqn{\alpha} parameter.
#' @param specialized whether specialized (named) Phi divergences should be returned for
#'   certain values of `alpha`.
#' @export
#' @rdname phi_divergence
alpha_divergence <- function (alpha, specialized = TRUE) {
  if (length(alpha) != 1L || !is.numeric(alpha) || !is.finite(alpha)) {
    abort("`alpha` must be a finite scalar")
  }
  if (isTRUE(specialized)) {
    if (abs(alpha - 1) < .Machine$double.eps) {
      return(phi_divergence("KL"))
    } else if (abs(alpha) < .Machine$double.eps) {
      return(phi_divergence("ReverseKL"))
    } else if (abs(alpha - 0.5) < .Machine$double.eps) {
      return(phi_divergence("Hellinger"))
    }
  }

  # Drop the constant 1 / alpha from w1 and w2 for numerical stability
  PhiDivergence$new(name = sprintf("alpha divergence (alpha=%f)", alpha),
                    phi  = \(x) (x^alpha - alpha * x - (1 - alpha)) / (alpha * (alpha - 1)),
                    w1   = \(x) -x^alpha / alpha,
                    w2   = \(x) x^alpha * (alpha - 1) / alpha,
                    phi_deriv_inf = if (alpha < 1) {
                      1 / (1 - alpha)
                    } else if (abs(alpha - 1) < .Machine$double.eps) 1 else Inf,
                    type = "alpha")
}

#' @description
#' The convenience function `power_divergence()` defines a specific
#' power divergence with \eqn{\lambda=1} and generator
#' \deqn{f(x) = \frac{1}{\alpha (\alpha + 1)} (x^{\alpha + 1} - 1).}
#'
#' @param alpha value for the \eqn{\alpha} power parameter.
#' @param specialized whether specialized (named) Phi divergences should be returned for
#'   certain values of `alpha`.
#' @export
#' @rdname phi_divergence
power_divergence <- function (alpha, specialized = TRUE) {
  if (length(alpha) != 1L || !is.numeric(alpha) || !is.finite(alpha)) {
    abort("`alpha` must be a finite scalar")
  }
  if (isTRUE(specialized)) {
    if (abs(alpha + 1) < .Machine$double.eps) {
      return(phi_divergence("ReverseKL"))
    } else if (abs(alpha) < .Machine$double.eps) {
      return(phi_divergence("KL"))
    } else if (abs(alpha + 0.5) < .Machine$double.eps) {
      return(phi_divergence("Hellinger"))
    }
  }

  # Drop the constant 1 / (alpha * (alpha + 1)) from w1 and w2
  # for numerical stability
  PhiDivergence$new(name = sprintf("power divergence (lambda=%f)", alpha),
                    phi  = \(x) (x^(alpha + 1) - 1)     / (alpha * (alpha + 1)),
                    w1   = \(x) -alpha * x^(alpha + 1)  / (alpha * (alpha + 1)),
                    w2   = \(x) alpha^2 * x^(alpha + 1) / (alpha * (alpha + 1)),
                    phi_deriv_inf = if (alpha < 0) 0 else Inf,
                    type = "power")
}

#' Define a Generalized Divergence
#'
#' The family of generalized divergences is defined as
#' \deqn{
#' D_{(\alpha,\lambda)}(\hat f \| g_\theta) =
#'     \frac{1}{1 + \lambda(1 + \alpha)} \int g_\theta^{1 + \alpha}
#'   - \frac{1+\alpha}{( + \lambda(1 + \alpha)) (\alpha - \lambda(1 - \alpha))} \int g_\theta^{\alpha - \lambda(1 - \alpha)} \hat f^{1 + \lambda(1 + \alpha)}
#'   + \frac{1}{\alpha - \lambda(1 - \alpha)} \int \hat f^{1 + \alpha}
#' }
#'
#' @param alpha,lambda the parameters for the generalized divergence.
#' @returns an object of type [GeneralizedDivergence] or [PhiDivergence], depending
#'   on the combination of parameters `alpha` and `lambda`.
#' @importFrom rlang abort
#' @export
generalized_divergence <- function (alpha, lambda) {
  if (length(alpha) != 1L || !all(is.numeric(alpha)) || !all(is.finite(alpha))) {
    abort("`alpha` must be a scalar, finite number")
  }
  if (length(lambda) != 1L || !all(is.numeric(lambda)) || !all(is.finite(lambda))) {
    abort("`lambda` must be a scalar, finite number")
  }

  if (alpha < -.Machine$double.eps) {
    abort("`alpha` must be >= 0")
  } else if (abs(alpha) < .Machine$double.eps) {
    return(power_divergence(lambda))
  } else if (abs(alpha - 1) < .Machine$double.eps) {
    lambda <- 1
  }

  GeneralizedDivergence$new(alpha = alpha, lambda = lambda)
}

#' @title GeneralizedDivergence Class
#' @description A class to encapsulate generalized divergences defined as
#' \deqn{
#' D_{\alpha,\lambda}(\hat f \| g_\theta) =
#'     \frac{1}{1 + \lambda(1 + \alpha)} \int \hat f^{1 + \alpha}
#'   - \frac{1+\alpha}{( + \lambda(1 + \alpha)) (\alpha - \lambda(1 - \alpha))} \int \hat f^{\alpha - \lambda(1 - \alpha)} g_\theta^{1 + \lambda(1 + \alpha)}
#'   + \frac{1}{\alpha - \lambda(1 - \alpha)} \int g_\theta^{1 + \alpha}
#' }
#' @importFrom R6 R6Class
GeneralizedDivergence <- R6Class(
  classname = "GeneralizedDivergence",

  public = list(
    #' @field alpha parameter of the divergence generator
    alpha = NULL,

    #' @field lambda parameter of the divergence generator
    lambda = NULL,

    #' @field const_A pre-computed `1 + lambda * (1 - alpha)`
    const_A = NULL,

    #' @field const_B pre-computed `alpha - lambda * (1 - alpha)`
    const_B = NULL,

    #' @description
    #' Define a new generalized divergence.
    #'
    #' @param alpha the \eqn{\alpha} parameter of the divergence (between 0 and 1)
    #' @param lambda the \eqn{\lambda} parameter of the divergence
    #' @importFrom rlang abort
    initialize = \(alpha, lambda) {
      self$alpha <- alpha[[1]]
      self$lambda <- lambda[[1]]

      if (self$alpha < .Machine$double.eps) {
        abort("`alpha` must be greater an 0.")
      }

      self$const_A <- 1 + self$lambda * (1 - self$alpha)
      self$const_B <- self$alpha - self$lambda * (1 - self$alpha)

      if (abs(self$const_A) < .Machine$double.eps) {
        abort("`1 + lambda * (1 - alpha)` must be != 0!")
      }
      if (abs(self$const_B) < .Machine$double.eps) {
        abort("`alpha - lambda * (1 - alpha)` must be != 0!")
      }
    }
  )
)

#' @title PhiDivergence Class
#' @description A class to encapsulate generalized phi divergences.
#' @importFrom R6 R6Class
PhiDivergence <- R6Class(
  classname = "PhiDivergence",

  public = list(
    #' @field name full name of the divergence measure
    name = NULL,

    #' @field type type of this phi divergence
    type = NULL,

    #' @field lambda exponent for the model density \eqn{q^\lambda}.
    lambda = NULL,

    #' @field phi divergence generator \eqn{\phi()}
    phi = NULL,

    #' @field w1 function to evaluate \eqn{\phi(x) - x \phi'(x)}
    w1 = NULL,

    #' @field w2 function to evaluate \eqn{x^2 \phi''(x)}
    w2 = NULL,

    #' @field phi_2nd_deriv_at_1 value of \eqn{\phi''(1)}
    phi_2nd_deriv_at_1 = NULL,

    #' @field phi_deriv_inf value of \eqn{\phi'(\infty)}
    phi_deriv_inf = NULL,

    #' @field lambda_neq_1 whether \eqn{\lambda \neq 1}
    lambda_neq_1 = NULL,

    #' @description
    #' Define a new phi divergence.
    #'
    #' @param name the name of phi divergence.
    #' @param phi the generator function.
    #'   The generator must be convex function \eqn{\phi\colon [0, \infty) \to (-\infty, \infty]}
    #'   and satisfy \eqn{\phi(1) = 1}, \eqn{\lim_{t \to 0^+} \phi(t) = \phi(0)}.
    #' @param w1 function to evaluate \eqn{\phi(x) - x \phi'(x)}
    #' @param w2 function to evaluate \eqn{\phi(x) - x \phi'(x) + x^2 \phi''(x)}
    #' @param phi_deriv_inf the value of \eqn{\phi'(\infty)} (can be `Inf`).
    #' @param lambda the model density exponent
    #' @param type a type description for this divergence
    #' @importFrom stringr str_to_lower
    #' @importFrom rlang abort
    initialize = \(name, phi, w1, w2, phi_deriv_inf, lambda = 1, type = 'user-defined') {
      self$name <- name

      self$phi <- match.fun(phi)
      self$w1 <- match.fun(w1)
      self$w2 <- match.fun(w2)
      self$lambda <- as.numeric(lambda[[1]])
      self$phi_deriv_inf <- as.numeric(phi_deriv_inf[[1]])
      lambda_neq_1 <- isTRUE(abs(self$lambda - 1) > .Machine$double.eps)

      if (!missing(type)) {
        self$type <- type[[1]]
      }

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

      self$phi_2nd_deriv_at_1 <- self$w2(1) - self$w1(1)

      if (!is.finite(self$phi_2nd_deriv_at_1)) {
        abort("Cannot evaluate phi''(1)")
      }
    },

    #' @description
    #' Check if this Phi divergence is of a specific type.
    #'
    #' @param type the id to check against
    is = \(type) {
      identical(type[[1]], self$type)
    }
  )
)

.phi_divergence_register$hellinger <-
  .phi_divergence_register$hd <-
  .phi_divergence_register$hellingerdistance <- PhiDivergence$new(
    name = "Hellinger",
    type = 'hd',
    phi  = \(x) 1 - sqrt(x),
    w1   = \(x) -0.5 * sqrt(x),  # + 1 for w1 AND w2 does not change the integral (Bartlett)
    w2   = \(x) -0.25 * sqrt(x),
    phi_deriv_inf = 0)

.phi_divergence_register$negativeexponential <-
  .phi_divergence_register$ned <- PhiDivergence$new(
  name = "NegativeExponential",
  type = 'ned',
  phi  = \(x) expm1(1 - x),
  w1   = \(x) (x + 1)       * exp(1 - x), # - 1 for w1 AND w2 does not change the integral (Bartlett)
  w2   = \(x) (1 + x + x^2) * exp(1 - x),
  phi_deriv_inf = 0)

.phi_divergence_register$kl <-
  .phi_divergence_register$kullbackleibler <- PhiDivergence$new(
    name = "KL",
    type = 'kl',
    phi  = \(x) {
      zeros <- which(x < .Machine$double.eps)
      x[] <- x * log(x)
      x[zeros] <- 0
      x
    },
    w1   = \(x) -x,
    w2   = \(x) 0,
    phi_deriv_inf = Inf)

.phi_divergence_register$reversekl <-
  .phi_divergence_register$revkl <-
  .phi_divergence_register$reversekullbackleibler <- PhiDivergence$new(
    name = "Reverse-KL",
    type = 'rkl',
    phi  = \(x) -log(x),
    w1   = \(x) 1 - log(x),
    w2   = \(x) 2 - log(x),
    phi_deriv_inf = 0)

lockEnvironment(.phi_divergence_register, bindings = TRUE)
