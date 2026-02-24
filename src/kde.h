#ifndef MDESURVEY_KDE_HPP_
#define MDESURVEY_KDE_HPP_

#include <RcppArmadillo.h>

namespace mdesurvey {
namespace r_interface {

//! Evaluate the Univariate KDE at Multiple Points
//
//! @param x      a numeric vector of observations
//! @param wgts   a numeric vector of observation weights
//! @param eval   an **ordered** numeric vector of evaluation points
//! @param bw     the bandwidth for the kernel (numeric scalar)
//! @param kernel a character indicating the kernel function. Supported choices are
//!  'e' ... Epanechnikov kernel
//!  'b' ... Biweight kernel
//!  'r' ... Rectangular kernel
//!  't' ... Triangular kernel
//! @return a numeric vector of KDE evaluations
SEXP KDEUnivariate(SEXP x, SEXP wgts, SEXP eval, SEXP bw, SEXP kernel) noexcept;

} // namespace r_interface

} // namespace mdesurvey

#endif // MDESURVEY_KDE_HPP_
