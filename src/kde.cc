#include <RcppArmadillo.h>
#include "r_interface_utilities.hpp"

namespace {
template <typename KernelFunc>
Rcpp::NumericVector EvaluateUnivariateKDE(const arma::vec& x, const arma::vec& wgts,
                                          const arma::vec& eval, const double bw, KernelFunc kernel) {
  // 1. Sort x and the accompanying weights
  const arma::uvec indices = arma::sort_index(x);
  const arma::vec x_sorted = x(indices);
  const arma::vec wgts_sorted = wgts(indices);

  const arma::uword n = x_sorted.n_elem;
  const arma::uword m = eval.n_elem;
  Rcpp::NumericVector out(Rcpp::no_init(m));

  // 2. Initialize sliding window pointers
  arma::uword left = 0;
  arma::uword right = 0;

  const double sum_wgts_bw = arma::accu(wgts) * bw;

  // 3. Iterate through sorted evaluation points
  for (arma::uword i = 0; i < m; ++i) {
    const double val = eval[i];
    const double lower_bound = val - bw;
    const double upper_bound = val + bw;

    // Advance 'right' pointer to the first element > upper_bound
    // This marks the exclusive end of our range
    while (right < n && x_sorted[right] <= upper_bound) {
      right++;
    }

    // Advance 'left' pointer to the first element >= lower_bound
    // This marks the inclusive start of our range
    while (left < right && x_sorted[left] < lower_bound) {
      left++;
    }

    // 4. Sum kernel contributions in the active window [left, right)
    double sum = 0.0;

    // Note: Using raw pointers here is slightly faster than x_sorted[j]
    // inside the hot loop, but indices are safer and usually optimized well.
    for (int j = left; j < right; ++j) {
      sum += wgts_sorted[j] * kernel(x_sorted[j] - val);
    }

    out[i] = sum / sum_wgts_bw;
  }
  return out;
}


} // namespace

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
//! @return a numeric vector of KDE evaluations
SEXP KDEUnivariate(SEXP r_x, SEXP r_wgts, SEXP r_eval, SEXP r_bw, SEXP r_kernel) noexcept {
  BEGIN_RCPP
  auto         x      = MakeNumericVectorView(r_x);
  auto         wgts   = MakeNumericVectorView(r_wgts);
  auto         eval   = MakeNumericVectorView(r_eval);
  const double bw     = Rcpp::as<double>(r_bw);
  const char   kernel = Rcpp::as<char>(r_kernel);

  Rcpp::NumericVector kde;

  if (kernel == 'b') {
    // Biweight kernel
    auto biweight_kernel = [bw](double d) -> double {
      d /= bw;
      if (d > 1.0) return 0.0;
      d = 1.0 - d * d;
      return 0.9375 * d * d; // note: 15/16 = 0.9375
    };
    kde = EvaluateUnivariateKDE(*x, *wgts, *eval, bw, biweight_kernel);
  } else if (kernel == 'r') {
    // Rectangular kernel
    auto rect_kernel = [bw](double d) -> double {
      d = std::abs(d / bw);
      if (d > 1.0) return 0.0;
      return 0.5;
    };
    kde = EvaluateUnivariateKDE(*x, *wgts, *eval, bw, rect_kernel);
  } else if (kernel == 't') {
    // Triangular kernel
    auto triangular_kernel = [bw](double d) -> double {
      d = std::abs(d / bw);
      if (d > 1.0) return 0.0;
      return 1.0 - d;
    };
    kde = EvaluateUnivariateKDE(*x, *wgts, *eval, bw, triangular_kernel);
  } else {
    // Epanechnikov kernel
    auto epanechnikov_kernel = [bw](double d) -> double {
      d /= bw;
      if (std::abs(d) > 1.0) return 0.0;
      return 0.75 * (1.0 - d * d);
    };
    kde = EvaluateUnivariateKDE(*x, *wgts, *eval, bw, epanechnikov_kernel);
  }

  return Rcpp::wrap(kde);
  END_RCPP
}

} // namespace r_interface
} // namespace mdesurvey
