#include "r_interface_utilities.hpp"

using arma::vec;

namespace mdesurvey {

//! Get an unsafe view to the given R vector without copying any data.
//!
//! @param numeric_vector a numeric R vector
std::unique_ptr<const arma::vec> MakeNumericVectorView(SEXP numeric_vector) noexcept {
  if (TYPEOF(numeric_vector) != REALSXP) {
    return std::make_unique<const vec>();
  }
  return std::make_unique<const vec>(REAL(numeric_vector), Rf_length(numeric_vector), false, true);
}

} // namespace mdesurvey

