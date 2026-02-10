#ifndef MDESURVEY_R_INTERFACE_UTILITIES_HPP_
#define MDESURVEY_R_INTERFACE_UTILITIES_HPP_

#include <RcppArmadillo.h>
#include <memory>

namespace mdesurvey {

//! Get an unsafe view to the given R vector without copying any data.
//!
//! @param numeric_vector a numeric R vector
std::unique_ptr<const arma::vec> MakeNumericVectorView(SEXP numeric_vector) noexcept;

} // namespace mdesurvey

#endif // MDESURVEY_R_INTERFACE_UTILITIES_HPP_
