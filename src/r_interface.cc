#include <R_ext/Rdynload.h>

#include "kde.hpp"

//! R initialzing function (must be in the global namespace).
extern "C" void R_init_mdesurvey(DllInfo *dll) noexcept;

using namespace mdesurvey::r_interface;

namespace {
//! Exported methods
const R_CallMethodDef kExportedCallMethods[] = {
  {"C_kde_univariate", (DL_FUNC) &KDEUnivariate, 5},
  {NULL, NULL, 0}
};
}  // namespace

extern "C" void R_init_mdesurvey(DllInfo *dll) noexcept {
  R_registerRoutines(dll, NULL, kExportedCallMethods, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
