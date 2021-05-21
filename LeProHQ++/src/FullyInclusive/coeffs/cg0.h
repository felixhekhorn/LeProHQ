#ifndef FullyInclusive_coeffs_cg0_H
#define FullyInclusive_coeffs_cg0_H

#include "config.h"

namespace FullyInclusive {

namespace coeffs {

#define cg0_proj(proj) cdbl cg0_##proj(cdbl m2, cdbl q2, cdbl s);
interateAllProj(cg0_proj)

} // namespace coeffs

} // namespace FullyInclusive

#endif // FullyInclusive_coeffs_cg0_H