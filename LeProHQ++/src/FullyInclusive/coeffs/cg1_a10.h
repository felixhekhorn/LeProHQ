#ifndef FullyInclusive_coeffs_cg1_a10_H
#define FullyInclusive_coeffs_cg1_a10_H

#include "config.h"
#include "FullyInclusive/PartonicVars.hpp"

namespace FullyInclusive {

namespace coeffs {

cdbl cg1_a10_g1(cdbl m2, cdbl q2);
cdbl cg1_a10_g2(cdbl m2, cdbl q2);

// full/combined
#define cg1_a10_proj(proj) \
cdbl cg1_a10_##proj##_OK(cdbl m2, cdbl q2);\
cdbl cg1_a10_##proj##_QED(cdbl m2, cdbl q2);
interateAllProj(cg1_a10_proj)

} // namespace coeffs

} // namespace FullyInclusive

#endif // FullyInclusive_coeffs_cg1_a10_H