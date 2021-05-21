#ifndef FullyInclusive_coeffs_cqBarF1_H
#define FullyInclusive_coeffs_cqBarF1_H

#include "config.h"

namespace FullyInclusive {

namespace coeffs {

cdbl cqBarF1_h1(cdbl m2, cdbl q2, cdbl s);

#define cqBarF1_proj(proj) cdbl cqBarF1_##proj(cdbl m2, cdbl q2, cdbl s);
interateAllProj(cqBarF1_proj)

} // namespace coeffs

} // namespace FullyInclusive

#endif // FullyInclusive_coeffs_cqBarF1_H