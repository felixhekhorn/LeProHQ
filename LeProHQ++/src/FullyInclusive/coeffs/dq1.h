#ifndef FullyInclusive_coeffs_dq1_H
#define FullyInclusive_coeffs_dq1_H

#include "config.h"

namespace FullyInclusive {

namespace coeffs {

cdbl dq1_h4(cdbl m2, cdbl q2, cdbl s);

#define dq1_proj(proj) cdbl dq1_##proj(cdbl m2, cdbl q2, cdbl s);
interateAllProj(dq1_proj)

} // namespace coeffs

} // namespace FullyInclusive

#endif // FullyInclusive_coeffs_dq1_H