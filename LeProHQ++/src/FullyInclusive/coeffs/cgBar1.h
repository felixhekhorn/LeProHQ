#ifndef FullyInclusive_coeffs_cgBar1_H
#define FullyInclusive_coeffs_cgBar1_H

#include "config.h"
#include "FullyInclusive/PartonicVars.hpp"

namespace FullyInclusive {

namespace coeffs {

cdbl cgBar1_h1(cdbl m2, cdbl q2, cdbl s);
cdbl cgBar1_h2(cdbl m2, cdbl q2, cdbl s);
cdbl cgBar1_h3(cdbl m2, cdbl q2, cdbl s);

#define cgBar1_proj(proj) cdbl cgBar1_##proj(cdbl m2, cdbl q2, cdbl s);
interateAllProj(cgBar1_proj)

} // namespace coeffs

} // namespace FullyInclusive

#endif // FullyInclusive_coeffs_cgBar1_H