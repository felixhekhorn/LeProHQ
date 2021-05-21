#ifndef FullyInclusive_coeffs_cg1_H
#define FullyInclusive_coeffs_cg1_H

#include "config.h"
#include "FullyInclusive/PartonicVars.hpp"

namespace FullyInclusive {

namespace coeffs {

// true threshold
cdbl cg1t_(cdbl m2, cdbl q2, cdbl s, fPtr3dbl cg0t, cdbl a12, cdbl a11, fPtr2dbl a10_OK, fPtr2dbl a10_QED);
#define cg1t_proj(proj) cdbl cg1t_##proj(cdbl m2, cdbl q2, cdbl s);
interateAllProj(cg1t_proj)


// improved threshold
cuint cg1tp_lnxi_num = 21;
cdbl cg1tp_lneta_max = log(1e-1);
cdbl cg1tp_(cdbl m2, cdbl q2, cdbl s, fPtr3dbl cg1t, str proj_cc);
#define cg1tp_proj(proj) cdbl cg1tp_##proj(cdbl m2, cdbl q2, cdbl s);
interateAllProj(cg1tp_proj)

// bulk
cuint cg1b_lnxi_num = 21;
cdbl cg1b_lneta_min = log(.9e-1);
cdbl cg1b_lneta_max = log(1e5);
cuint cg1b_lneta_num = 51;
cdbl cg1b_(cdbl m2, cdbl q2, cdbl s, str proj_cc);
#define cg1b_proj(proj) cdbl cg1b_##proj(cdbl m2, cdbl q2, cdbl s);
interateAllProj(cg1b_proj)

// full/combined
cdbl cg1_(cdbl m2, cdbl q2, cdbl s, fPtr3dbl cg1tp, fPtr3dbl cg1b);
#define cg1_proj(proj) cdbl cg1_##proj(cdbl m2, cdbl q2, cdbl s);
interateAllProj(cg1_proj)

} // namespace coeffs

} // namespace FullyInclusive

#endif // FullyInclusive_coeffs_cg1_H