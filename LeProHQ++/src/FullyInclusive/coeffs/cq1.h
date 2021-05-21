#ifndef FullyInclusive_coeffs_cq1_H
#define FullyInclusive_coeffs_cq1_H

#include "config.h"
#include "FullyInclusive/PartonicVars.hpp"

namespace FullyInclusive {

namespace coeffs {

// true threshold
cdbl cq1t_(cdbl m2, cdbl q2, cdbl s, fPtr3dbl cg0t, cdbl a11, cdbl a10);
#define cq1t_proj(proj) cdbl cq1t_##proj(cdbl m2, cdbl q2, cdbl s);
interateAllProj(cq1t_proj)

// improved threshold
cuint cq1tp_lnxi_num = 21;
cdbl cq1tp_lneta_max = log(1e0);
cdbl cq1tp_(cdbl m2, cdbl q2, cdbl s, fPtr3dbl cq1t, str proj_cc);
#define cq1tp_proj(proj) cdbl cq1tp_##proj(cdbl m2, cdbl q2, cdbl s);
interateAllProj(cq1tp_proj)

// bulk
cuint cq1b_lnxi_num = 21;
cdbl cq1b_lneta_min = log(.9e0);
cdbl cq1b_lneta_max = log(1e5);
cuint cq1b_lneta_num = 51;
cdbl cq1b_(cdbl m2, cdbl q2, cdbl s, str proj_cc);
#define cq1b_proj(proj) cdbl cq1b_##proj(cdbl m2, cdbl q2, cdbl s);
interateAllProj(cq1b_proj)

// full/combined
cdbl cq1_(cdbl m2, cdbl q2, cdbl s, fPtr3dbl cq1tp, fPtr3dbl cq1b);
#define cq1_proj(proj) cdbl cq1_##proj(cdbl m2, cdbl q2, cdbl s);
interateAllProj(cq1_proj)

} // namespace coeffs

} // namespace FullyInclusive

#endif // FullyInclusive_coeffs_cq1_H