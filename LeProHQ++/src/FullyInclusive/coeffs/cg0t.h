#ifndef FullyInclusive_coeffs_cg0t_H
#define FullyInclusive_coeffs_cg0t_H

#include "config.h"

namespace FullyInclusive {

namespace coeffs {
    
#define cg0t_proj(proj) cdbl cg0t_##proj(cdbl m2, cdbl q2, cdbl s);
interateAllProj(cg0t_proj)

} // namespace coeffs

} // namespace FullyInclusive

#endif // FullyInclusive_coeffs_cg0t_H