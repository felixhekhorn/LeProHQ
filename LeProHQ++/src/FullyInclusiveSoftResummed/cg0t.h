#ifndef FullyInclusiveSoftResummed_coeffs_cg0t_H
#define FullyInclusiveSoftResummed_coeffs_cg0t_H

#include "config.h"
#include "FullyInclusiveSoftResummed/MellinVars.hpp"

namespace FullyInclusiveSoftResummed {

namespace coeffs {
    
#define cg0t_proj(proj) cdcmplx cg0t_##proj(cdbl m2, cdbl q2, cdcmplx N);
interateAllProj(cg0t_proj)

} // namespace coeffs

} // namespace FullyInclusiveSoftResummed

#endif // FullyInclusiveSoftResummed_coeffs_cg0t_H