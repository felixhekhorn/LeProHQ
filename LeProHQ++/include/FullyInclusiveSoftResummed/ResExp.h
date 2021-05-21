#ifndef FullyInclusiveSoftResummed_ResExp_H
#define FullyInclusiveSoftResummed_ResExp_H

#include "../config.h"
#include "MellinVars.hpp"

namespace FullyInclusiveSoftResummed {

/**
 * @brief Leading-logarithmic resummation exponent
 * @param Ag1 leading-order gluon operator
 * @param beta0 leading-order beta function
 * @param lambda resummation parameter
 * @return \f$g^{\text{LL}}(\lambda)\f$
 */
cdcmplx gLL(dbl Ag1, dbl beta0, cdcmplx lambda);

} // namespace FullyInclusiveSoftResummed

#endif // FullyInclusiveSoftResummed_ResExp_H