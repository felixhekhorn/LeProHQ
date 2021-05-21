#ifndef FullyInclusiveSoftResummed_MellinVars_HPP
#define FullyInclusiveSoftResummed_MellinVars_HPP

#include <complex>

#include "config.h"

namespace FullyInclusiveSoftResummed {

/** @brief complex number */
typedef complex<dbl> dcmplx;

/** @brief constant complex number */
typedef const dcmplx cdcmplx;

/** matrix element */
typedef cdcmplx (*fPtr2dblN)(cdbl m2, cdbl q2, cdcmplx N);

} // namespace FullyInclusiveSoftResummed

#endif // FullyInclusiveSoftResummed_MellinVars_HPP