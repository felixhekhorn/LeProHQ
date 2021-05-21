#include "FullyInclusiveSoftResummed/ResExp.h"

using FullyInclusiveSoftResummed::cdcmplx;

cdcmplx FullyInclusiveSoftResummed::gLL(dbl Ag1, dbl beta0, cdcmplx lambda) {
    return Ag1 / M_PI / beta0 * (1. + (1. -2.*lambda)*log(1.-2.*lambda)/(2.*lambda));
}