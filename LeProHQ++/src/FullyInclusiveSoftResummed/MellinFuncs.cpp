#include "FullyInclusiveSoftResummed/MellinFuncs.h"
#include <gsl/gsl_sf.h>

using FullyInclusiveSoftResummed::cdcmplx;

cdcmplx FullyInclusiveSoftResummed::lngamma(cdcmplx N) {
    gsl_sf_result lnr;
    gsl_sf_result arg;
    gsl_sf_lngamma_complex_e(N.real(), N.imag(), &lnr, &arg);
    return dcmplx(lnr.val, arg.val);
}