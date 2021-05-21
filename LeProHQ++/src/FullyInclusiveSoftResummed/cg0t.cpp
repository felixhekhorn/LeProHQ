#include "cg0t.h"
#include "FullyInclusive/PartonicVars.hpp"
#include "FullyInclusiveSoftResummed/MellinFuncs.h"

using FullyInclusiveSoftResummed::cdcmplx;

// cast MMa as Macro
#define Power(x,y) pow(x,y)
// Mellin[beta^v](N) = EulerBeta(N,1+v/2)
#define MellinBeta sqrt(M_PI)/2. * exp(lngamma(N) - lngamma(N+1.5)) // Mellin[beta](N) = EulerBeta(N,3/2)
#define MellinBeta3 sqrt(M_PI)/2.*(3./2.) * exp(lngamma(N) - lngamma(N+2.5)); // Mellin[beta^3](N) = EulerBeta(N,5/2)
//#define MellinBeta5 sqrt(M_PI)/2.*(3./2.)*(5./2.) * exp(lngamma(N) - lngamma(N+3.5)); // Mellin[beta^5](N) = EulerBeta(N,7/2)

cdcmplx FullyInclusiveSoftResummed::coeffs::cg0t_F2_VV(cdbl m2, cdbl q2, cdcmplx N) {
    partonic_rhoq(m2,q2)
    return (M_PI)/(2.) * sqrt(-rhoq/(1 - rhoq)) * MellinBeta;
}

cdcmplx FullyInclusiveSoftResummed::coeffs::cg0t_F2_AA(cdbl m2, cdbl q2, cdcmplx N) {
    partonic_rhoq(m2,q2)
    return (1. - 2*rhoq) * cg0t_F2_VV(m2,q2,N);
}

cdcmplx FullyInclusiveSoftResummed::coeffs::cg0t_FL_VV(cdbl m2, cdbl q2, cdcmplx N) {
    partonic_rhoq(m2,q2)
    return (4*M_PI)/(3.) * sqrt((-rhoq)/Power(1 - rhoq,3)) * MellinBeta3;
    //partonic_betaq(m2,q2)
    //cdcmplx lo = (4*M_PI)/(3.) * sqrt((-rhoq)/Power(1 - rhoq,3)) * MellinBeta3;
    //cdcmplx nlo = 2.*M_PI/5.*(5.-6.*pow(betaq,2))/pow(betaq,3)/sqrt(-rhoq)*MellinBeta5;
    //return lo + nlo;
}

cdcmplx FullyInclusiveSoftResummed::coeffs::cg0t_FL_AA(cdbl m2, cdbl q2, cdcmplx N) {
    partonic_rhoq(m2,q2)
    return M_PI * sqrt((Power(-rhoq,3))/(1 - rhoq)) * MellinBeta;
}

cdcmplx FullyInclusiveSoftResummed::coeffs::cg0t_x2g1_VV(cdbl m2, cdbl q2, cdcmplx N) {
    return cg0t_F2_VV(m2,q2,N);
    //partonic_betaq(m2,q2)
    //return cg0t_F2_VV(m2,q2,N); + M_PI/12.*(11. - 10.*pow(betaq,2))/betaq/sqrt(-rhoq) * MellinBeta3;
}

cdcmplx FullyInclusiveSoftResummed::coeffs::cg0t_x2g1_AA(cdbl m2, cdbl q2, cdcmplx N) {
    return cg0t_x2g1_VV(m2,q2,N);
}

cdcmplx FullyInclusiveSoftResummed::coeffs::cg0t_xF3_VA(cdbl m2, cdbl q2, cdcmplx N) { return 0.; }
cdcmplx FullyInclusiveSoftResummed::coeffs::cg0t_g4_VA(cdbl m2, cdbl q2, cdcmplx N) { return 0.; }
cdcmplx FullyInclusiveSoftResummed::coeffs::cg0t_gL_VA(cdbl m2, cdbl q2, cdcmplx N) { return 0.; }
