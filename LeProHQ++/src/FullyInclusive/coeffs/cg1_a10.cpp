#include "cg1_a10.h"

#include <gsl/gsl_sf_dilog.h>

#include "FullyInclusive/PartonicVars.hpp"

// cast MMa as Macro
#define Power(x,y) pow(x,y)
#define pi M_PI
#define ln log
#define Li2 gsl_sf_dilog


cdbl FullyInclusive::coeffs::cg1_a10_g1(cdbl m2, cdbl q2) {
partonic_chiq(m2,q2)
return (-Power(pi,2)/2. + Li2(-(rhoq/(-2 + rhoq))) - (2*(1 - rhoq)*ln(chiq))/betaq - (3*Power(ln(chiq),2))/2. + Power(ln(rhoq/(-2 + rhoq)),2)/2.)/8.;
}

cdbl FullyInclusive::coeffs::cg1_a10_g2(cdbl m2, cdbl q2) {
partonic_chiq(m2,q2)
return (12.5 - 15*ln(2) + 9*Power(ln(2),2) + Power(ln(chiq),2) - Power(ln(rhoq/(2.*(-1 + rhoq))),2))/4.;
}


cdbl FullyInclusive::coeffs::cg1_a10_F2_VV_OK(cdbl m2, cdbl q2) {
partonic_chiq(m2,q2)
cdbl g1 = FullyInclusive::coeffs::cg1_a10_g1(m2,q2);
cdbl g2 = FullyInclusive::coeffs::cg1_a10_g2(m2,q2);
return g1 + g2 - (Power(pi,2)*rhoq)/(32.*(-1 + rhoq)) + ((-5 + (7 - 2*rhoq)*rhoq)*ln(chiq))/(8.*betaq*(-1 + rhoq)) - (rhoq*Power(ln(chiq),2))/(32.*(-1 + rhoq)) - ((-3 + rhoq)*ln(rhoq/(2.*(-1 + rhoq))))/(4.*(-2 + rhoq));
}

cdbl FullyInclusive::coeffs::cg1_a10_F2_VV_QED(cdbl m2, cdbl q2) {
partonic_chiq(m2,q2)
cdbl g1 = FullyInclusive::coeffs::cg1_a10_g1(m2,q2);

return -(g1/(-1 + rhoq)) + (Power(pi,2)*(-4 + 3*rhoq))/(96.*(-1 + rhoq)) - (-9 + 5*rhoq)/(8.*(-2 + rhoq)) - ln(chiq)/(8.*betaq) + ((-4 + rhoq)*Power(ln(chiq),2))/(32.*(-1 + rhoq)) + ((-3 + (5 - 2*rhoq)*rhoq)*ln(rhoq/(2.*(-1 + rhoq))))/(4.*Power(-2 + rhoq,2)*(-1 + rhoq));
}

cdbl FullyInclusive::coeffs::cg1_a10_FL_VV_OK(cdbl m2, cdbl q2) {
partonic_chiq(m2,q2)
cdbl g1 = FullyInclusive::coeffs::cg1_a10_g1(m2,q2);
cdbl g2 = FullyInclusive::coeffs::cg1_a10_g2(m2,q2);
return 0.6805555555555556 + g2 - (Power(pi,2)*(-4 + rhoq)*rhoq)/(24.*Power(-1 + rhoq,2)) + (g1*(1 + 2*rhoq))/Power(-1 + rhoq,2) - ln(2) - ((-4 + rhoq)*rhoq*Power(ln(chiq),2))/(8.*Power(-1 + rhoq,2)) + ((-1 + rhoq*(-3 - 2*(-3 + rhoq)*rhoq))*ln(rhoq/(2.*(-1 + rhoq))))/(4.*(-2 + rhoq)*Power(-1 + rhoq,2));
}

cdbl FullyInclusive::coeffs::cg1_a10_FL_VV_QED(cdbl m2, cdbl q2) {
partonic_chiq(m2,q2)
cdbl g1 = FullyInclusive::coeffs::cg1_a10_g1(m2,q2);

return (3 - 2*rhoq)/(8.*(-2 + rhoq)) - (g1*(-1 + 6*rhoq))/Power(-1 + rhoq,2) - (Power(pi,2)*(-1 + 6*rhoq))/(24.*Power(-1 + rhoq,2)) + ((-6 + rhoq + Power(rhoq,2))*ln(chiq))/(8.*betaq*(-2 + rhoq)) - ((-1 + 6*rhoq)*Power(ln(chiq),2))/(8.*Power(-1 + rhoq,2)) + ((3 + 2*rhoq*(5 + (-5 + rhoq)*rhoq))*ln(rhoq/(2.*(-1 + rhoq))))/(4.*Power(-2 + rhoq,2)*(-1 + rhoq));
}

cdbl FullyInclusive::coeffs::cg1_a10_x2g1_VV_OK(cdbl m2, cdbl q2) {
partonic_chiq(m2,q2)
cdbl g1 = FullyInclusive::coeffs::cg1_a10_g1(m2,q2);
cdbl g2 = FullyInclusive::coeffs::cg1_a10_g2(m2,q2);
return g1 + g2 - (Power(pi,2)*rhoq)/(32.*(-1 + rhoq)) + ((-5 + (7 - 2*rhoq)*rhoq)*ln(chiq))/(8.*betaq*(-1 + rhoq)) - (rhoq*Power(ln(chiq),2))/(32.*(-1 + rhoq)) - ((-3 + rhoq)*ln(rhoq/(2.*(-1 + rhoq))))/(4.*(-2 + rhoq));
}

cdbl FullyInclusive::coeffs::cg1_a10_x2g1_VV_QED(cdbl m2, cdbl q2) {
partonic_chiq(m2,q2)
cdbl g1 = FullyInclusive::coeffs::cg1_a10_g1(m2,q2);

return -(g1/(-1 + rhoq)) + (Power(pi,2)*(-4 + 3*rhoq))/(96.*(-1 + rhoq)) - (-9 + 5*rhoq)/(8.*(-2 + rhoq)) - ln(chiq)/(8.*betaq) + ((-4 + rhoq)*Power(ln(chiq),2))/(32.*(-1 + rhoq)) + ((-3 + (5 - 2*rhoq)*rhoq)*ln(rhoq/(2.*(-1 + rhoq))))/(4.*Power(-2 + rhoq,2)*(-1 + rhoq));
}

cdbl FullyInclusive::coeffs::cg1_a10_F2_AA_OK(cdbl m2, cdbl q2) {
partonic_chiq(m2,q2)
cdbl g1 = FullyInclusive::coeffs::cg1_a10_g1(m2,q2);
cdbl g2 = FullyInclusive::coeffs::cg1_a10_g2(m2,q2);
return g2 - (Power(pi,2)*rhoq*(3 + rhoq*(-25 + 18*rhoq)))/(96.*Power(-1 + rhoq,2)*(-1 + 2*rhoq)) + (g1*(-1 - rhoq*(-4 + rhoq + Power(rhoq,2))))/(Power(-1 + rhoq,2)*(-1 + 2*rhoq)) + ((-3 + 2*rhoq*(2 + rhoq))*ln(chiq))/(8.*betaq*(-1 + 2*rhoq)) - (rhoq*(1 + rhoq*(-19 + 14*rhoq))*Power(ln(chiq),2))/(32.*Power(-1 + rhoq,2)*(-1 + 2*rhoq)) - ((-1 + (-1 + rhoq)*Power(rhoq,2))*ln(rhoq/(2.*(-1 + rhoq))))/(4.*(-2 + rhoq)*(-1 + rhoq)*(-1 + 2*rhoq));
}

cdbl FullyInclusive::coeffs::cg1_a10_F2_AA_QED(cdbl m2, cdbl q2) {
partonic_chiq(m2,q2)
cdbl g1 = FullyInclusive::coeffs::cg1_a10_g1(m2,q2);

return (-5 - 4*(-3 + rhoq)*rhoq)/(8.*(-2 + rhoq)*(-1 + 2*rhoq)) + (g1*(-1 + (-2 + rhoq)*rhoq*(-2 + 3*rhoq)))/(Power(-1 + rhoq,2)*(-1 + 2*rhoq)) + (Power(pi,2)*(-4 + rhoq*(19 + rhoq*(-41 + 18*rhoq))))/(96.*Power(-1 + rhoq,2)*(-1 + 2*rhoq)) - ((1 + 6*(-1 + rhoq)*rhoq)*ln(chiq))/(8.*betaq*(-1 + 2*rhoq)) + ((-4 + rhoq*(17 + 7*rhoq*(-5 + 2*rhoq)))*Power(ln(chiq),2))/(32.*Power(-1 + rhoq,2)*(-1 + 2*rhoq)) + ((1 + rhoq)*(-1 + rhoq*(7 + rhoq*(-7 + 2*rhoq)))*ln(rhoq/(2.*(-1 + rhoq))))/(4.*Power(-2 + rhoq,2)*(-1 + rhoq)*(-1 + 2*rhoq));
}

cdbl FullyInclusive::coeffs::cg1_a10_FL_AA_OK(cdbl m2, cdbl q2) {
partonic_chiq(m2,q2)
cdbl g1 = FullyInclusive::coeffs::cg1_a10_g1(m2,q2);
cdbl g2 = FullyInclusive::coeffs::cg1_a10_g2(m2,q2);
return g2 + (Power(pi,2)*(4 - 9*(-1 + rhoq)*rhoq))/(96.*Power(-1 + rhoq,2)) + (g1*(4 - rhoq*(1 + rhoq)))/(2.*Power(-1 + rhoq,2)) + ((2 - 3*rhoq + Power(rhoq,3))*ln(chiq))/(8.*betaq*Power(-1 + rhoq,2)) + ((4 - 7*(-1 + rhoq)*rhoq)*Power(ln(chiq),2))/(32.*Power(-1 + rhoq,2)) + ((-4 + 5*rhoq - Power(rhoq,3))*ln(rhoq/(2.*(-1 + rhoq))))/(8.*(-2 + rhoq)*Power(-1 + rhoq,2));
}

cdbl FullyInclusive::coeffs::cg1_a10_FL_AA_QED(cdbl m2, cdbl q2) {
partonic_chiq(m2,q2)
cdbl g1 = FullyInclusive::coeffs::cg1_a10_g1(m2,q2);

return -(-5 + 2*rhoq)/(8.*(-2 + rhoq)) + (g1*rhoq*(-7 + 3*rhoq))/(2.*Power(-1 + rhoq,2)) + (Power(pi,2)*rhoq*(-17 + 9*rhoq))/(96.*Power(-1 + rhoq,2)) + ((-0.5 + rhoq - (3*Power(rhoq,2))/8.)*ln(chiq))/(betaq*(-2 + rhoq)) + (rhoq*(-15 + 7*rhoq)*Power(ln(chiq),2))/(32.*Power(-1 + rhoq,2)) + ((12 + rhoq*(-7 + rhoq*(-3 + 2*rhoq)))*ln(rhoq/(2.*(-1 + rhoq))))/(8.*Power(-2 + rhoq,2)*(-1 + rhoq));
}

cdbl FullyInclusive::coeffs::cg1_a10_x2g1_AA_OK(cdbl m2, cdbl q2) {
partonic_chiq(m2,q2)
cdbl g1 = FullyInclusive::coeffs::cg1_a10_g1(m2,q2);
cdbl g2 = FullyInclusive::coeffs::cg1_a10_g2(m2,q2);
return g2 + g1/Power(-1 + rhoq,2) + (Power(pi,2)*(11 - 7*rhoq)*rhoq)/(96.*Power(-1 + rhoq,2)) + (3*ln(chiq))/(8.*betaq) + ((9 - 5*rhoq)*rhoq*Power(ln(chiq),2))/(32.*Power(-1 + rhoq,2)) + ((1 + rhoq*(-5 - 2*(-3 + rhoq)*rhoq))*ln(rhoq/(2.*(-1 + rhoq))))/(4.*(-2 + rhoq)*Power(-1 + rhoq,2));
}

cdbl FullyInclusive::coeffs::cg1_a10_x2g1_AA_QED(cdbl m2, cdbl q2) {
partonic_chiq(m2,q2)
cdbl g1 = FullyInclusive::coeffs::cg1_a10_g1(m2,q2);

return (5 - 2*rhoq)/(8.*(-2 + rhoq)) + (g1*(1 + (-4 + rhoq)*rhoq))/Power(-1 + rhoq,2) + (Power(pi,2)*(4 + rhoq*(-19 + 7*rhoq)))/(96.*Power(-1 + rhoq,2)) + ((-2 + (5 - 2*rhoq)*rhoq)*ln(chiq))/(8.*betaq*(-2 + rhoq)) + ((4 + rhoq*(-17 + 5*rhoq))*Power(ln(chiq),2))/(32.*Power(-1 + rhoq,2)) + ((1 + (-2 + rhoq)*rhoq*(-3 + 2*rhoq))*ln(rhoq/(2.*(-1 + rhoq))))/(4.*Power(-2 + rhoq,2)*(-1 + rhoq));
}
