#include "cg0.h"

#include "FullyInclusive/PartonicVars.hpp"

// cast MMa as Macro
#define Power(x,y) pow(x,y)
#define pi M_PI
#define ln log


cdbl FullyInclusive::coeffs::cg0_F2_VV(cdbl m2, cdbl q2, cdbl s) {
partonic_chi(m2,s) partonic_rhoq(m2,q2)
return (beta*pi*rho*rhoq*(Power(rho,2) + Power(rhoq,2) + rho*rhoq*(6 + rhoq)))/(2.*Power(rho - rhoq,3)) + (pi*rho*rhoq*(2*Power(rhoq,2) + 2*rho*Power(rhoq,2) + Power(rho,2)*(2 - (-4 + rhoq)*rhoq))*ln(chi))/(4.*Power(rho - rhoq,3));
}

cdbl FullyInclusive::coeffs::cg0_F2_AA(cdbl m2, cdbl q2, cdbl s) {
partonic_chi(m2,s) partonic_rhoq(m2,q2)
return (beta*pi*rho*rhoq*(Power(rho,2) + Power(rhoq,2) + rho*rhoq*(6 + rhoq)))/(2.*Power(rho - rhoq,3)) + (pi*rho*rhoq*(6*rho*Power(rhoq,2) - 2*(-1 + rhoq)*Power(rhoq,2) + Power(rho,2)*(2 - (-2 + rhoq)*rhoq))*ln(chi))/(4.*Power(rho - rhoq,3));
}

cdbl FullyInclusive::coeffs::cg0_FL_VV(cdbl m2, cdbl q2, cdbl s) {
partonic_chi(m2,s) partonic_rhoq(m2,q2)
return (2*beta*pi*Power(rho,2)*Power(rhoq,2))/Power(rho - rhoq,3) + (pi*Power(rho,3)*Power(rhoq,2)*ln(chi))/Power(rho - rhoq,3);
}

cdbl FullyInclusive::coeffs::cg0_FL_AA(cdbl m2, cdbl q2, cdbl s) {
partonic_chi(m2,s) partonic_rhoq(m2,q2)
return (beta*pi*Power(rho,2)*Power(rhoq,2)*(2 + rhoq))/Power(rho - rhoq,3) - (pi*rho*Power(rhoq,2)*(Power(rho,2)*(-1 + rhoq) - 4*rho*rhoq + Power(rhoq,2))*ln(chi))/(2.*Power(rho - rhoq,3));
}

cdbl FullyInclusive::coeffs::cg0_x2g1_VV(cdbl m2, cdbl q2, cdbl s) {
partonic_chi(m2,s) partonic_rhoq(m2,q2)
return (beta*pi*rho*rhoq*(rho + 3*rhoq))/(2.*Power(rho - rhoq,2)) + (pi*rho*rhoq*(rho + rhoq)*ln(chi))/(2.*Power(rho - rhoq,2));
}

cdbl FullyInclusive::coeffs::cg0_x2g1_AA(cdbl m2, cdbl q2, cdbl s) {
partonic_chi(m2,s) partonic_rhoq(m2,q2)
return (beta*pi*rho*rhoq*(rho + 3*rhoq))/(2.*Power(rho - rhoq,2)) + (pi*rho*rhoq*(rho + rhoq)*ln(chi))/(2.*Power(rho - rhoq,2));
}

cdbl FullyInclusive::coeffs::cg0_xF3_VA(cdbl m2, cdbl q2, cdbl s) {

return 0;
}

cdbl FullyInclusive::coeffs::cg0_g4_VA(cdbl m2, cdbl q2, cdbl s) {

return 0;
}

cdbl FullyInclusive::coeffs::cg0_gL_VA(cdbl m2, cdbl q2, cdbl s) {

return 0;
}
