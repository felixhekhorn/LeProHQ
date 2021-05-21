#include "cg0t.h"

#include "FullyInclusive/PartonicVars.hpp"

// cast MMa as Macro
#define Power(x,y) pow(x,y)
#define pi M_PI
#define ln log


cdbl FullyInclusive::coeffs::cg0t_F2_VV(cdbl m2, cdbl q2, cdbl s) {
partonic_beta(m2,s) partonic_rhoq(m2,q2)
return (beta*pi*rhoq)/(2.*(-1 + rhoq));
}

cdbl FullyInclusive::coeffs::cg0t_F2_AA(cdbl m2, cdbl q2, cdbl s) {
partonic_beta(m2,s) partonic_rhoq(m2,q2)
return (beta*pi*(1 - 2*rhoq)*rhoq)/(2.*(-1 + rhoq));
}

cdbl FullyInclusive::coeffs::cg0t_FL_VV(cdbl m2, cdbl q2, cdbl s) {
partonic_beta(m2,s) partonic_rhoq(m2,q2)
return (-4*Power(beta,3)*pi*Power(rhoq,2))/(3.*Power(-1 + rhoq,3));
}

cdbl FullyInclusive::coeffs::cg0t_FL_AA(cdbl m2, cdbl q2, cdbl s) {
partonic_beta(m2,s) partonic_rhoq(m2,q2)
return -((beta*pi*Power(rhoq,2))/(-1 + rhoq));
}

cdbl FullyInclusive::coeffs::cg0t_x2g1_VV(cdbl m2, cdbl q2, cdbl s) {
partonic_beta(m2,s) partonic_rhoq(m2,q2)
return (beta*pi*rhoq)/(2.*(-1 + rhoq));
}

cdbl FullyInclusive::coeffs::cg0t_x2g1_AA(cdbl m2, cdbl q2, cdbl s) {
partonic_beta(m2,s) partonic_rhoq(m2,q2)
return (beta*pi*rhoq)/(2.*(-1 + rhoq));
}

cdbl FullyInclusive::coeffs::cg0t_xF3_VA(cdbl m2, cdbl q2, cdbl s) {

return 0;
}

cdbl FullyInclusive::coeffs::cg0t_g4_VA(cdbl m2, cdbl q2, cdbl s) {

return 0;
}

cdbl FullyInclusive::coeffs::cg0t_gL_VA(cdbl m2, cdbl q2, cdbl s) {

return 0;
}
