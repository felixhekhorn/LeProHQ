#include "dq1.h"

#include <gsl/gsl_sf_dilog.h>

#include "FullyInclusive/PartonicVars.hpp"

// cast MMa as Macro
#define Power(x,y) pow(x,y)
#define pi M_PI
#define ln log
#define Li2 gsl_sf_dilog

cdbl FullyInclusive::coeffs::dq1_h4(cdbl m2, cdbl q2, cdbl s) {
partonic_chi(m2,s) partonic_chip(m2,q2,s)
return -Li2(((1 + chi)*chip)/(1 + chip)) + Li2(((1 + chi)*chip)/(chi*(1 + chip))) + Li2((1 + chip)/(1 + chi)) - Li2((chi*(1 + chip))/(1 + chi)) + Power(ln(chi),2)/2. + ln(chi)*(ln(1 + chi) - ln(chi - chip) + ln(1 + chip) - ln(1 - chi*chip));
}


cdbl FullyInclusive::coeffs::dq1_F2_VV(cdbl m2, cdbl q2, cdbl s) {
partonic_chi(m2,s) partonic_chip(m2,q2,s)
cdbl h4 = FullyInclusive::coeffs::dq1_h4(m2,q2,s);
return (beta*(2*Power(rho,2)*(718 + 5*rho) - 4*rho*(758 + 91*rho)*rhop + (2488 + 5108*rho)*Power(rhop,2) - 6456*Power(rhop,3)) + h4*(-288*Power(rho,2) + 288*rho*rhop + 36*(-4 + 3*Power(rho,2))*Power(rhop,2) - 972*rho*Power(rhop,3) + 972*Power(rhop,4)) + (27*Power(rho,2)*(-8 + Power(rho,2)) + 18*rho*(8 + 3*Power(rho,2))*rhop + 486*Power(rho,2)*Power(rhop,2) - 972*rho*Power(rhop,3))*ln(chi) + betap*(912*Power(rho,2) - 24*rho*(56 + 23*rho)*rhop + 48*(17 + 52*rho)*Power(rhop,2) - 2256*Power(rhop,3))*ln((chi - chip)/(1 - chi*chip)))/(5184.*pi*rho);
}

cdbl FullyInclusive::coeffs::dq1_FL_VV(cdbl m2, cdbl q2, cdbl s) {
partonic_chi(m2,s) partonic_chip(m2,q2,s)
cdbl h4 = FullyInclusive::coeffs::dq1_h4(m2,q2,s);
return (rhop*(beta*(4*rho*(-38 + 23*rho) + (200 + 532*rho)*rhop - 744*Power(rhop,2)) + h4*(-108*rho*Power(rhop,2) + 108*Power(rhop,3)) + (18*Power(rho,3) + 54*Power(rho,2)*rhop - 108*rho*Power(rhop,2))*ln(chi) + betap*(-48*rho + (48 + 264*rho)*rhop - 264*Power(rhop,2))*ln((chi - chip)/(1 - chi*chip))))/(864.*pi*rho);
}

cdbl FullyInclusive::coeffs::dq1_x2g1_VV(cdbl m2, cdbl q2, cdbl s) {
partonic_chi(m2,s) partonic_chip(m2,q2,s)
cdbl h4 = FullyInclusive::coeffs::dq1_h4(m2,q2,s);
return (beta*(2*Power(rho,2)*(718 + 5*rho) - 4*rho*(530 + 229*rho)*rhop + 8*(218 + 205*rho)*Power(rhop,2) - 1200*Power(rhop,3)) + h4*(-288*Power(rho,2) + 288*rho*rhop + 36*(-4 + 3*Power(rho,2))*Power(rhop,2) - 108*rho*Power(rhop,3)) + (27*Power(rho,2)*(-8 + Power(rho,2)) + 9*rho*(16 - 6*Power(rho,2))*rhop + 108*Power(rho,2)*Power(rhop,2))*ln(chi) + betap*(912*Power(rho,2) - 24*rho*(44 + 23*rho)*rhop + (672 + 912*rho)*Power(rhop,2) - 600*Power(rhop,3))*ln((chi - chip)/(1 - chi*chip)))/(5184.*pi*rho);
}

cdbl FullyInclusive::coeffs::dq1_F2_AA(cdbl m2, cdbl q2, cdbl s) {
return FullyInclusive::coeffs::dq1_F2_VV(m2,q2,s);
}

cdbl FullyInclusive::coeffs::dq1_FL_AA(cdbl m2, cdbl q2, cdbl s) {
return FullyInclusive::coeffs::dq1_FL_VV(m2,q2,s);
}

cdbl FullyInclusive::coeffs::dq1_x2g1_AA(cdbl m2, cdbl q2, cdbl s) {
return FullyInclusive::coeffs::dq1_x2g1_VV(m2,q2,s);
}

cdbl FullyInclusive::coeffs::dq1_xF3_VA(cdbl m2, cdbl q2, cdbl s) {
return FullyInclusive::coeffs::dq1_x2g1_VV(m2,q2,s);
}

cdbl FullyInclusive::coeffs::dq1_g4_VA(cdbl m2, cdbl q2, cdbl s) {
return -FullyInclusive::coeffs::dq1_F2_VV(m2,q2,s);
}

cdbl FullyInclusive::coeffs::dq1_gL_VA(cdbl m2, cdbl q2, cdbl s) {
return -FullyInclusive::coeffs::dq1_FL_VV(m2,q2,s);
}
