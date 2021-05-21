#include "cqBarF1.h"

#include <gsl/gsl_sf_dilog.h>

#include "FullyInclusive/PartonicVars.hpp"
#include "cgBar1.h"

// cast MMa as Macro
#define Power(x,y) pow(x,y)
#define pi M_PI
#define ln log
#define Li2 gsl_sf_dilog



cdbl FullyInclusive::coeffs::cqBarF1_F2_VV(cdbl m2, cdbl q2, cdbl s) {
partonic_chi(m2,s) partonic_chiq(m2,q2)
cdbl h1 = FullyInclusive::coeffs::cgBar1_h1(m2,q2,s);
return (256*chi*chiq*(3*h1*rho*(rho - rhoq)*(-1 + rhoq)*(-2*rhoq + rho*(4 + rhoq)) - 2*beta*(-1 + rhoq)*rhoq*(-7*rho*rhoq - Power(rhoq,2) + Power(rho,2)*(16 + rhoq)) + rho*(-1 + rhoq)*(-21*rho*rhoq + 3*Power(rhoq,2) + Power(rho,2)*(14 + (-6 + rhoq)*rhoq))*ln(chi) + betaq*(-rho + rhoq)*(-2*rho*rhoq*(4 + rhoq) + Power(rhoq,2)*(4 + rhoq) + Power(rho,2)*(-14 + 19*rhoq))*ln((chi + chiq)/(1 + chi*chiq))))/(9.*(-1 + beta)*Power(1 + beta,5)*Power(1 + betaq,8)*Power(chi + chiq,3)*Power(1 + chi*chiq,3)*Power(1 - Power(chiq,2),2)*pi);
}

cdbl FullyInclusive::coeffs::cqBarF1_FL_VV(cdbl m2, cdbl q2, cdbl s) {
partonic_chi(m2,s) partonic_chiq(m2,q2)

return (256*chi*chiq*(-2*beta*(-1 + rhoq)*rhoq*(2*rho*Power(rhoq,2) + (2 - 3*rhoq)*Power(rhoq,2) + Power(rho,2)*(-6 + 5*rhoq)) - 4*Power(rho,2)*Power(-1 + rhoq,2)*(rho*(-3 + rhoq) + 3*rhoq)*ln(chi) + betaq*(-rho + rhoq)*(2*rho*Power(rhoq,2) + Power(rhoq,3)*(-4 + 3*rhoq) + Power(rho,2)*(12 + rhoq*(-22 + 9*rhoq)))*ln((chi + chiq)/(1 + chi*chiq))))/(9.*(-1 + beta)*Power(1 + beta,5)*Power(1 + betaq,8)*Power(chi + chiq,3)*Power(1 + chi*chiq,3)*Power(1 - Power(chiq,2),2)*pi*(-1 + rhoq));
}

cdbl FullyInclusive::coeffs::cqBarF1_x2g1_VV(cdbl m2, cdbl q2, cdbl s) {
partonic_chi(m2,s) partonic_chiq(m2,q2)
cdbl h1 = FullyInclusive::coeffs::cgBar1_h1(m2,q2,s);
return (256*chi*chiq*(-36*beta*rho*(rho - rhoq)*Power(-1 + rhoq,2)*rhoq + 6*h1*rho*Power(-1 + rhoq,2)*(2*Power(rho,2) - 3*rho*rhoq + Power(rhoq,2)) - 3*rho*(rho - rhoq)*Power(-1 + rhoq,2)*(rho*(-6 + rhoq) + 7*rhoq)*ln(chi) - 6*betaq*rho*Power(-1 + rhoq,2)*(3*Power(rho,2) - 4*rho*rhoq + Power(rhoq,2))*ln((chi + chiq)/(1 + chi*chiq))))/(9.*(-1 + beta)*Power(1 + beta,5)*Power(1 + betaq,8)*Power(chi + chiq,3)*Power(1 + chi*chiq,3)*Power(1 - Power(chiq,2),2)*pi*(-1 + rhoq));
}

cdbl FullyInclusive::coeffs::cqBarF1_F2_AA(cdbl m2, cdbl q2, cdbl s) {
partonic_chi(m2,s) partonic_chiq(m2,q2)
cdbl h1 = FullyInclusive::coeffs::cgBar1_h1(m2,q2,s);
return (256*chi*chiq*(-3*h1*rho*(rho - rhoq)*(-1 + rhoq)*(rho*(-4 + rhoq) - 2*(-1 + rhoq)*rhoq) + beta*(-1 + rhoq)*rhoq*(Power(rho,2)*(-32 + rhoq) + rho*(14 - 3*rhoq)*rhoq + 2*Power(rhoq,2)) + (rho*(-1 + rhoq)*(6*Power(rhoq,2)*(1 + rhoq) + 3*rho*rhoq*(-14 + (-6 + rhoq)*rhoq) - Power(rho,2)*(-28 + Power(rhoq,2)))*ln(chi))/2. + betaq*(Power(rho,3)*(14 - 13*rhoq) + 3*Power(rho,2)*(-2 + rhoq)*rhoq + (4 - 5*rhoq)*Power(rhoq,3) + 3*rho*Power(rhoq,2)*(-4 + 5*rhoq))*ln((chi + chiq)/(1 + chi*chiq))))/(9.*(-1 + beta)*Power(1 + beta,5)*Power(1 + betaq,8)*Power(chi + chiq,3)*Power(1 + chi*chiq,3)*Power(1 - Power(chiq,2),2)*pi);
}

cdbl FullyInclusive::coeffs::cqBarF1_FL_AA(cdbl m2, cdbl q2, cdbl s) {
partonic_chi(m2,s) partonic_chiq(m2,q2)
cdbl h1 = FullyInclusive::coeffs::cgBar1_h1(m2,q2,s);
return (128*chi*chiq*(12*h1*rho*(rho - rhoq)*(-1 + rhoq)*Power(rhoq,2) - 2*beta*(-1 + rhoq)*rhoq*(-4*Power(rhoq,2) + 3*rho*Power(rhoq,2) + Power(rho,2)*(12 + rhoq)) + rho*(-1 + rhoq)*(6*Power(rhoq,3) + 3*rho*rhoq*(-8 + (-6 + rhoq)*rhoq) + Power(rho,2)*(24 + (-4 + rhoq)*rhoq))*ln(chi) - 8*betaq*(3*Power(rho,2)*rhoq - 3*rho*Power(rhoq,3) + Power(rhoq,4) + Power(rho,3)*(-3 + 2*rhoq))*ln((chi + chiq)/(1 + chi*chiq))))/(9.*(-1 + beta)*Power(1 + beta,5)*Power(1 + betaq,8)*Power(chi + chiq,3)*Power(1 + chi*chiq,3)*Power(1 - Power(chiq,2),2)*pi);
}

cdbl FullyInclusive::coeffs::cqBarF1_x2g1_AA(cdbl m2, cdbl q2, cdbl s) {
partonic_chi(m2,s) partonic_chiq(m2,q2)
cdbl h1 = FullyInclusive::coeffs::cgBar1_h1(m2,q2,s);
return (256*chi*chiq*(-36*beta*rho*(rho - rhoq)*Power(-1 + rhoq,2)*rhoq + 6*h1*rho*Power(-1 + rhoq,2)*(2*Power(rho,2) - 3*rho*rhoq + Power(rhoq,2)) - 3*rho*(rho - rhoq)*Power(-1 + rhoq,2)*(rho*(-6 + rhoq) + 7*rhoq)*ln(chi) - 6*betaq*rho*Power(-1 + rhoq,2)*(3*Power(rho,2) - 4*rho*rhoq + Power(rhoq,2))*ln((chi + chiq)/(1 + chi*chiq))))/(9.*(-1 + beta)*Power(1 + beta,5)*Power(1 + betaq,8)*Power(chi + chiq,3)*Power(1 + chi*chiq,3)*Power(1 - Power(chiq,2),2)*pi*(-1 + rhoq));
}

cdbl FullyInclusive::coeffs::cqBarF1_xF3_VA(cdbl m2, cdbl q2, cdbl s) {
return 0.;
}

cdbl FullyInclusive::coeffs::cqBarF1_g4_VA(cdbl m2, cdbl q2, cdbl s) {
return 0.;
}

cdbl FullyInclusive::coeffs::cqBarF1_gL_VA(cdbl m2, cdbl q2, cdbl s) {
return 0.;
}
