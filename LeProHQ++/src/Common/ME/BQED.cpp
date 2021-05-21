#include "BQED.h"

// cast MMa as Macro
#define Power(x,y) pow(x,y)

cdbl Common::ME::BQED_F2_VV(cdbl m2, cdbl q2, cdbl sp, cdbl t1) {
cdbl u1 = - sp - t1;
return (-6*Power(q2,2))/Power(sp,2) - (6*q2)/sp - (m2*q2)/Power(t1,2) + ((6*m2*q2 + Power(q2,2))/u1 + u1/2.)/t1 - (m2*q2)/Power(u1,2) - (2*Power(m2,2)*Power(sp,2))/(Power(t1,2)*Power(u1,2)) + ((2*m2 + q2)*sp)/(t1*u1) + t1/(2.*u1);
}

cdbl Common::ME::BQED_F2_AA(cdbl m2, cdbl q2, cdbl sp, cdbl t1) {
cdbl u1 = - sp - t1;
return (-6*Power(q2,2))/Power(sp,2) - (6*q2)/sp - (m2*q2)/Power(t1,2) + Power(sp,2)*((-2*Power(m2,2))/(Power(t1,2)*Power(u1,2)) - (2*m2)/(q2*t1*u1)) + ((6*m2*q2 + Power(q2,2))/u1 + u1/2.)/t1 - (m2*q2)/Power(u1,2) + ((2*m2 + q2)*sp)/(t1*u1) + t1/(2.*u1);
}

cdbl Common::ME::BQED_FL_VV(cdbl m2, cdbl q2, cdbl sp, cdbl t1) {
cdbl u1 = - sp - t1;
return (-4*Power(q2,2))/Power(sp,2) - (4*q2)/sp + (4*m2*q2)/(t1*u1);
}

cdbl Common::ME::BQED_FL_AA(cdbl m2, cdbl q2, cdbl sp, cdbl t1) {
cdbl u1 = - sp - t1;
return (-4*Power(q2,2))/Power(sp,2) - (4*q2)/sp + Power(sp,2)*((-4*Power(m2,2))/(Power(t1,2)*Power(u1,2)) - (2*m2)/(q2*t1*u1)) + (8*m2*q2)/(t1*u1) + (4*m2*sp)/(t1*u1);
}

cdbl Common::ME::BQED_x2g1_VV(cdbl m2, cdbl q2, cdbl sp, cdbl t1) {
cdbl u1 = - sp - t1;
return sp*(m2/Power(t1,2) + m2/Power(u1,2)) - t1/(2.*u1) - u1/(2.*t1) + (-((q2*t1)/u1) - (q2*u1)/t1)/sp;
}

cdbl Common::ME::BQED_x2g1_AA(cdbl m2, cdbl q2, cdbl sp, cdbl t1) {
cdbl u1 = - sp - t1;
return sp*(m2/Power(t1,2) + m2/Power(u1,2)) - t1/(2.*u1) - u1/(2.*t1) + (-((q2*t1)/u1) - (q2*u1)/t1)/sp;
}

cdbl Common::ME::BQED_xF3_VA(cdbl m2, cdbl q2, cdbl sp, cdbl t1) {
cdbl u1 = - sp - t1;
cdbl r = -(q2/t1) + (-(Power(q2,2)/t1) + Power(q2,2)/u1)/sp + sp*((-0.5 - (m2*q2)/Power(u1,2))/t1 + 1/(2.*u1) + (m2*q2)/(Power(t1,2)*u1)) + q2/u1;
return r;
}

cdbl Common::ME::BQED_g4_VA(cdbl m2, cdbl q2, cdbl sp, cdbl t1) {
cdbl u1 = - sp - t1;
return -(q2/t1) + sp*(-1/(2.*t1) + 1/(2.*u1)) + Power(sp,2)*(-(m2/(t1*Power(u1,2))) + m2/(Power(t1,2)*u1)) + q2/u1;
}

cdbl Common::ME::BQED_gL_VA(cdbl m2, cdbl q2, cdbl sp, cdbl t1) {
//cdbl u1 = - sp - t1;
return 0;
}
