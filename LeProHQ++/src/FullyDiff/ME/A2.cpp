/// auto-generated by build.py on 2018/09/07 15:42:55
#include "A2.h"

#define Power(a,b) pow(a,b)

#define init1 cdbl s5 = sp + tp + up;
#define init2 cdbl s = sp + q2;

cdbl FullyDiff::ME::A2_F2_VV(cdbl m2, cdbl q2, cdbl sp, cdbl t1, cdbl u1, cdbl tp, cdbl up) {
init1
init2
return ((q2 + 2*s)/Power(s,2) + (18*q2)/Power(sp,2))*tp + s5*((-6*q2*(-2*m2 + sp - 2*t1)*tp)/(Power(sp,2)*Power(up,2)) + ((-2/s - (18*q2)/Power(sp,2))*tp)/up) + (6*q2*Power(s5,2)*tp)/(Power(sp,2)*Power(up,2)) + ((q2 - (12*q2*(m2*sp + (sp - t1)*t1))/Power(sp,2))*tp)/Power(up,2) + (tp*((12*q2*(-3*m2 + sp - 2*t1))/Power(sp,2) - (2*(2*(m2 + t1) + u1))/s))/up + (tp*((-6*q2*(-6*m2 + sp - 2*t1))/Power(sp,2) + (-Power(q2,2) + 2*m2*(q2 + 2*s) + s*sp - Power(sp,2) + 2*s*t1 + s*u1 - sp*u1 + q2*(s + 2*t1 + u1))/Power(s,2)) + (2*q2*(m2*sp + (sp - 6*t1)*t1)*tp)/(sp*Power(up,2)) - (tp*(Power(q2,2)*Power(sp,2) - 4*m2*(6*q2*s*sp + Power(sp,3)) + Power(sp,2)*(Power(sp,2) + 2*sp*t1 + 4*Power(t1,2) - s*u1 + 3*sp*u1 + 4*t1*u1 + 2*Power(u1,2)) + q2*(2*s*(Power(sp,2) - 6*sp*t1 + 6*Power(t1,2)) - Power(sp,2)*(4*sp + 6*t1 + 3*u1))))/(s*Power(sp,2)*up) + ((-6*q2)/Power(sp,2) - (q2 + sp)/Power(s,2))*tp*up)/s5 + ((q2*tp*(-12*m2*Power(s,2) + sp*(Power(s,2) + 2*Power(sp,2) + q2*(s - t1) + 3*sp*t1 + 2*Power(t1,2) - 3*s*(sp + t1) - 4*s*u1 + 4*(sp + t1)*u1 + 2*Power(u1,2))))/(Power(s,2)*sp) + (2*q2*Power(t1,2)*tp)/Power(up,2) - (tp*(2*m2*s*(q2 + sp) + q2*t1*(q2 + 3*s - 3*sp - 4*(t1 + u1))))/(s*up) + ((-12*m2*q2)/Power(sp,2) - (2*m2*(q2 + sp))/Power(s,2))*tp*up)/Power(s5,2);
}

cdbl FullyDiff::ME::A2_F2_AA(cdbl m2, cdbl q2, cdbl sp, cdbl t1, cdbl u1, cdbl tp, cdbl up) {
init1
init2
return ((2*sp*(6*Power(s,2) - 6*s*sp + Power(sp,2)) + 3*q2*(6*Power(s,2) - 4*s*sp + Power(sp,2)))*tp)/(Power(s,2)*Power(sp,2)) + s5*((-6*q2*(-2*m2 + sp - 2*t1)*tp)/(Power(sp,2)*Power(up,2)) + (2*(3*q2*(-3*s + sp) + sp*(-3*s + 2*sp))*tp)/(s*Power(sp,2)*up)) + (6*q2*Power(s5,2)*tp)/(Power(sp,2)*Power(up,2)) + ((q2 - (12*q2*(m2*sp + (sp - t1)*t1))/Power(sp,2))*tp)/Power(up,2) + (tp*(-3*Power(q2,2)*sp + 4*m2*(3*q2*(-3*s + sp) + sp*(-3*s + 2*sp)) + 3*q2*(s*(5*sp - 8*t1) + sp*(-2*sp + 4*t1 + u1)) + sp*(3*s*(sp - 4*t1 - u1) + sp*(-3*sp + 8*t1 + u1))))/(s*Power(sp,2)*up) + ((tp*(2*Power(q2,2)*(3*s - sp)*sp + m2*(4*sp*(6*Power(s,2) - 6*s*sp + Power(sp,2)) + 6*q2*(6*Power(s,2) - 4*s*sp + Power(sp,2))) + q2*(12*Power(s,2)*(-sp + t1) + 3*s*sp*(4*sp - 5*t1 - 2*u1) + Power(sp,2)*(-2*sp + 6*t1 + 3*u1)) + sp*(s*sp*(6*sp - 17*t1 - 7*u1) + Power(sp,2)*(-2*sp + 4*t1 + u1) + Power(s,2)*(-4*sp + 15*t1 + 6*u1))))/(Power(s,2)*Power(sp,2)) + (2*q2*(m2*sp + (sp - 6*t1)*t1)*tp)/(sp*Power(up,2)) + (tp*(4*m2*sp*(6*q2*s - 3*q2*sp + 3*s*sp - 2*Power(sp,2)) - Power(q2,2)*sp*(sp + 9*t1) - sp*(Power(sp,3) + 6*s*t1*(2*t1 + u1) + Power(sp,2)*(5*t1 + 3*u1) - sp*(3*s*t1 + 8*Power(t1,2) + s*u1 + 2*t1*u1 - 2*Power(u1,2))) + q2*(s*(-2*Power(sp,2) + 21*sp*t1 - 12*Power(t1,2)) + sp*(4*Power(sp,2) + 3*sp*(-2*t1 + u1) + 6*t1*(2*t1 + u1)))))/(s*Power(sp,2)*up) - ((6*q2*Power(s,2) + 6*s*(-q2 + s)*sp + (3*q2 - 8*s)*Power(sp,2) + 3*Power(sp,3))*tp*up)/(Power(s,2)*Power(sp,2)))/s5 + ((tp*(-4*m2*(3*q2*Power(s,2) + 3*s*(-q2 + s)*sp + (q2 - 4*s)*Power(sp,2) + Power(sp,3)) + q2*(Power(s,2)*(sp - 9*t1) + q2*(-2*sp*t1 + s*(sp + 9*t1)) + s*sp*(-3*sp + 7*t1 - 4*u1) + 2*sp*(Power(sp,2) + Power(t1 + u1,2) + sp*(t1 + 2*u1)))))/(Power(s,2)*sp) + (2*q2*Power(t1,2)*tp)/Power(up,2) - (tp*(2*m2*s*sp*(q2 + sp) + q2*t1*(3*s*(sp - 2*t1) + q2*(sp + 6*t1) + sp*(-3*sp + 2*t1 - 4*u1))))/(s*sp*up) - (2*m2*(6*q2*Power(s,2) + 6*s*(-q2 + s)*sp + (3*q2 - 8*s)*Power(sp,2) + 3*Power(sp,3))*tp*up)/(Power(s,2)*Power(sp,2)))/Power(s5,2);
}

cdbl FullyDiff::ME::A2_FL_VV(cdbl m2, cdbl q2, cdbl sp, cdbl t1, cdbl u1, cdbl tp, cdbl up) {
init1

return (12*q2*tp)/Power(sp,2) + s5*((-4*q2*(-2*m2 + sp - 2*t1)*tp)/(Power(sp,2)*Power(up,2)) - (12*q2*tp)/(Power(sp,2)*up)) + (4*q2*Power(s5,2)*tp)/(Power(sp,2)*Power(up,2)) - (8*q2*(m2*sp + (sp - t1)*t1)*tp)/(Power(sp,2)*Power(up,2)) + (8*q2*(-3*m2 + sp - 2*t1)*tp)/(Power(sp,2)*up) + ((-4*q2*(-6*m2 + sp - 2*t1)*tp)/Power(sp,2) - (8*q2*Power(t1,2)*tp)/(sp*Power(up,2)) + (8*q2*(2*m2*sp + (sp - t1)*t1)*tp)/(Power(sp,2)*up) - (4*q2*tp*up)/Power(sp,2))/s5 + ((-8*m2*q2*tp)/sp - (8*m2*q2*tp*up)/Power(sp,2))/Power(s5,2);
}

cdbl FullyDiff::ME::A2_FL_AA(cdbl m2, cdbl q2, cdbl sp, cdbl t1, cdbl u1, cdbl tp, cdbl up) {
init1
init2
return (12*q2*tp)/Power(sp,2) + s5*((-4*q2*(-2*m2 + sp - 2*t1)*tp)/(Power(sp,2)*Power(up,2)) - (12*q2*tp)/(Power(sp,2)*up)) + (4*q2*Power(s5,2)*tp)/(Power(sp,2)*Power(up,2)) - (8*q2*(m2*sp + (sp - t1)*t1)*tp)/(Power(sp,2)*Power(up,2)) - (2*q2*(12*m2*s + sp*(q2 - 5*s + sp) + 8*s*t1)*tp)/(s*Power(sp,2)*up) + ((q2*(24*m2*Power(s,2) - sp*(8*Power(s,2) - 5*s*sp + Power(sp,2) + q2*(-4*s + sp)) + 8*Power(s,2)*t1)*tp)/(Power(s,2)*Power(sp,2)) - (8*q2*Power(t1,2)*tp)/(sp*Power(up,2)) - (2*q2*(-8*m2*s*sp + t1*(sp*(3*q2 - 7*s + 3*sp) + 4*s*t1))*tp)/(s*Power(sp,2)*up) - (4*q2*tp*up)/Power(sp,2))/s5 + ((-8*m2*q2*tp)/sp - (8*m2*q2*tp*up)/Power(sp,2))/Power(s5,2);
}

cdbl FullyDiff::ME::A2_x2g1_VV(cdbl m2, cdbl q2, cdbl sp, cdbl t1, cdbl u1, cdbl tp, cdbl up) {
init1
init2
return ((-2*s*sp + q2*(-4*s + sp))*tp)/(Power(s,2)*sp) + s5*((-2*q2*tp)/(sp*Power(up,2)) + (2*q2*tp)/(s*sp*up)) + (q2*(-4*m2 + sp - 4*t1)*tp)/(sp*Power(up,2)) + (tp*(4*m2*q2 - Power(q2,2) - (s + sp)*u1 + q2*(5*s - sp + 4*t1 + u1)))/(s*sp*up) + ((q2*tp*(Power(s,2)*sp - 3*s*Power(sp,2) + 2*Power(sp,3) - 4*m2*s*(q2 + sp) - 3*Power(s,2)*t1 + 3*Power(sp,2)*t1 + 2*sp*Power(t1,2) + q2*(s*sp + 3*s*t1 - sp*t1) + 4*sp*(-s + sp + t1)*u1 + 2*sp*Power(u1,2)))/(Power(s,2)*sp) + (2*q2*Power(t1,2)*tp)/Power(up,2) - (tp*(2*m2*s*sp*(q2 + sp) + q2*t1*(s*(3*sp - 2*t1) + q2*(sp + 2*t1) - sp*(3*sp + 2*t1 + 4*u1))))/(s*sp*up) + (m2*(4*s - 2*sp)*(q2 + sp)*tp*up)/(Power(s,2)*sp))/Power(s5,2) + (-((tp*(Power(s,2)*sp - 2*s*Power(sp,2) + Power(sp,3) + Power(q2,2)*(-2*s + sp) + m2*(8*q2*s - 2*q2*sp + 4*s*sp) - Power(s,2)*t1 + 3*s*sp*t1 - 2*Power(s,2)*u1 + s*sp*u1 + Power(sp,2)*u1 + q2*(4*Power(s,2) - sp*(2*t1 + u1) + s*(-4*sp + 5*t1 + 2*u1))))/(Power(s,2)*sp)) + (2*q2*(m2*sp + (sp - 2*t1)*t1)*tp)/(sp*Power(up,2)) + (tp*(4*m2*(2*q2*s - q2*sp + s*sp) - Power(q2,2)*(sp + 3*t1) + q2*(-2*s*sp + 4*Power(sp,2) + 7*s*t1 + 2*sp*t1 + 4*Power(t1,2) + 3*sp*u1) - sp*(Power(sp,2) - s*(t1 + u1) + 3*sp*(t1 + u1) + 2*u1*(2*t1 + u1))))/(s*sp*up) + ((2*s - sp)*(q2 + sp)*tp*up)/(Power(s,2)*sp))/s5;
}

cdbl FullyDiff::ME::A2_x2g1_AA(cdbl m2, cdbl q2, cdbl sp, cdbl t1, cdbl u1, cdbl tp, cdbl up) {
init1
init2
return ((-2*s*sp + q2*(-4*s + sp))*tp)/(Power(s,2)*sp) + s5*((-2*q2*tp)/(sp*Power(up,2)) + (2*q2*tp)/(s*sp*up)) + (q2*(-4*m2 + sp - 4*t1)*tp)/(sp*Power(up,2)) + (tp*(4*m2*q2 - Power(q2,2) - (s + sp)*u1 + q2*(5*s - sp + 4*t1 + u1)))/(s*sp*up) + (-((tp*(Power(s,2)*sp - 2*s*Power(sp,2) + Power(sp,3) + Power(q2,2)*(-2*s + sp) + m2*(8*q2*s - 2*q2*sp + 4*s*sp) - Power(s,2)*t1 + 3*s*sp*t1 - 2*Power(s,2)*u1 + s*sp*u1 + Power(sp,2)*u1 + q2*(4*Power(s,2) - sp*(2*t1 + u1) + s*(-4*sp + 5*t1 + 2*u1))))/(Power(s,2)*sp)) + (2*q2*(m2*sp + (sp - 2*t1)*t1)*tp)/(sp*Power(up,2)) + (tp*(4*m2*(2*q2*s - q2*sp + s*sp) - Power(q2,2)*(sp + 3*t1) + q2*(-2*s*sp + 4*Power(sp,2) + 7*s*t1 + 2*sp*t1 + 4*Power(t1,2) + 3*sp*u1) - sp*(Power(sp,2) - s*(t1 + u1) + 3*sp*(t1 + u1) + 2*u1*(2*t1 + u1))))/(s*sp*up) + ((2*s - sp)*(q2 + sp)*tp*up)/(Power(s,2)*sp))/s5 + ((q2*tp*(Power(s,2)*sp - 3*s*Power(sp,2) + 2*Power(sp,3) - 4*m2*s*(q2 + sp) - 3*Power(s,2)*t1 + 3*Power(sp,2)*t1 + 2*sp*Power(t1,2) + q2*(s*sp + 3*s*t1 - sp*t1) + 4*sp*(-s + sp + t1)*u1 + 2*sp*Power(u1,2)))/(Power(s,2)*sp) + (2*q2*Power(t1,2)*tp)/Power(up,2) - (tp*(2*m2*s*sp*(q2 + sp) + q2*t1*(s*(3*sp - 2*t1) + q2*(sp + 2*t1) - sp*(3*sp + 2*t1 + 4*u1))))/(s*sp*up) + (2*m2*(2*s - sp)*(q2 + sp)*tp*up)/(Power(s,2)*sp))/Power(s5,2);
}

cdbl FullyDiff::ME::A2_xF3_VA(cdbl m2, cdbl q2, cdbl sp, cdbl t1, cdbl u1, cdbl tp, cdbl up) {
init1
init2
return ((-2*s*sp + q2*(-4*s + sp))*tp)/(Power(s,2)*sp) + s5*((-2*q2*tp)/(sp*Power(up,2)) + (2*q2*tp)/(s*sp*up)) + (q2*(-4*m2 + sp - 4*t1)*tp)/(sp*Power(up,2)) + (tp*(4*m2*q2 - Power(q2,2) - (s + sp)*u1 + q2*(5*s - sp + 4*t1 + u1)))/(s*sp*up) + (-((tp*(Power(s,2)*sp - 2*s*Power(sp,2) + Power(sp,3) + Power(q2,2)*(-2*s + sp) + m2*(8*q2*s - 2*q2*sp + 4*s*sp) - Power(s,2)*t1 + 3*s*sp*t1 - 2*Power(s,2)*u1 + s*sp*u1 + Power(sp,2)*u1 + q2*(4*Power(s,2) - sp*(2*t1 + u1) + s*(-4*sp + 5*t1 + 2*u1))))/(Power(s,2)*sp)) + (2*q2*(m2*sp + (sp - 2*t1)*t1)*tp)/(sp*Power(up,2)) + (tp*(4*m2*(2*q2*s - q2*sp + s*sp) - Power(q2,2)*(sp + 3*t1) + q2*(-2*s*sp + 4*Power(sp,2) + 7*s*t1 + 2*sp*t1 + 4*Power(t1,2) + 3*sp*u1) - sp*(Power(sp,2) - s*(t1 + u1) + 3*sp*(t1 + u1) + 2*u1*(2*t1 + u1))))/(s*sp*up) + ((2*s - sp)*(q2 + sp)*tp*up)/(Power(s,2)*sp))/s5 + ((q2*tp*(Power(s,2)*sp - 3*s*Power(sp,2) + 2*Power(sp,3) - 4*m2*s*(q2 + sp) - 3*Power(s,2)*t1 + 3*Power(sp,2)*t1 + 2*sp*Power(t1,2) + q2*(s*sp + 3*s*t1 - sp*t1) + 4*sp*(-s + sp + t1)*u1 + 2*sp*Power(u1,2)))/(Power(s,2)*sp) + (2*q2*Power(t1,2)*tp)/Power(up,2) - (tp*(2*m2*s*sp*(q2 + sp) + q2*t1*(s*(3*sp - 2*t1) + q2*(sp + 2*t1) - sp*(3*sp + 2*t1 + 4*u1))))/(s*sp*up) + (2*m2*(2*s - sp)*(q2 + sp)*tp*up)/(Power(s,2)*sp))/Power(s5,2);
}

cdbl FullyDiff::ME::A2_g4_VA(cdbl m2, cdbl q2, cdbl sp, cdbl t1, cdbl u1, cdbl tp, cdbl up) {
init1
init2
return -(((2*sp*(6*Power(s,2) - 6*s*sp + Power(sp,2)) + 3*q2*(6*Power(s,2) - 4*s*sp + Power(sp,2)))*tp)/(Power(s,2)*Power(sp,2))) + s5*((6*q2*(-2*m2 + sp - 2*t1)*tp)/(Power(sp,2)*Power(up,2)) + (2*(9*q2*s - 3*q2*sp + 3*s*sp - 2*Power(sp,2))*tp)/(s*Power(sp,2)*up)) - (6*q2*Power(s5,2)*tp)/(Power(sp,2)*Power(up,2)) + (q2*(-1 + (12*(m2*sp + (sp - t1)*t1))/Power(sp,2))*tp)/Power(up,2) + (tp*(3*Power(q2,2)*sp + 4*m2*(9*q2*s - 3*q2*sp + 3*s*sp - 2*Power(sp,2)) - 3*q2*(s*(5*sp - 8*t1) + sp*(-2*sp + 4*t1 + u1)) + sp*(sp*(3*sp - 8*t1 - u1) + 3*s*(-sp + 4*t1 + u1))))/(s*Power(sp,2)*up) + (-((tp*(2*Power(q2,2)*(3*s - sp)*sp + m2*(4*sp*(6*Power(s,2) - 6*s*sp + Power(sp,2)) + 6*q2*(6*Power(s,2) - 4*s*sp + Power(sp,2))) + q2*(12*Power(s,2)*(-sp + t1) + 3*s*sp*(4*sp - 5*t1 - 2*u1) + Power(sp,2)*(-2*sp + 6*t1 + 3*u1)) + sp*(s*sp*(6*sp - 17*t1 - 7*u1) + Power(sp,2)*(-2*sp + 4*t1 + u1) + Power(s,2)*(-4*sp + 15*t1 + 6*u1))))/(Power(s,2)*Power(sp,2))) - (2*q2*(m2*sp + (sp - 6*t1)*t1)*tp)/(sp*Power(up,2)) + (tp*(4*m2*sp*(3*q2*(-2*s + sp) + sp*(-3*s + 2*sp)) + Power(q2,2)*sp*(sp + 9*t1) + sp*(Power(sp,3) + 6*s*t1*(2*t1 + u1) + Power(sp,2)*(5*t1 + 3*u1) - sp*(3*s*t1 + 8*Power(t1,2) + s*u1 + 2*t1*u1 - 2*Power(u1,2))) + q2*(s*(2*Power(sp,2) - 21*sp*t1 + 12*Power(t1,2)) - sp*(4*Power(sp,2) + 3*sp*(-2*t1 + u1) + 6*t1*(2*t1 + u1)))))/(s*Power(sp,2)*up) + ((6*q2*Power(s,2) + 6*s*(-q2 + s)*sp + (3*q2 - 8*s)*Power(sp,2) + 3*Power(sp,3))*tp*up)/(Power(s,2)*Power(sp,2)))/s5 + ((tp*(4*m2*(3*q2*Power(s,2) + 3*s*(-q2 + s)*sp + (q2 - 4*s)*Power(sp,2) + Power(sp,3)) - q2*(Power(s,2)*(sp - 9*t1) + q2*(-2*sp*t1 + s*(sp + 9*t1)) + s*sp*(-3*sp + 7*t1 - 4*u1) + 2*sp*(Power(sp,2) + Power(t1 + u1,2) + sp*(t1 + 2*u1)))))/(Power(s,2)*sp) - (2*q2*Power(t1,2)*tp)/Power(up,2) + (tp*(2*m2*s*sp*(q2 + sp) + q2*t1*(3*s*(sp - 2*t1) + q2*(sp + 6*t1) + sp*(-3*sp + 2*t1 - 4*u1))))/(s*sp*up) + (2*m2*(6*q2*Power(s,2) + 6*s*(-q2 + s)*sp + (3*q2 - 8*s)*Power(sp,2) + 3*Power(sp,3))*tp*up)/(Power(s,2)*Power(sp,2)))/Power(s5,2);
}

cdbl FullyDiff::ME::A2_gL_VA(cdbl m2, cdbl q2, cdbl sp, cdbl t1, cdbl u1, cdbl tp, cdbl up) {
init1
init2
return (-12*q2*tp)/Power(sp,2) + s5*((4*q2*(-2*m2 + sp - 2*t1)*tp)/(Power(sp,2)*Power(up,2)) + (12*q2*tp)/(Power(sp,2)*up)) - (4*q2*Power(s5,2)*tp)/(Power(sp,2)*Power(up,2)) + (8*q2*(m2*sp + (sp - t1)*t1)*tp)/(Power(sp,2)*Power(up,2)) + (2*q2*(12*m2*s + sp*(q2 - 5*s + sp) + 8*s*t1)*tp)/(s*Power(sp,2)*up) + ((q2*(-24*m2*Power(s,2) + sp*(8*Power(s,2) - 5*s*sp + Power(sp,2) + q2*(-4*s + sp)) - 8*Power(s,2)*t1)*tp)/(Power(s,2)*Power(sp,2)) + (8*q2*Power(t1,2)*tp)/(sp*Power(up,2)) + (2*q2*(-8*m2*s*sp + t1*(sp*(3*q2 - 7*s + 3*sp) + 4*s*t1))*tp)/(s*Power(sp,2)*up) + (4*q2*tp*up)/Power(sp,2))/s5 + ((8*m2*q2*tp)/sp + (8*m2*q2*tp*up)/Power(sp,2))/Power(s5,2);
}