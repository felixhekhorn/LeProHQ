#include "R.h"
#include "ROK.h"
#include "RQED.h"
#include "RCounterX.h"
#include "RCounterY.h"
#include "Common/Color.hpp"
#include "FullyDiff/KinematicVars.hpp"

#define build(proj) cdbl FullyDiff::ME::R_##proj(cdbl m2, cdbl q2, cdbl sp, cdbl x, cdbl y, cdbl Theta1, cdbl Theta2) {\
    if (1. - x < 1e-7) return RCounterX_##proj(m2,q2,sp,y,Theta1,Theta2);\
    if (1. + y < 1e-7) return RCounterY_##proj(m2,q2,sp,x,Theta1,Theta2);\
    const KinematicVars vs(m2,q2,sp,x,y,Theta1,Theta2);\
    return Color::CA * ROK_##proj(m2,q2,sp,vs.t1,vs.u1,vs.tp,vs.up) + 2.*Color::CF * RQED_##proj(m2,q2,sp,vs.t1,vs.u1,vs.tp,vs.up);\
}

interateAllProj(build)