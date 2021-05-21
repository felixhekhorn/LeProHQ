#ifndef FullyDiff_ME_RCounterX_H
#define FullyDiff_ME_RCounterX_H

#include "config.h"

namespace FullyDiff {
    
namespace ME {

cdbl EikonalFactorOKCounterX(cdbl m2, cdbl q2, cdbl sp, cdbl y, cdbl Theta1, cdbl Theta2);    
cdbl EikonalFactorQEDCounterX(cdbl m2, cdbl q2, cdbl sp, cdbl y, cdbl Theta1, cdbl Theta2);

#define RCounterX(proj) cdbl RCounterX_##proj(cdbl m2, cdbl q2, cdbl sp, cdbl y, cdbl Theta1, cdbl Theta2);

interateAllProj(RCounterX)
    
} // namespace ME

} // namespace FullyDiff

#endif // FullyDiff_ME_RCounterX_H