#ifndef FullyDiff_ME_RCounterY_H
#define FullyDiff_ME_RCounterY_H

#include "config.h"

namespace FullyDiff {
    
namespace ME {
    
#define RCounterY(proj) cdbl RCounterY_##proj(cdbl m2, cdbl q2, cdbl sp, cdbl x, cdbl Theta1, cdbl Theta2);

interateAllProj(RCounterY)
    
} // namespace ME

} // namespace FullyDiff

#endif // FullyDiff_ME_RCounterY_H