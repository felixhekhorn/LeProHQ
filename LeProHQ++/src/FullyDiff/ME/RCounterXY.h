#ifndef FullyDiff_ME_RCounterXY_H
#define FullyDiff_ME_RCounterXY_H

#include "config.h"

namespace FullyDiff {
    
namespace ME {
    
#define RCounterXY(proj) cdbl RCounterXY_##proj(cdbl m2, cdbl q2, cdbl sp, cdbl t1);

interateAllProj(RCounterXY)
    
} // namespace ME

} // namespace FullyDiff

#endif // FullyDiff_ME_RCounterXY_H