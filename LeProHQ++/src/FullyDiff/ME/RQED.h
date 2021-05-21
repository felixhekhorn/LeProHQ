#ifndef FullyDiff_ME_RQED_H
#define FullyDiff_ME_RQED_H

#include "config.h"

namespace FullyDiff {
    
namespace ME {
    
#define RQED(proj) cdbl RQED_##proj(cdbl m2, cdbl q2, cdbl sp, cdbl t1, cdbl u1, cdbl tp, cdbl up);

interateAllProj(RQED)
    
} // namespace ME

} // namespace FullyDiff

#endif // FullyDiff_ME_RQED_H