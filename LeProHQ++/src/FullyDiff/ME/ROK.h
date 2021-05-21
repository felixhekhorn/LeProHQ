#ifndef FullyDiff_ME_ROK_H
#define FullyDiff_ME_ROK_H

#include "config.h"

namespace FullyDiff {
    
namespace ME {
    
#define ROK(proj) cdbl ROK_##proj(cdbl m2, cdbl q2, cdbl sp, cdbl t1, cdbl u1, cdbl tp, cdbl up);

interateAllProj(ROK)
    
} // namespace ME

} // namespace FullyDiff

#endif // FullyDiff_ME_ROK_H