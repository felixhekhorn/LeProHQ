#ifndef FullyDiff_ME_A2_H
#define FullyDiff_ME_A2_H

#include "config.h"

namespace FullyDiff {
    
namespace ME {
    
#define A2(proj) cdbl A2_##proj(cdbl m2, cdbl q2, cdbl sp, cdbl t1, cdbl u1, cdbl tp, cdbl up);

interateAllProj(A2)
    
} // namespace ME

} // namespace FullyDiff

#endif // FullyDiff_ME_A2_H