#ifndef FullyDiff_ME_A3_H
#define FullyDiff_ME_A3_H

#include "config.h"

namespace FullyDiff {
    
namespace ME {
    
#define A3(proj) cdbl A3_##proj(cdbl m2, cdbl q2, cdbl sp, cdbl t1, cdbl u1, cdbl tp, cdbl up);

interateAllProj(A3)
    
} // namespace ME

} // namespace FullyDiff

#endif // FullyDiff_ME_A3_H