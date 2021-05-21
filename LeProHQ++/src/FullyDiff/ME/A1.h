#ifndef FullyDiff_ME_A1_H
#define FullyDiff_ME_A1_H

#include "config.h"

namespace FullyDiff {
    
namespace ME {
    
#define A1(proj) cdbl A1_##proj(cdbl m2, cdbl q2, cdbl sp, cdbl t1, cdbl u1, cdbl tp, cdbl up);

interateAllProj(A1)
    
} // namespace ME

} // namespace FullyDiff

#endif // FullyDiff_ME_A1_H