#ifndef FullyDiff_ME_A1Counter_H
#define FullyDiff_ME_A1Counter_H

#include "config.h"

namespace FullyDiff {
    
namespace ME {
    
#define A1Counter(proj) cdbl A1Counter_##proj(cdbl m2, cdbl q2, cdbl sp, cdbl x, cdbl Theta1, cdbl Theta2);

interateAllProj(A1Counter)
    
} // namespace ME

} // namespace FullyDiff

#endif // FullyDiff_ME_A1Counter_H