#ifndef Inclusive_ME_IntA2_H
#define Inclusive_ME_IntA2_H

#include "config.h"

namespace Inclusive {

namespace ME {
    
#define IntA2(proj) cdbl IntA2_##proj(cdbl m2, cdbl q2, cdbl sp, cdbl t1, cdbl s4);

interateAllProj(IntA2)

} // namespace ME

} // namespace Inclusive

#endif // Inclusive_ME_IntA2_H