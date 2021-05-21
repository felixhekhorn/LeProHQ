#ifndef Inclusive_ME_IntA3_H
#define Inclusive_ME_IntA3_H

#include "config.h"

namespace Inclusive {

namespace ME {
    
#define IntA3(proj) cdbl IntA3_##proj(cdbl m2, cdbl q2, cdbl sp, cdbl t1, cdbl s4);

interateAllProj(IntA3)

} // namespace ME

} // namespace Inclusive

#endif // Inclusive_ME_IntA3_H
