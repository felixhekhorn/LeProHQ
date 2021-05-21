#ifndef Inclusive_ME_SV_H
#define Inclusive_ME_SV_H

#include "config.h"

namespace Inclusive {

namespace ME {
    
#define SVOK(proj)  cdbl SVOK_##proj(cdbl m2, cdbl q2, cdbl sp, cdbl t1);
#define SVQED(proj) cdbl SVQED_##proj(cdbl m2, cdbl q2, cdbl sp, cdbl t1);
#define SV(proj) SVOK(proj) SVQED(proj)

interateAllProj(SV)

} // namespace ME

} // namespace Inclusive

#endif // Inclusive_ME_SV_H