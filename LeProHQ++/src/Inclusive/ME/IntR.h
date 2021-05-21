#ifndef Inclusive_ME_IntR_H
#define Inclusive_ME_IntR_H

#include "config.h"

namespace Inclusive {

namespace ME {
    
#define IntROK(proj)  cdbl IntROK_##proj(cdbl m2, cdbl q2, cdbl sp, cdbl t1, cdbl s4);
#define IntRQED(proj) cdbl IntRQED_##proj(cdbl m2, cdbl q2, cdbl sp, cdbl t1, cdbl s4);
#define IntR(proj) IntROK(proj) IntRQED(proj)

interateAllProj(IntR)

} // namespace ME

} // namespace Inclusive

#endif // Inclusive_ME_IntR_H