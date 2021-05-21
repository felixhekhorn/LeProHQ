#ifndef FullyDiff_ME_SV_H
#define FullyDiff_ME_SV_H

#include "config.h"

namespace FullyDiff {
    
namespace ME {

#define SV(proj)       cdbl SV_##proj(cdbl m2, cdbl q2, cdbl sp, cdbl t1, cdbl betaTilde);
#define SVOK(proj)   cdbl SVOK_##proj(cdbl m2, cdbl q2, cdbl sp, cdbl t1, cdbl betatil);
#define SVQED(proj) cdbl SVQED_##proj(cdbl m2, cdbl q2, cdbl sp, cdbl t1, cdbl betatil);
#define SVtot(proj) SV(proj) SVOK(proj) SVQED(proj)

interateAllProj(SVtot)

} // namespace ME

} // namespace FullyDiff

#endif // FullyDiff_ME_SV_H