#ifndef Common_ME_BQED_H
#define Common_ME_BQED_H

#include "config.h"

namespace Common {

namespace ME {

cdbl BQED_F2_VV(cdbl m2, cdbl q2, cdbl sp, cdbl t1);
cdbl BQED_F2_AA(cdbl m2, cdbl q2, cdbl sp, cdbl t1);
cdbl BQED_FL_VV(cdbl m2, cdbl q2, cdbl sp, cdbl t1);
cdbl BQED_FL_AA(cdbl m2, cdbl q2, cdbl sp, cdbl t1);
cdbl BQED_x2g1_VV(cdbl m2, cdbl q2, cdbl sp, cdbl t1);
cdbl BQED_x2g1_AA(cdbl m2, cdbl q2, cdbl sp, cdbl t1);

cdbl BQED_xF3_VA(cdbl m2, cdbl q2, cdbl sp, cdbl t1);
cdbl BQED_g4_VA(cdbl m2, cdbl q2, cdbl sp, cdbl t1);
cdbl BQED_gL_VA(cdbl m2, cdbl q2, cdbl sp, cdbl t1);

} // namespace ME

} // namespace Common

#endif // Common_ME_BQED_H