#ifndef FullyInclusive_PartonicVars_HPP
#define FullyInclusive_PartonicVars_HPP

#include "../config.h"

namespace FullyInclusive {
    
typedef cdbl (*fPtr2dbl)(cdbl m2, cdbl q2);
typedef cdbl (*fPtr3dbl)(cdbl m2, cdbl q2, cdbl sp);
    
#define partonic_rho(m2,s) cdbl rho = 4.0 * m2 / s;
#define partonic_beta(m2,s) partonic_rho(m2,s) cdbl beta = sqrt(1.0 - rho);
#define partonic_chi(m2,s) partonic_beta(m2,s) cdbl chi = (1.0 - beta) / (1.0 + beta);
#define partonic_eta(m2,s) partonic_rho(m2,s) cdbl eta = 1./rho - 1.;

#define partonic_rhoq(m2,q2) cdbl rhoq = 4.0 * m2 / q2;
#define partonic_betaq(m2,q2) partonic_rhoq(m2,q2) cdbl betaq = sqrt(1.0 - rhoq);
#define partonic_chiq(m2,q2) partonic_betaq(m2,q2) cdbl chiq = (betaq - 1.0) / (betaq + 1.0);

#define partonic_rhop(m2,q2,s) cdbl rhop = 4.0 * m2 / (s - q2);
#define partonic_betap(m2,q2,s) partonic_rhop(m2,q2,s) cdbl betap = sqrt(1.0 - rhop);
#define partonic_chip(m2,q2,s) partonic_betap(m2,q2,s) cdbl chip = (1.0 - betap) / (1.0 + betap);

} // namespace FullyInclusive

#endif // FullyInclusive_PartonicVars_HPP