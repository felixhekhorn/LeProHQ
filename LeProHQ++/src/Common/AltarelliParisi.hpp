#ifndef AltarelliParisi_HPP
#define AltarelliParisi_HPP

#include "config.h"
#include "Common/Color.hpp"

using namespace Color;

/**
 * @brief collects all Altarelli-Parisi kernels
 * 
 * gluon splitting another gluon:
 * \f[(\Delta)P_{\Pg\Pg}(z) = (\Delta) P_{\Pg\Pg}^{H}(z) + (\Delta) P_{\Pg\Pg}^S(z)\delta(1-z)\f]
 * \f[(\Delta)P_{\Pg\Pg}^H(z) = (\Delta) P_{\Pg\Pg}^{H,(0)}(z) + \epsilon (\Delta) P_{\Pg\Pg}^{H,(1)}(z)\f]
 * \f[(\Delta)P_{\Pg\Pg}^S(z) = \frac {\beta_0}{2} + \epsilon (\Delta) P_{\Pg\Pg}^{S,(1)}\f]
 * 
 * quark splitting a gluon:
 * \f[(\Delta)P_{\Pg\Pq}(z) = (\Delta)P_{\Pg\Pq}^{(0)}(z) + \epsilon (\Delta)P_{\Pg\Pq}^{(1)}(z)\f]
 */
namespace AltarelliParisi {
    
/**
 * @brief quark-splitting-gluon AP-kernel
 * @param x momentum fraction
 * @return \f$P_{\Pg\Pq}^{(0)}(x)\f$
 */
inline cdbl Pgq0 (cdbl x) { return CF*(1.+pow(1.-x,2))/x; }

/**
 * @brief quark-splitting-gluon AP-kernel
 * @param x momentum fraction
 * @return \f$P_{\Pg\Pq}^{(1)}(x)\f$
 */
inline cdbl Pgq1 (cdbl x) { return CF*x/2.; }
    
/**
 * @brief polarized quark-splitting-gluon AP-kernel
 * @param x momentum fraction
 * @return \f$\Delta P_{\Pg\Pq}^{(0)}(x)\f$
 */
inline cdbl DeltaPgq0 (cdbl x) { return CF*(2.-x); }

/**
 * @brief polarized quark-splitting-gluon AP-kernel
 * @param x momentum fraction
 * @return \f$\Delta P_{\Pg\Pq}^{(1)}(x)\f$
 */
inline cdbl DeltaPgq1 (cdbl x) { return CF*(-(1.-x)); }

/**
 * @brief unregulated hard part of gluon-splitting-gluon AP-kernel
 * @param x momentum fraction
 * @return \f$P_{\Pg\Pg}^{H,(0)}(x)\f$
 */
inline cdbl PggH0 (cdbl x) { return 2.*CA*(x/(1.-x) + (1.-x)/x + x*(1.-x)); }

/**
 * @brief unregulated hard part of gluon-splitting-gluon AP-kernel
 * @param x momentum fraction
 * @return \f$P_{\Pg\Pg}^{H,(1)}(x)\f$
 */
inline cdbl PggH1 (cdbl x) { return 0.; }

/**
 * @brief polarized unregulated hard part of gluon-splitting-gluon AP-kernel
 * @param x momentum fraction
 * @return \f$\Delta P_{\Pg\Pg}^{H,(0)}(x)\f$
 */
inline cdbl DeltaPggH0 (cdbl x) { return 2.*CA*(1./(1.-x) - 2.*x + 1.); }

/**
 * @brief polarized unregulated hard part of gluon-splitting-gluon AP-kernel
 * @param x momentum fraction
 * @return \f$\Delta P_{\Pg\Pg}^{H,(1)}(x)\f$
 */
inline cdbl DeltaPggH1 (cdbl x) { return 2.*CA*(-(1.-x)); }

/**
 * @brief soft part of gluon-splitting-gluon AP-kernel
 * @return \f$P_{\Pg\Pg}^{S,(1)}\f$
 */
inline cdbl PggS1 () { return 0.; }

/**
 * @brief polarized soft part of gluon-splitting-gluon AP-kernel
 * @return \f$\Delta P_{\Pg\Pg}^{S,(1)}\f$
 */
inline cdbl DeltaPggS1 () { return -((cdbl)NC)/12.; }

} // namespace AltarelliParisi

#endif // AltarelliParisi_HPP