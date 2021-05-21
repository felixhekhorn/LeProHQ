#ifndef FullyDiff_KinematicVars_HPP_
#define FullyDiff_KinematicVars_HPP_

#include "config.h"

/** @brief forces x to be min <= x <= max  */
#define mima(x,mi,ma) min(max(x,mi),ma)

namespace FullyDiff {

/**
 * @class KinematicVars
 * @brief computes all 2-to-3-particle phasespace variables
 */
struct KinematicVars {
    
    /** @name momentum variables computed in the cms of the heavy quark pair and q || z */
    ///@{
    dbl q0,k10,absq,cosPsi,sinPsi,beta5;
    ///@}
    
    /** @name mandelstam variables (Lorentz invariants) */
    ///@{
    dbl s5,tp,up,t1,u1;
    ///@}
    
    /**
     * @brief constructor
     * @param m2
     * @param q2
     * @param sp
     * @param x
     * @param y
     * @param Theta1
     * @param Theta2
     */
    KinematicVars(cdbl m2, cdbl q2, cdbl sp, cdbl x, cdbl y, cdbl Theta1, cdbl Theta2) {
        if (!isfinite(m2) || m2 <= 0.) throw domain_error("m2 has to be strictly positive!");
        if (!isfinite(q2) || q2 >= 0.) throw domain_error("q2 has to be strictly negative!");
        cdbl s = sp+q2;
        this->s5 = max(q2 + sp*x,4.*m2);
        this->tp = -1./2.*sp*(sp/s)*(1.-x)*(1.+y);
        this->up = q2*(s5/s)-1./2.*sp*(sp/s)*(1.-x)*(1.-y);
        
        cdbl Sqrts5 = sqrt(s5);
        this->q0 = (s+up)/(2.*Sqrts5);
        this->k10 = (s5-up)/(2.*Sqrts5);
        this->absq = 1./(2.*Sqrts5)*sqrt(pow(s+up,2)-4.*s5*q2);
        this->cosPsi = mima((-sp+2.*k10*q0)/(2.*k10*absq),-1.,1.);
        this->sinPsi = sqrt(1.-cosPsi*cosPsi);
        this->beta5 = sqrt(1.-4.*m2/s5);
        
        this->t1 = -1./2.*(s5-up)*(1. + beta5*cos(Theta2)*sin(Theta1)*sinPsi + beta5*cos(Theta1)*cosPsi);
        this->u1 = q2 - 1./2.*(s + up + 2.*Sqrts5*beta5*absq*cos(Theta1));
    }
};


} // namespace FullyDiff


#endif // FullyDiff_KinematicVars_HPP_