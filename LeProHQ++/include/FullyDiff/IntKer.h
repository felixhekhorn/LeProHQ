#ifndef FullyDiff_IntKer_H
#define FullyDiff_IntKer_H

#include "../Common/AbstractIntKer.h"
#include "PhasespacePoint.h"
#include "PhasespaceValues.hpp"
#include "histT.hpp"

/**
 * @brief namespace for fully differential lepto production of heavy quarks
 */
namespace FullyDiff {

/**
 * @brief fully differential integration kernel
 */
class IntKer : public Common::AbstractIntKer {
/** @name integration variables */
///@{
    
/** @brief Theta1 */
    dbl Theta1 = dblNaN;
/** @brief jacobian=volume for Theta1 trafo */
    cdbl V_Theta1 = M_PI;
    
/** @brief soft regulation parameter \f$\tilde\rho\f$ */
    dbl rhoTilde = dblNaN;
/** @brief current x-for-event */
    dbl xE = dblNaN;
/** @brief volume of x-for-event space */
    dbl V_xE = 0.;
/** @brief current x-for-counter-event */
    dbl xC = dblNaN;
/** @brief volume of x-for-counter-event space */
    dbl V_xC = 0.;
    
/** @brief current y-for-event */
    dbl yE = dblNaN;
/** @brief volume of y-for-event space */
    dbl V_yE = 0.;
/** @brief current y-for-counter-event */
    dbl yC = dblNaN;
/** @brief volume of y-for-counter-event space */
    dbl V_yC = 0.;
    
/** @brief Theta2 */
    dbl Theta2 = dblNaN;
/** @brief volume for Theta2 trafo */
    cdbl V_Theta2 = M_PI;
    
/** @brief z */
    dbl z = dblNaN;
/** @brief volume for z trafo */
    dbl V_z = 0.;
    
/**
 * @brief sets Theta1
 * @param a integration variable
 */
    void setTheta1(cdbl a);
    
/**
 * @brief sets x
 * @param a integration variable
 */
    void setX(cdbl a);
    
/**
 * @brief sets Theta2
 * @param a integration variable
 */
    void setTheta2(cdbl a);
    
/**
 * @brief sets y
 * @param a integration variable
 */
    void setY(cdbl a);
    
/**
 * @brief sets z
 * @param a integration variable
 */
    void setZ(cdbl a);
    
///@}

/** @name hadronic variables */
///@{
    
/**
 * @brief integration weight determined by VEGAS
 */
    cdbl* vegasWeight = 0;
    
///@}
    
/** @brief define 5-point-function shortcut */
    typedef cdbl (*fPtr5dbl)(cdbl m2, cdbl q2, cdbl sp, cdbl t1, cdbl betaTilde);
/** @brief define 6-point-function shortcut */
    typedef cdbl (*fPtr6dbl)(cdbl m2, cdbl q2, cdbl sp, cdbl x_or_y, cdbl Theta1, cdbl Theta2);
/** @brief define 7-point-function shortcut */
    typedef cdbl (*fPtr7dbl)(cdbl m2, cdbl q2, cdbl sp, cdbl t1, cdbl u1, cdbl tp, cdbl up);
    
/**
 * @brief sets correct pointers to A1
 * @param fVV vector-vector part
 * @param fVA vector-axial part
 * @param fAA axial-axial part
 */
    void getA1(fPtr7dbl &fVV, fPtr7dbl &fVA, fPtr7dbl &fAA) const;
    
/**
 * @brief sets correct pointers to A1Counter
 * @param fVV vector-vector part
 * @param fVA vector-axial part
 * @param fAA axial-axial part
 */
    void getA1Counter(fPtr6dbl &fVV, fPtr6dbl &fVA, fPtr6dbl &fAA) const;
    
/**
 * @brief sets correct pointers to A2
 * @param fVV vector-vector part
 * @param fVA vector-axial part
 * @param fAA axial-axial part
 */
    void getA2(fPtr7dbl &fVV, fPtr7dbl &fVA, fPtr7dbl &fAA) const;
    
/**
 * @brief sets correct pointers to A3
 * @param fVV vector-vector part
 * @param fVA vector-axial part
 * @param fAA axial-axial part
 */
    void getA3(fPtr7dbl &fVV, fPtr7dbl &fVA, fPtr7dbl &fAA) const;
    
/**
 * @brief sets correct pointers to R
 * @param fVV vector-vector part
 * @param fVA vector-axial part
 * @param fAA axial-axial part
 */
    void getR(fPtr7dbl &fVV, fPtr7dbl &fVA, fPtr7dbl &fAA) const;
    
/**
 * @brief sets correct pointers to RCounterX
 * @param fVV vector-vector part
 * @param fVA vector-axial part
 * @param fAA axial-axial part
 */
    void getRCounterX(fPtr6dbl &fVV, fPtr6dbl &fVA, fPtr6dbl &fAA) const;
    
/**
 * @brief sets correct pointers to RCounterY
 * @param fVV vector-vector part
 * @param fVA vector-axial part
 * @param fAA axial-axial part
 */
    void getRCounterY(fPtr6dbl &fVV, fPtr6dbl &fVA, fPtr6dbl &fAA) const;
    
/**
 * @brief sets correct pointers to RCounterXY
 * @param fVV vector-vector part
 * @param fVA vector-axial part
 * @param fAA axial-axial part
 */
    void getRCounterXY(fPtr4dbl &fVV, fPtr4dbl &fVA, fPtr4dbl &fAA) const;
    
/**
 * @brief sets correct pointers to SV
 * @param fVV vector-vector part
 * @param fVA vector-axial part
 * @param fAA axial-axial part
 */
    void getSV(fPtr5dbl &fVV, fPtr5dbl &fVA, fPtr5dbl &fAA) const;

/** @name partonic coefficient functions */
///@{
    
/** 
 * @brief computes cg0
 * @return cg0
 */
    cdbl cg0() const;
    
/** 
 * @brief computes cg1
 * @return cg1
 */
    const FullyDiff::PhasespaceValues cg1() const;
    
/** 
 * @brief computes cgBarF1
 * @return cgBarF1
 */
    const FullyDiff::PhasespaceValues cgBarF1() const;
    
/** 
 * @brief computes cgBarR1
 * @return cgBarR1
 */
    const FullyDiff::PhasespaceValues cgBarR1() const;
    
/** 
 * @brief computes cq1
 * @return cq1
 */
    const FullyDiff::PhasespaceValues cq1() const;
    
/** 
 * @brief computes cqBarF1
 * @return cqBarF1
 */
    const FullyDiff::PhasespaceValues cqBarF1() const;
    
/** 
 * @brief computes dq1
 * @return dq1
 */
    const FullyDiff::PhasespaceValues dq1() const;
    
/** 
 * @brief computes oq1
 * @return oq1
 */
    const FullyDiff::PhasespaceValues oq1() const;
    
///@}
    
/**
 * @brief computes a single partonic function
 * @param a1 integration variable
 * @param a2 integration variable
 * @param a3 integration variable
 * @param a4 integration variable
 * @return a single partonic function
 */
    cdbl runPartonic(cdbl a1, cdbl a2, cdbl a3, cdbl a4);
    
/**
 * @brief computes a hadronic function
 * @param a1 integration variable
 * @param a2 integration variable
 * @param a3 integration variable
 * @param a4 integration variable
 * @param a5 integration variable
 * @return a single hadronic function
 */
    cdbl runHadronic(cdbl a1, cdbl a2, cdbl a3, cdbl a4, cdbl a5);
    
/**
 * @brief computes a hadronic function
 * @param a1 integration variable
 * @param a2 integration variable
 * @param a3 integration variable
 * @param a4 integration variable
 * @param a5 integration variable
 * @return a single leptonic function (currently only leptonic xs)
 */
    cdbl runLeptonic(cdbl a1, cdbl a2, cdbl a3, cdbl a4, cdbl a5);

/** @name hadronic structure functions */
///@{
    
/**
 * @brief combines all parts to hadronic NLOg
 * @param x
 * @param y
 * @param cg1
 * @param cgBarR1
 * @param cgBarF1
 * @return total kernel
 */
    cdbl combineNLOg(cdbl x, cdbl y, cdbl cg1, cdbl cgBarR1, cdbl cgBarF1);
    
/**
 * @brief combines all parts of hadronic NLOq
 * @param x
 * @param y
 * @param cq1
 * @param cqBarF1
 * @param dq1_VV
 * @param dq1_VA
 * @param dq1_AA
 * @param oq1_VV
 * @param oq1_VA
 * @param oq1_AA
 * @return total kernel
 */
    cdbl combineNLOq(cdbl x, cdbl y, cdbl cq1, cdbl cqBarF1, cdbl dq1_VV, cdbl dq1_VA, cdbl dq1_AA, cdbl oq1_VV, cdbl oq1_VA, cdbl oq1_AA);
    
/**
 * @brief scales all avtive histograms
 * @param s factor
 */
    void scaleHistograms(cdbl s) const;

/**
 * @brief fills all avtive histograms available in all orders
 * These histograms may only depend on LO variables such as z, p1, p2, k1, q
 * @param p current phase space point
 * @param i integrand
 */
    void fillAllOrderHistograms(const PhasespacePoint& p, cdbl i) const;

/**
 * @brief fills all avtive histograms available in NLO
 * @param p current phase space point
 * @param i integrand
 */
    void fillNLOHistograms(const PhasespacePoint& p, cdbl i) const;
    
///@}

/** @name additional kernel modes */
///@{
    static cuint Mode_oq1_VV = 25;
    static cuint Mode_oq1_VA = 26;
    static cuint Mode_oq1_AA = 27;
///@}
    
public:

/** @brief factor to soft regulation parameter \f$\tilde\rho = 1-\tilde x(1-\rho^*) \f$ */
    dbl xTilde = dblNaN;
    
/** @brief collinear regulation parameter \f$\omega\f$ */
    dbl omega = dblNaN;
    
/** @brief offset to upper integration bound in x \f$\delta_x \f$ */
    dbl deltax = 0.;
    
/** @brief offset to lower integration bound in y \f$\delta_y \f$ */
    dbl deltay = 0.;
    
/** @brief pointer to map of histograms */
    histMapT histMap;

/** @brief constructor */
    IntKer();
    
/** @brief destructor */
    ~IntKer();
    
/**
 * @brief called function in gsl kernel
 * @param a1 integration variable
 * @param a2 integration variable
 * @param a3 integration variable
 * @param a4 integration variable
 * @param a5 integration variable
 * @return kernel
 */
    cdbl operator()(cdbl a1, cdbl a2 = 0., cdbl a3 = 0., cdbl a4 = 0., cdbl a5 = 0.);
    
/**
 * @brief called function in Dvegas kernel
 * @param x adapted continuous integration variables
 * @param k discrete integration variables
 * @param weight integration weight
 * @param aux unadapted continuous integration variables
 * @param f output
 */
    void operator()(const double x[], const int k[], const double& weight, const double aux[], double f[]) override;
    
/**
 * @see AbstractIntKer::Dvegas_init
 */
    void Dvegas_init() const override;
    
/**
 * @see AbstractIntKer::Dvegas_final
 */
    void Dvegas_final(cuint iterations) const override;
};

}

#endif // FullyDiff_IntKer_H