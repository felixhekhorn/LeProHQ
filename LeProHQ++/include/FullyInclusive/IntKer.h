#ifndef FullyInclusive_IntKer_H
#define FullyInclusive_IntKer_H

#include "../Common/AbstractIntKer.h"

/**
 * @brief namespace for fully inclusive lepto production of heavy quarks
 */
namespace FullyInclusive {

/**
 * @brief Fully inclusive integration kernel
 */
class IntKer : public Common::AbstractIntKer {
/** @name integration variables */
///@{
    
/** @brief parton moment variable xi */
    dbl xi = dblNaN;
    
/** @brief jacobian for z-integration */
    dbl jac_z = dblNaN;
    
/** @brief sets xi, partonic s */
    void setPartonicVars(cdbl a);
    
///@}
    
/** @brief define shortcut */
    typedef cdbl (*fPtr3dbl)(cdbl m2, cdbl q2, cdbl sp);

/** @name partonic coefficient functions */
///@{
    
/** 
 * @brief computes full cg0
 * @return cg0
 */
    cdbl cg0() const;
    
/** 
 * @brief computes full cg0
 * @return cg0
 */
    cdbl cg0t() const;

/** 
 * @brief computes full cq1
 * @return cq1
 */
    cdbl cq1() const;

/** 
 * @brief computes full cqBarF1
 * @return cqBarF1
 */
    cdbl cqBarF1() const;

/** 
 * @brief computes cgBar1
 * @return cgBar1
 */
    cdbl cgBar1() const;

/** 
 * @brief computes cgBarF1
 * @return cgBarF1
 */
    cdbl cgBarF1() const;

/** 
 * @brief computes cgBarR1
 * @return cgBarR1
 */
    cdbl cgBarR1() const;

/** 
 * @brief computes cg1
 * @return cg1
 */
    cdbl cg1() const;

///@}
    
    /** @brief partonic wrapper */
    cdbl runPartonic();
    
    /** @brief hadronic wrapper */
    cdbl runHadronic(cdbl a1);

/** @name hadronic structure functions */
///@{
    
/** 
 * @brief computes Fg0
 * @return Fg0
 */
    cdbl Fg0() const;
    
/** 
 * @brief computes Fg1
 * @return Fg1
 */
    cdbl Fg1() const;
    
/** 
 * @brief computes Fq1
 * @return Fq1
 */
    cdbl Fq1() const;
    
///@}
    
public:

/** @brief constructor */
    IntKer();
    
/** @brief destructor */
    ~IntKer();
    
/**
 * @brief called function in gsl kernel
 * @param a1 integration variable
 * @param a2 integration variable
 * @param a3 integration variable
 * @return kernel
 */
    cdbl operator()(cdbl a1 = 0., cdbl a2 = 0., cdbl a3 = 0.);
    
/**
 * @brief called function in Dvegas kernel
 * @param x adapted continuous integration variables
 * @param k discrete integration variables
 * @param weight integration weight
 * @param aux unadapted continuous integration variables
 * @param f output
 */
    void operator()(const double x[], const int k[], const double& weight, const double aux[], double f[]);

/** @name additional kernel modes */
///@{
    static cuint Mode_cg0t_VV = 41;
    static cuint Mode_cg0t_VA = 42;
    static cuint Mode_cg0t_AA = 43;
    static cuint Mode_cq1t_VV = 44;
    static cuint Mode_cq1t_VA = 45;
    static cuint Mode_cq1t_AA = 46;
    static cuint Mode_cq1tp_VV = 47;
    static cuint Mode_cq1tp_VA = 48;
    static cuint Mode_cq1tp_AA = 49;
    static cuint Mode_cq1b_VV = 50;
    static cuint Mode_cq1b_VA = 51;
    static cuint Mode_cq1b_AA = 52;
    static cuint Mode_cg1t_VV = 53;
    static cuint Mode_cg1t_VA = 54;
    static cuint Mode_cg1t_AA = 55;
    static cuint Mode_cg1tp_VV = 56;
    static cuint Mode_cg1tp_VA = 57;
    static cuint Mode_cg1tp_AA = 58;
    static cuint Mode_cg1b_VV = 59;
    static cuint Mode_cg1b_VA = 60;
    static cuint Mode_cg1b_AA = 61;
///@}
};

}

#endif // FullyInclusive_IntKer_H