#ifndef FullyInclusiveSoftResummed_IntKer_H
#define FullyInclusiveSoftResummed_IntKer_H

#include "../Common/AbstractIntKer.h"
#include "Interpolation.h"

/**
 * @brief namespace for fully inclusive soft resummed lepto production of heavy quarks
 */
namespace FullyInclusiveSoftResummed {

/**
 * @brief Fully inclusive soft resummed integration kernel
 */
class IntKer : public Common::AbstractIntKer {
/** @name integration variables */
///@{
    
/** @brief jacobian for the path */
    dcmplx jac_path;

///@}
    
/**
 * @brief set path and its jacobian
 * @param a integration variable
 */
    void setPath(cdbl a);
    
/**
 * @brief (full) LO coefficient function
 * @return \f$c_{\Pg}^{(0),thresh}\f$
 */
    cdcmplx cg0t() const;
    
/** 
 * @brief partonic wrapper
 * @param a1 integration variable
 */
    cdbl runPartonic(cdbl a1);
    
/** 
 * @brief hadronic wrapper
 * @param a1 integration variable
 */
    cdbl runHadronic(cdbl a1);
    
public:

/** @brief constructor */
    IntKer();
    
/** @brief destructor */
    ~IntKer();
    
/** @brief interpolator to bridge between x and N space */
    Interpolation::Interpolator* interpolator = 0;
    
/** @brief path */
    dcmplx N;
    
/** @name path variables */
///@{
    
/** @brief inversion path offset on real axis */
    dbl path_offset = 2.;
    
/** @brief inversion path angle towards the imaginary axis */
    dbl path_angle = M_PI_4;
    
/** @brief inversion path length */
    dbl path_length = 30.;

///@}
    
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
///@}

};

}

#endif // FullyInclusiveSoftResummed_IntKer_H