#ifndef FullyInclusiveSoftResummedLeptoProduction_H
#define FullyInclusiveSoftResummedLeptoProduction_H

#include "Common/AbstractLeptoProduction.h"

#include "FullyInclusiveSoftResummed/IntKer.h"

/**
 * @class FullyInclusiveSoftResummedLeptoProduction
 * @brief fully inclusive soft resummed lepto production of heavy quarks
 */
class FullyInclusiveSoftResummedLeptoProduction : public Common::AbstractLeptoProduction {
    
/** @brief cast ker to my type */
    #define FISRker ((FullyInclusiveSoftResummed::IntKer*)this->ker)
    
/**
 * @brief wrapper to 1D integration
 * @return integral
 */
    cdbl int1D() const;
    
/**
 * @brief wrapper to 3D integration
 * @return integral
 */
    cdbl int3D() const;
    
public:

/**
 * @brief constructor
 * @param nlf number of light flavours
 * @param m2 mass^2 of heavy flavours
 */
    FullyInclusiveSoftResummedLeptoProduction(cuint nlf, cdbl m2);
    
/** @brief destructor */
    ~FullyInclusiveSoftResummedLeptoProduction();
    
/** @name global getter and setter */
///@{
    
/** @see AbstractLeptoProduction::getIntegrationConfig */
    Common::IntegrationConfig* getIntegrationConfig(const str& method) const override;
    
/**
 * @brief Sets the inversion path parameters
 * @param offset intersection with real axis
 * @param angle angle with imaginary axis
 * @param length path distance
 */
    void setPath(cdbl offset, cdbl angle, cdbl length);

///@}

/** @name hadronic setter */
///@{

/**
 * @brief Initialize interpolation setup
 * @param xgrid grid of support points in x-space
 * @param degree degree of interpolation polynomial
 * @param isLog interpolation in log(x)?
 */
    void setInterpolation(FullyInclusiveSoftResummed::Interpolation::vdbl xgrid, cuint degree, bool isLog = true);
    
///@}

/** @name partonic functions */
///@{
    
/** @brief threshold limit of cg0_VV() */
    cdbl cg0t_VV() const;
/** @brief threshold limit of cg0_VA() */
    cdbl cg0t_VA() const;
/** @brief threshold limit of cg0_AA() */
    cdbl cg0t_AA() const;

///@}

/** @name hadronic functions */
///@{

/** @see Common::AbstractLeptoProduction::F() */
    cdbl F() const override;
    
///@}
    
    
/** @name leptonic cross sections */
///@{

/** @see Common::AbstractLeptoProduction::sigma() */
    cdbl sigma() const override;

///@}
};

#endif // FullyInclusiveSoftResummedLeptoProduction_H
