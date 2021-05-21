#ifndef InclusiveLeptoProduction_H
#define InclusiveLeptoProduction_H

#include "Common/AbstractFixedOrderLeptoProduction.h"

#include "Inclusive/IntKer.h"

/**
 * @class InclusiveLeptoProduction
 * @brief inclusive lepto production of heavy quarks
 */
class InclusiveLeptoProduction : public Common::AbstractFixedOrderLeptoProduction {
    
/** @brief cast ker to my type */
    #define kker ((Inclusive::IntKer*)this->ker)
    
/**
 * @brief runs the integration
 * @return dF_dHAQX
 */
    cdbl rundF_dHAQX() const;
    
/**
 * @brief wrapper to 1D integration
 * @return integral
 */
    cdbl int1D() const;
    
/**
 * @brief wrapper to 2D integration
 * @return integral
 */
    cdbl int2D() const;
    
/**
 * @brief wrapper to 3D integration
 * @return integral
 */
    cdbl int3D() const;
    
/**
 * @brief wrapper to 4D integration
 * @return integral
 */
    cdbl int4D() const;
    
/**
 * @brief wrapper to 5D integration
 * @return integral
 */
    cdbl int5D() const;
    
public:

/**
 * @brief constructor
 * @param nlf number of light flavours
 * @param m2 mass^2 of heavy flavours
 * @param Delta soft cut-off \f$\Delta\f$
 */
    InclusiveLeptoProduction(cuint nlf, cdbl m2, cdbl Delta);
    
/** @brief destructor */
    ~InclusiveLeptoProduction();
    
/** @name global getter and setter */
///@{

/**
 * @brief sets soft cut-off \f$\Delta\f$
 * @param Delta soft cut-off \f$\Delta\f$
 */
    inline void setDelta(cdbl Delta) {
        if (Delta < 0.)
            throw domain_error("Delta has to be positive!");
        kker->Delta = Delta;
    };
    
/** @see AbstractLeptoProduction::getIntegrationConfig */
    Common::IntegrationConfig* getIntegrationConfig(const str& method) const override;
        
///@}

/** @name partonic functions */
///@{
    
/** @see AbstractLeptoProduction::cg0_VV() */
    cdbl cg0_VV() const override;
/** @see AbstractLeptoProduction::cg0_VA() */
    cdbl cg0_VA() const override;
/** @see AbstractLeptoProduction::cg0_AA() */
    cdbl cg0_AA() const override;
    
/** @see AbstractLeptoProduction::dq1_VV() */
    cdbl dq1_VV() const override;
/** @see AbstractLeptoProduction::dq1_VA() */
    cdbl dq1_VA() const override;
/** @see AbstractLeptoProduction::dq1_AA() */
    cdbl dq1_AA() const override;
    
/** @see AbstractLeptoProduction::cq1_VV() */
    cdbl cq1_VV() const override;
/** @see AbstractLeptoProduction::cq1_VA() */
    cdbl cq1_VA() const override;
/** @see AbstractLeptoProduction::cq1_AA() */
    cdbl cq1_AA() const override;
    
/** @see AbstractLeptoProduction::cqBarF1_VV() */
    cdbl cqBarF1_VV() const override;
/** @see AbstractLeptoProduction::cqBarF1_VA() */
    cdbl cqBarF1_VA() const override;
/** @see AbstractLeptoProduction::cqBarF1_AA() */
    cdbl cqBarF1_AA() const override;
    
/** @see AbstractLeptoProduction::cg1_VV() */
    cdbl cg1_VV() const override;
/** @see AbstractLeptoProduction::cg1_VA() */
    cdbl cg1_VA() const override;
/** @see AbstractLeptoProduction::cg1_AA() */
    cdbl cg1_AA() const override;
    
/** @see AbstractLeptoProduction::cgBarF1_VV() */
    cdbl cgBarF1_VV() const override;
/** @see AbstractLeptoProduction::cgBarF1_VA() */
    cdbl cgBarF1_VA() const override;
/** @see AbstractLeptoProduction::cgBarF1_AA() */
    cdbl cgBarF1_AA() const override;
    
/** @see AbstractLeptoProduction::cgBarR1_VV() */
    cdbl cgBarR1_VV() const override;
/** @see AbstractLeptoProduction::cgBarR1_VA() */
    cdbl cgBarR1_VA() const override;
/** @see AbstractLeptoProduction::cgBarR1_AA() */
    cdbl cgBarR1_AA() const override;
    
/** @see AbstractLeptoProduction::cgBar1_VV() */
    cdbl cgBar1_VV() const override;
/** @see AbstractLeptoProduction::cgBar1_VA() */
    cdbl cgBar1_VA() const override;
/** @see AbstractLeptoProduction::cgBar1_AA() */
    cdbl cgBar1_AA() const override;
    
///@}

/** @name hadronic functions */
///@{

/** @see Common::AbstractLeptoProduction::F() */
    cdbl F() const override;
    
/**
 * @brief derivative of structure function to transverse momentum of heavy anti quark
 * @param HAQTransverseMomentum transverse momentum  of heavy anti quark
 * @return dF_dHAQTransverseMomentum(pt)
 */
    cdbl dF_dHAQTransverseMomentum(cdbl HAQTransverseMomentum) const;
    
/**
 * @brief derivative of structure function to transverse momentum scaling of heavy anti quark
 * @param HAQTransverseMomentumScaling transverse momentum scaling of heavy anti quark
 * @see dF_dHAQTransverseMomentum
 * @return dF_dHAQTransverseMomentumScaling(xt)
 */
    cdbl dF_dHAQTransverseMomentumScaling(cdbl HAQTransverseMomentumScaling) const;
    
/**
 * @brief derivative of structure function to rapidity of heavy anti quark
 * @param HAQRapidity rapidity of heavy anti quark
 * @return dF_dHAQRapidity(y)
 */
    cdbl dF_dHAQRapidity(cdbl HAQRapidity) const;
    
/**
 * @brief derivative of structure function to Feynman-x of heavy anti quark
 * @param HAQFeynmanX Feynman-x of heavy anti quark
 * @return dF_dHAQFeynmanX(xF)
 */
    cdbl dF_dHAQFeynmanX(cdbl HAQFeynmanX) const;
    
///@}
    
    
/** @name leptonic cross sections */
///@{

/** @see Common::AbstractLeptoProduction::sigma() */
    cdbl sigma() const override;

///@}
};

#endif // InclusiveLeptoProduction_H
