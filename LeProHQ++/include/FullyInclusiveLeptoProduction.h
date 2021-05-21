#ifndef FullyInclusiveLeptoProduction_H
#define FullyInclusiveLeptoProduction_H

#include "Common/AbstractFixedOrderLeptoProduction.h"

#include "FullyInclusive/IntKer.h"

/**
 * @class FullyInclusiveLeptoProduction
 * @brief inclusive lepto production of heavy quarks
 */
class FullyInclusiveLeptoProduction : public Common::AbstractFixedOrderLeptoProduction {
    
/** @brief cast ker to my type */
    #define FIker ((FullyInclusive::IntKer*)this->ker)
    
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
    FullyInclusiveLeptoProduction(cuint nlf, cdbl m2);
    
/** @brief destructor */
    ~FullyInclusiveLeptoProduction();
    
/** @name global getter and setter */
///@{
    
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
    
/** @brief threshold limit of cg0_VV() */
    cdbl cg0t_VV() const;
/** @brief threshold limit of cg0_VA() */
    cdbl cg0t_VA() const;
/** @brief threshold limit of cg0_AA() */
    cdbl cg0t_AA() const;
    
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
    
/** @brief threshold limit of cq1_VV() */
    cdbl cq1t_VV() const;
/** @brief threshold limit of cq1_VA() */
    cdbl cq1t_VA() const;
/** @brief threshold limit of cq1_AA() */
    cdbl cq1t_AA() const;
    
/** @brief improved threshold limit of cq1_VV() */
    cdbl cq1tp_VV() const;
/** @brief improved threshold limit of cq1_VA() */
    cdbl cq1tp_VA() const;
/** @brief improved threshold limit of cq1_AA() */
    cdbl cq1tp_AA() const;
    
/** @brief bulk of cq1_VV() */
    cdbl cq1b_VV() const;
/** @brief bulk of cq1_VA() */
    cdbl cq1b_VA() const;
/** @brief bulk of cq1_AA() */
    cdbl cq1b_AA() const;
    
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
    
/** @brief threshold limit of cg1_VV() */
    cdbl cg1t_VV() const;
/** @brief threshold limit of cg1_VA() */
    cdbl cg1t_VA() const;
/** @brief threshold limit of cg1_AA() */
    cdbl cg1t_AA() const;
    
/** @brief improved threshold limit of cg1_VV() */
    cdbl cg1tp_VV() const;
/** @brief improved threshold limit of cg1_VA() */
    cdbl cg1tp_VA() const;
/** @brief improved threshold limit of cg1_AA() */
    cdbl cg1tp_AA() const;
    
/** @brief bulk of cg1_VV() */
    cdbl cg1b_VV() const;
/** @brief bulk of cg1_VA() */
    cdbl cg1b_VA() const;
/** @brief bulk of cg1_AA() */
    cdbl cg1b_AA() const;
    
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
    
///@}
    
    
/** @name leptonic cross sections */
///@{

/** @see Common::AbstractLeptoProduction::sigma() */
    cdbl sigma() const override;

///@}
};

#endif // FullyInclusiveLeptoProduction_H
