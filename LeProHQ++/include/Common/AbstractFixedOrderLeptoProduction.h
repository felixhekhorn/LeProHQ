#ifndef Common_AbstractFixedOrderLeptoProduction_H_
#define Common_AbstractFixedOrderLeptoProduction_H_

#include "config.h"
#include "Common/AbstractLeptoProduction.h"
#include "Common/AbstractIntKer.h"
#include "Common/IntegrationMeta.hpp"

#include <vector>

namespace Common {

/**
 * @class AbstractFixedOrderLeptoProduction
 * @brief abstract base class for fixed order computations, that provide all partonic coefficient functions
 */
class AbstractFixedOrderLeptoProduction : public AbstractLeptoProduction {
protected:
    
public:

/**
 * @brief constructor
 * @param ker integration kernel
 * @param nlf number of light flavours
 * @param m2 mass^2 of heavy flavours
 */
    AbstractFixedOrderLeptoProduction(AbstractIntKer* ker, cuint nlf, cdbl m2);
    
/** @name partonic coefficient functions */
///@{

/**
 * @brief computes vector-vector part of leading order partonic gluon contribution \f$c_{\Pg}^{(0),VV}\f$
 * @return \f$c_{\Pg}^{(0),VV}\f$
 */
    virtual cdbl cg0_VV() const = 0;

/**
 * @brief computes vector-axial part of leading order partonic gluon contribution \f$c_{\Pg}^{(0),VA}\f$
 * @return \f$c_{\Pg}^{(0),VA}\f$
 */
    virtual cdbl cg0_VA() const = 0;

/**
 * @brief computes axial-axial part of leading order partonic gluon contribution \f$c_{\Pg}^{(0),AA}\f$
 * @return \f$c_{\Pg}^{(0),AA}\f$
 */
    virtual cdbl cg0_AA() const = 0;

/**
 * @brief computes vector-vector part of next-to-leading order partonic quark contribution with light quark charge \f$d_{\Pq}^{(1),VV}\f$
 * @return \f$d_{\Pq}^{(1),VV}\f$
 */
    virtual cdbl dq1_VV() const = 0;

/**
 * @brief computes vector-axial part of next-to-leading order partonic quark contribution with light quark charge \f$d_{\Pq}^{(1),VA}\f$
 * @return \f$d_{\Pq}^{(1),VA}\f$
 */
    virtual cdbl dq1_VA() const = 0;

/**
 * @brief computes axial-axial part of next-to-leading order partonic quark contribution with light quark charge \f$d_{\Pq}^{(1),AA}\f$
 * @return \f$d_{\Pq}^{(1),AA}\f$
 */
    virtual cdbl dq1_AA() const = 0;

/**
 * @brief computes vector-vector part of next-to-leading order partonic quark contribution with heavy quark charge \f$c_{\Pq}^{(1),VV}\f$
 * @return \f$c_{\Pq}^{(1),VV}\f$
 */
    virtual cdbl cq1_VV() const = 0;

/**
 * @brief computes vector-axial part of next-to-leading order partonic quark contribution with heavy quark charge \f$c_{\Pq}^{(1),VA}\f$
 * @return \f$c_{\Pq}^{(1),VA}\f$
 */
    virtual cdbl cq1_VA() const = 0;

/**
 * @brief computes axial-axial part of next-to-leading order partonic quark contribution with heavy quark charge \f$c_{\Pq}^{(1),AA}\f$
 * @return \f$c_{\Pq}^{(1),AA}\f$
 */
    virtual cdbl cq1_AA() const = 0;

/**
 * @brief computes vector-vector part of next-to-leading order partonic quark contribution with factorisation scale \f$\bar c_{\Pq}^{F,(1),VV}\f$
 * @return \f$\bar c_{\Pq}^{F,(1),VV}\f$
 */
    virtual cdbl cqBarF1_VV() const = 0;

/**
 * @brief computes vector-axial part of next-to-leading order partonic quark contribution with factorisation scale \f$\bar c_{\Pq}^{F,(1),VA}\f$
 * @return \f$\bar c_{\Pq}^{F,(1),VA}\f$
 */
    virtual cdbl cqBarF1_VA() const = 0;

/**
 * @brief computes axial-axial part of next-to-leading order partonic quark contribution with factorisation scale \f$\bar c_{\Pq}^{F,(1),AA}\f$
 * @return \f$\bar c_{\Pq}^{F,(1),AA}\f$
 */
    virtual cdbl cqBarF1_AA() const = 0;

/**
 * @brief computes vector-vector part of next-to-leading order gluon contribution \f$c_{\Pg}^{(1),VV}\f$
 * @return \f$c_{\Pg}^{(1),VV}\f$
 */
    virtual cdbl cg1_VV() const = 0;

/**
 * @brief computes vector-axial part of next-to-leading order gluon contribution \f$c_{\Pg}^{(1),VA}\f$
 * @return \f$c_{\Pg}^{(1),VA}\f$
 */
    virtual cdbl cg1_VA() const = 0;

/**
 * @brief computes axial-axial part of next-to-leading order gluon contribution \f$c_{\Pg}^{(1),AA}\f$
 * @return \f$c_{\Pg}^{(1),AA}\f$
 */
    virtual cdbl cg1_AA() const = 0;

/**
 * @brief computes vector-vector part of next-to-leading order partonic gluon contribution with factorisation scale \f$\bar c_{\Pg}^{F,(1),VV}\f$
 * @return \f$\bar c_{\Pg}^{F,(1),VV}\f$
 */
    virtual cdbl cgBarF1_VV() const = 0;

/**
 * @brief computes vector-axial part of next-to-leading order partonic gluon contribution with factorisation scale \f$\bar c_{\Pg}^{F,(1),VA}\f$
 * @return \f$\bar c_{\Pg}^{F,(1),VA}\f$
 */
    virtual cdbl cgBarF1_VA() const = 0;

/**
 * @brief computes axial-axial part of next-to-leading order partonic gluon contribution with factorisation scale \f$\bar c_{\Pg}^{F,(1),AA}\f$
 * @return \f$\bar c_{\Pg}^{F,(1),AA}\f$
 */
    virtual cdbl cgBarF1_AA() const = 0;

/**
 * @brief computes vector-vector part of next-to-leading order partonic gluon contribution with renormalization scale \f$\bar c_{\Pg}^{R,(1),VV}\f$
 * @return \f$\bar c_{\Pg}^{R,(1),VV}\f$
 */
    virtual cdbl cgBarR1_VV() const = 0;

/**
 * @brief computes vector-axial part of next-to-leading order partonic gluon contribution with renormalization scale \f$\bar c_{\Pg}^{R,(1),VA}\f$
 * @return \f$\bar c_{\Pg}^{R,(1),VA}\f$
 */
    virtual cdbl cgBarR1_VA() const = 0;

/**
 * @brief computes axial-axial part of next-to-leading order partonic gluon contribution with renormalization scale \f$\bar c_{\Pg}^{R,(1),AA}\f$
 * @return \f$\bar c_{\Pg}^{R,(1),AA}\f$
 */
    virtual cdbl cgBarR1_AA() const = 0;

/**
 * @brief computes vector-vector part of next-to-leading order partonic gluon contribution with common scale \f$\bar c_{\Pg}^{(1),VV} = \bar c_{\Pg}^{F,(1),VV} + \bar c_{\Pg}^{R,(1),VV}\f$
 * @return \f$\bar c_{\Pg}^{R,(1),VV}\f$
 */
    virtual cdbl cgBar1_VV() const = 0;

/**
 * @brief computes vector-axial part of next-to-leading order partonic gluon contribution with common scale \f$\bar c_{\Pg}^{(1),VA} = \bar c_{\Pg}^{F,(1),VA} + \bar c_{\Pg}^{R,(1),VA}\f$
 * @return \f$\bar c_{\Pg}^{R,(1),VA}\f$
 */
    virtual cdbl cgBar1_VA() const = 0;

/**
 * @brief computes axial-axial part of next-to-leading order partonic gluon contribution with common scale \f$\bar c_{\Pg}^{(1),AA} = \bar c_{\Pg}^{F,(1),AA} + \bar c_{\Pg}^{R,(1),AA}\f$
 * @return \f$\bar c_{\Pg}^{R,(1),AA}\f$
 */
    virtual cdbl cgBar1_AA() const = 0;

///@}
        
};

} // namespace Common

#endif // Common_AbstractFixedOrderLeptoProduction_H_
