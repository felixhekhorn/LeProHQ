#include "FullyInclusiveLeptoProduction.h"

#include <gsl/gsl_integration.h>
#include <boost/algorithm/string.hpp>

#include "gslpp/gslpp.Functor.hpp"

#include "Common/Integration.hpp"

FullyInclusiveLeptoProduction::FullyInclusiveLeptoProduction(cuint nlf, cdbl m2) :
    AbstractFixedOrderLeptoProduction(new FullyInclusive::IntKer(), nlf, m2) {
    // setup default values for IntegrationConfig
    this->intConfigs.resize(2);
    Common::IntegrationConfig* i1 = new Common::IntegrationConfig;
    i1->method = "gsl_integration_qag";
    i1->calls = 10000;
    i1->GslQag_epsabs = 5e-6;
    i1->GslQag_epsrel = 1e-4;
    i1->GslQag_key = GSL_INTEG_GAUSS41;
    this->intConfigs[0] = i1;
    Common::IntegrationConfig* i2 = new Common::IntegrationConfig;
    i2->method = "gsl_monte_vegas_integrate";
    i2->MC_warmupCalls = 5000;
    i2->calls = 50000;
    this->intConfigs[1] = i2;
}

FullyInclusiveLeptoProduction::~FullyInclusiveLeptoProduction() {
}

Common::IntegrationConfig* FullyInclusiveLeptoProduction::getIntegrationConfig(const str& method) const {
    // partonic functions
    if (boost::starts_with(method,"cg0_") || boost::starts_with(method,"cg0t_") || 
        boost::starts_with(method,"cg1_") || boost::starts_with(method,"cg1t_") || boost::starts_with(method,"cg1tp_") || boost::starts_with(method,"cg1b_") ||
        boost::starts_with(method,"cgBarF1_") || boost::starts_with(method,"cgBar1_") || boost::starts_with(method,"cgBarR1_") ||
        boost::starts_with(method,"cq1_") || boost::starts_with(method,"cq1t_") || boost::starts_with(method,"cq1tp_") || boost::starts_with(method,"cq1b_") ||
        boost::starts_with(method,"cqBarF1_") || boost::starts_with(method,"dq1_"))
        throw invalid_argument("partonic coefficient functions are (almost) exact!");
    // hadronic functions
    if (boost::iequals(method,"F"))
        return this->intConfigs.at(0);
    // leptonic functions
    if (boost::iequals(method,"sigma"))
        return this->intConfigs.at(1);
    throw invalid_argument(str("unknown method: ")+method);
}

#define intND(N) cdbl FullyInclusiveLeptoProduction::int##N##D() const {\
    this->ker->dim = N;\
    return Common::integrate##N##D<FullyInclusive::IntKer>(FIker,*this->intConfigs.at(N-1),this->intOut);\
}
intND(1)
intND(3)

#define implementCoeffsRaw(ns,n) \
cdbl FullyInclusiveLeptoProduction::n##_VV() const { initPartonicVV this->ker->mode = ns::Mode_##n##_VV; return (*FIker)(); }\
cdbl FullyInclusiveLeptoProduction::n##_VA() const { initPartonicVA this->ker->mode = ns::Mode_##n##_VA; return (*FIker)(); }\
cdbl FullyInclusiveLeptoProduction::n##_AA() const { initPartonicAA this->ker->mode = ns::Mode_##n##_AA; return (*FIker)(); }

#define implementCoeffsCommon(cf) implementCoeffsRaw(Common::AbstractIntKer,cf)
#define implementCoeffsLocal(cf) implementCoeffsRaw(FullyInclusive::IntKer,cf)

implementCoeffsCommon(cg0)
implementCoeffsLocal(cg0t)

implementCoeffsCommon(dq1)
implementCoeffsCommon(cq1)
implementCoeffsLocal(cq1t)
implementCoeffsLocal(cq1tp)
implementCoeffsLocal(cq1b)
implementCoeffsCommon(cqBarF1)

implementCoeffsCommon(cg1)
implementCoeffsLocal(cg1t)
implementCoeffsLocal(cg1tp)
implementCoeffsLocal(cg1b)
implementCoeffsCommon(cgBarF1)
implementCoeffsCommon(cgBarR1)
implementCoeffsCommon(cgBar1)
       
#define checkMu if (0. != this->ker->muF2.cHAQTransverseMomentum || 0. != this->ker->muF2.cHQPairTransverseMomentum  || \
                    0. != this->ker->muR2.cHAQTransverseMomentum || 0. != this->ker->muR2.cHQPairTransverseMomentum)\
        throw domain_error("scales for fully inclusive computation may not depend on any variable scale!");

cdbl FullyInclusiveLeptoProduction::F() const {
    initF
    checkMu
    this->ker->mode = Common::AbstractIntKer::Mode_F;
    dbl r = this->int1D();
    // unset integration variables
    this->ker->s = dblNaN;
    return r;
}

cdbl FullyInclusiveLeptoProduction::sigma() const {
    return 0.;
}