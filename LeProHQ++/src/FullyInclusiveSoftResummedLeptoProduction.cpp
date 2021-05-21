#include "FullyInclusiveSoftResummedLeptoProduction.h"

#include <gsl/gsl_integration.h>
#include <boost/algorithm/string.hpp>

#include "gslpp/gslpp.Functor.hpp"

#include "Common/Integration.hpp"

FullyInclusiveSoftResummedLeptoProduction::FullyInclusiveSoftResummedLeptoProduction(cuint nlf, cdbl m2) :
    AbstractLeptoProduction(new FullyInclusiveSoftResummed::IntKer(), nlf, m2) {
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

FullyInclusiveSoftResummedLeptoProduction::~FullyInclusiveSoftResummedLeptoProduction() {
}

Common::IntegrationConfig* FullyInclusiveSoftResummedLeptoProduction::getIntegrationConfig(const str& method) const {
    // partonic functions
    if (boost::starts_with(method,"cg0_") || boost::starts_with(method,"cgBarR1_") ||
        boost::starts_with(method,"cg1_") || boost::starts_with(method,"cgBarF1_") || boost::starts_with(method,"cgBar1_") || 
        boost::starts_with(method,"cq1_") || boost::starts_with(method,"cqBarF1_") || boost::starts_with(method,"dq1_"))
        throw invalid_argument("partonic coefficient function is not available!");
    if (boost::iequals(method,"cg0t_"))
        return this->intConfigs.at(0);
    // hadronic functions
    if (boost::iequals(method,"F"))
        return this->intConfigs.at(0);
    // leptonic functions
    if (boost::iequals(method,"sigma"))
        return this->intConfigs.at(1);
    throw invalid_argument(str("unknown method: ")+method);
}

void FullyInclusiveSoftResummedLeptoProduction::setPath(cdbl offset, cdbl angle, cdbl length) {
    FISRker->path_offset = offset;
    FISRker->path_angle = angle;
    FISRker->path_length = length;
}

#define intND(N) cdbl FullyInclusiveSoftResummedLeptoProduction::int##N##D() const {\
    this->ker->dim = N;\
    return Common::integrate##N##D<FullyInclusiveSoftResummed::IntKer>(FISRker,*this->intConfigs.at(N-1),this->intOut);\
}
intND(1)
intND(3)

void FullyInclusiveSoftResummedLeptoProduction::setInterpolation(FullyInclusiveSoftResummed::Interpolation::vdbl xgrid, cuint degree, bool isLog) {
    if (FISRker->interpolator)
        delete FISRker->interpolator;
    FISRker->interpolator = new FullyInclusiveSoftResummed::Interpolation::Interpolator(xgrid, degree, isLog);
}

#define checkMu if (0. != this->ker->muF2.cHAQTransverseMomentum || 0. != this->ker->muF2.cHQPairTransverseMomentum  || \
                    0. != this->ker->muR2.cHAQTransverseMomentum || 0. != this->ker->muR2.cHQPairTransverseMomentum)\
        throw domain_error("scales for fully inclusive computation may not depend on any variable scale!");
        
#define initInterpolator if (!FISRker->interpolator) throw domain_error("Need an interpolation config!");\
    FISRker->interpolator->loadPdfData(this->ker->pdf, 21, this->ker->getScale(this->ker->muF2, 0.0, 0.0));

cdbl FullyInclusiveSoftResummedLeptoProduction::cg0t_VV() const {
    initPartonicVV
    this->ker->mode = FullyInclusiveSoftResummed::IntKer::Mode_cg0t_VV;
    dbl r = this->int1D();
    // unset integration variables
    FISRker->N = dblNaN;
    return r;
}

cdbl FullyInclusiveSoftResummedLeptoProduction::F() const {
    initF
    checkMu
    if (!isParityConservingProj(this->ker->proj))
        return 0.;
    initInterpolator
    
    this->ker->mode = Common::AbstractIntKer::Mode_F;
    dbl r = this->int1D();
    // unset integration variables
    FISRker->N = dblNaN;
    return r;
}

cdbl FullyInclusiveSoftResummedLeptoProduction::sigma() const {
    return 0.;
}