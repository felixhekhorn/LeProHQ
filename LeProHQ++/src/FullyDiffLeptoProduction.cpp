#include "FullyDiffLeptoProduction.h"

#include <gsl/gsl_integration.h>
#include <boost/algorithm/string.hpp>

#include "gslpp/gslpp.Functor.hpp"

#include "Common/Integration.hpp"

FullyDiffLeptoProduction::FullyDiffLeptoProduction(cuint nlf, cdbl m2, cdbl xTilde, cdbl omega, cdbl deltax, cdbl deltay) :
    AbstractFixedOrderLeptoProduction(new FullyDiff::IntKer(), nlf, m2) {
    this->setXTilde(xTilde);
    this->setOmega(omega);
    this->setDeltax(deltax);
    this->setDeltay(deltay);
    
    // setup default values for IntegrationConfig
    this->intConfigs.resize(5);
    Common::IntegrationConfig* i1 = new Common::IntegrationConfig;
    i1->method = "gsl_integration_qag";
    i1->calls = 10000;
    i1->GslQag_epsabs = 5e-6;
    i1->GslQag_epsrel = 1e-4;
    i1->GslQag_key = GSL_INTEG_GAUSS41;
    this->intConfigs[0] = i1;
    Common::IntegrationConfig* i2 = new Common::IntegrationConfig;
    i2->method = "Dvegas";
    i2->MC_warmupCalls = 3000;
    i2->calls = 30000;
    this->intConfigs[1] = i2;
    Common::IntegrationConfig* i3 = new Common::IntegrationConfig;
    i3->method = "Dvegas";
    i3->MC_warmupCalls = 5000;
    i3->calls = 50000;
    this->intConfigs[2] = i3;
    Common::IntegrationConfig* i4 = new Common::IntegrationConfig;
    i4->method = "Dvegas";
    i4->MC_warmupCalls = 7000;
    i4->calls = 70000;
    this->intConfigs[3] = i4;
    Common::IntegrationConfig* i5 = new Common::IntegrationConfig;
    i5->method = "Dvegas";
    i5->MC_warmupCalls = 10000;
    i5->calls = 100000;
    this->intConfigs[4] = i5;
}

FullyDiffLeptoProduction::~FullyDiffLeptoProduction() {
}

void FullyDiffLeptoProduction::setXTilde(cdbl xTilde) {
    if (xTilde <= 0. || xTilde >= 1.)
        throw domain_error((boost::format("xTilde (%e) has to be within (0,1)!")%xTilde).str());
    FDker->xTilde = xTilde;
}

void FullyDiffLeptoProduction::setOmega(cdbl omega) {
    if (omega <= 0. || omega >= 2.)
        throw domain_error((boost::format("omega (%e) has to be within (0,2)!")%omega).str());
    if (FDker->deltay >= omega)
        throw domain_error((boost::format("omega (%e) has to be bigger than deltay (%e)!")%omega%(FDker->deltay)).str());
    FDker->omega = omega;
}

void FullyDiffLeptoProduction::setDeltax(cdbl deltax) {
    if (deltax <= 0 || deltax >= 1)
        throw domain_error((boost::format("deltax (%e) has to be positive and smaller then 1!")%deltax).str());
    FDker->deltax = deltax;
}

void FullyDiffLeptoProduction::setDeltay(cdbl deltay) {
    if (deltay <= 0 || deltay >= 2 || deltay >= FDker->omega)
        throw domain_error((boost::format("deltay (%e) has to be positive, smaller than 2 and smaller than omega (%e)!")%deltay%(FDker->omega)).str());
    FDker->deltay = deltay;
}

Common::IntegrationConfig* FullyDiffLeptoProduction::getIntegrationConfig(const str& method) const {
    // partonic functions
    if (boost::starts_with(method,"cg0_") || boost::starts_with(method,"cgBarR1_"))
        return this->intConfigs.at(1-1);
    if (boost::starts_with(method,"cgBarF1_") || boost::starts_with(method,"cgBar1_") || 
        boost::starts_with(method,"cqBarF1_"))
        return this->intConfigs.at(2-1);
    if (boost::starts_with(method,"cg1_") || boost::starts_with(method,"cq1_") || boost::starts_with(method,"dq1_"))
        return this->intConfigs.at(4-1);
    // hadronic functions
    bool isNLO = this->ker->flags.useNextToLeadingOrder;
    if (boost::iequals(method,"F"))
        return isNLO ? this->intConfigs.at(5-1) : this->intConfigs.at(2-1);
    /*// leptonic functions
    if (boost::iequals(method,"sigma"))
        return isNLO ? this->intConfigs.at(7-1) : this->intConfigs.at(4-1);
    */throw invalid_argument(str("unknown method: ")+method);
}

#define intND(N) cdbl FullyDiffLeptoProduction::int##N##D() const {\
    this->ker->dim = N;\
    return Common::integrate##N##D<FullyDiff::IntKer>(FDker,*this->intConfigs.at(N-1),this->intOut);\
}
intND(1)
intND(2)
//intND(3)
intND(4)
intND(5)

// implement partonic coefficient functions with d-dim integration
#define implementCoeffs(n,d) \
cdbl FullyDiffLeptoProduction::n##_VV() const { initPartonicVV this->ker->mode = Common::AbstractIntKer::Mode_##n##_VV; return this->int##d##D(); }\
cdbl FullyDiffLeptoProduction::n##_VA() const { initPartonicVA this->ker->mode = Common::AbstractIntKer::Mode_##n##_VA; return this->int##d##D(); }\
cdbl FullyDiffLeptoProduction::n##_AA() const { initPartonicAA this->ker->mode = Common::AbstractIntKer::Mode_##n##_AA; return this->int##d##D(); }

implementCoeffs(cg0,1)
implementCoeffs(cgBarR1,1)

implementCoeffs(cgBarF1,2)
implementCoeffs(cgBar1,2)
implementCoeffs(cqBarF1,2)

implementCoeffs(dq1,4)
implementCoeffs(cq1,4)
implementCoeffs(cg1,4)

void FullyDiffLeptoProduction::activateHistogram(const FullyDiff::histT t, cuint size, const str& path, cdbl min, cdbl max) {
    // assert existance of path
    boost::filesystem::path fp (path);
    boost::filesystem::path par = fp.parent_path();
    if(!boost::filesystem::exists(par))
        throw ios::failure("histogram path \""+par.string()+"\" does not exist!");
    // create
    gslpp::Histogram* h = new gslpp::Histogram(size);
    if (!isnan(min) && !isnan(max))
        h->setRangesUniform(min,max);
    h->setPath(path);
    // insert
    FDker->histMap.insert({t,h});
}
    
void FullyDiffLeptoProduction::setupHistograms() const {
    // hadronic S
    cdbl Sh = this->ker->getHadronicS();
    for (FullyDiff::histMapT::const_iterator it = FDker->histMap.cbegin(); it != FDker->histMap.cend(); ++it) {
        if (it->second->isInitialized())
            continue;
        switch (it->first) {
            case FullyDiff::histT::log10z:                 it->second->setRangesLog10(this->ker->xBj,this->ker->getZMax());      break;
            case FullyDiff::histT::log10xi:                it->second->setRangesLog10(this->ker->xBj/this->ker->getZMax(),1.);   break;
            case FullyDiff::histT::Theta1:
            case FullyDiff::histT::Theta2:                 it->second->setRangesUniform(0.,M_PI);                           break;
            
            case FullyDiff::histT::HQPairInvMass:          it->second->setRangesUniform(2.*sqrt(this->ker->m2),sqrt(Sh));   break;
            case FullyDiff::histT::HQPairDeltaPhi:         it->second->setRangesUniform(0.,M_PI+1e-5);                      break;
            case FullyDiff::histT::HQPairTransverseMomentum:
                /** @todo defalt upper limit of HQPairTransverseMomentum is actually smaller, than 2p_{2,t,max} */
                it->second->setRangesUniform(0.,2.*this->ker->getHAQTransverseMomentumMax());
                break;
            case FullyDiff::histT::HQPairConeSizeVariable:
                /** @todo determine default upper limit of HQPairConeSizeVariable */
                throw domain_error("default upper limit of HQPairConeSizeVariable currently unknown!");
                break;
            
            case FullyDiff::histT::HAQRapidity:
                {cdbl y0 = this->ker->getHAQRapidityMax();
                it->second->setRangesUniform(-y0,y0);}
                break;
            case FullyDiff::histT::HAQTransverseMomentum:         it->second->setRangesUniform(0.,this->ker->getHAQTransverseMomentumMax());        break;
            case FullyDiff::histT::HAQTransverseMomentumScaling:  it->second->setRangesLog10(1.e-3,1.);                        break;
            case FullyDiff::histT::HAQFeynmanX:                   it->second->setRangesUniform(-1.,1.+1e-5);                   break;
            
            case FullyDiff::histT::x:                      it->second->setRangesLog10(this->ker->xBj/this->ker->getZMax(),1.+1e-5); break;
            case FullyDiff::histT::y:                      it->second->setRangesUniform(-1.,1.+1e-5);                               break;
            default: continue;
        }
    }
}

cdbl FullyDiffLeptoProduction::F() const { 
    initF
    this->ker->aS->setOrderQCD(1 + (this->ker->flags.useNextToLeadingOrder ? 1 : 0 ));
    this->ker->mode = Common::AbstractIntKer::Mode_F;
    // Dvegas?
    if (!FDker->histMap.empty() && "Dvegas" != this->getIntegrationConfig("F")->method)
        throw std::invalid_argument("To generate histograms, integration method have to supply integration weights! i.e. currently only Dvegas is available");
    this->setupHistograms();
    // run
    dbl r = dblNaN;
    if (this->ker->flags.useNextToLeadingOrder) {
        r = this->int5D();
    } else {
        r = this->int2D();
    }
    // write histograms
    for (FullyDiff::histMapT::const_iterator it = FDker->histMap.cbegin(); it != FDker->histMap.cend(); ++it)
        it->second->writeToFile();

    // unset integration variables
    this->ker->s = dblNaN;
    return r;
}

cdbl FullyDiffLeptoProduction::sigma() const { return 0.;/*
    checkAlphaS(this->ker->aS)
    checkAlphaEM(this->ker->alphaEM)
    if (0 == this->ker->pdf)
        throw domain_error("needs PDF for leptonic cross section!");
    checkMu
    checkLeptonicS(this->ker->Sl)
    if (!isfinite(this->ker->Q2min) && !isfinite(this->ker->q2minHVQDIS))
        throw domain_error("needs cut on Q2 for leptonic cross section!");
    // kinematic phase space available for y?
    cdbl yBjMax = 1.;
    cdbl yBjMin = 4.*this->ker->m2 / this->ker->Sl;
    if (yBjMin >= yBjMax)
        return 0.;
    // kinematic phase space available for Q2?
    cdbl Q2max = this->ker->Sl - 4.*this->ker->m2;
    if (isfinite(this->ker->Q2min) && this->ker->Q2min >= Q2max)
        return 0.;
    else if (isfinite(this->ker->q2minHVQDIS)) {
        cdbl Q2min = this->ker->q2minHVQDIS * yBjMin*yBjMin/(1. - yBjMin)/(1. - yBjMin);
        if (Q2min >= Q2max)
            return 0.;
    }
    this->ker->mode = Common::AbstractIntKer::Mode_sigma;
    dbl r = dblNaN;
    if (this->flags().useNextToLeadingOrder) {
        r = this->int7D();
    } else {
        r = this->int4D();
    }
    // unset integration variables
    this->ker->xBj = dblNaN;
    this->ker->Q2 = dblNaN;
    this->ker->s = dblNaN;
    return r;*/
}