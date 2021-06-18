#include "Common/AbstractLeptoProduction.h"

Common::AbstractLeptoProduction::AbstractLeptoProduction(AbstractIntKer* ker, cuint nlf, cdbl m2){
    this->ker = ker;
    this->setNumberOfLightFlavours(nlf);
    this->setM2(m2);
    this->intOut = new IntegrationOutput();
}

Common::AbstractLeptoProduction::~AbstractLeptoProduction() {
    if (0 != this->ker)
        delete this->ker;
    for (uint j = 0; j < this->intConfigs.size(); ++j)
        delete this->intConfigs[j];
    if (0 != this->intOut)
        delete this->intOut;
}

Common::IntegrationOutput Common::AbstractLeptoProduction::getIntegrationOutput() const {
    return Common::IntegrationOutput(*this->intOut);
}

void Common::AbstractLeptoProduction::setNumberOfLightFlavours(cuint nlf) {
    if (nlf < 3 || nlf > 5)
        throw domain_error("number of light flavours has to be between 3 and 5!");
    this->ker->nlf = nlf;
}

void Common::AbstractLeptoProduction::setM2(cdbl m2) {
    if (m2 <= 0.)
        throw domain_error("mass^2 of heavy flavours has to be positive!");
    this->ker->m2 = m2;
}

void Common::AbstractLeptoProduction::setQ2(cdbl Q2) {
    checkQ2(Q2)
    this->ker->Q2 = Q2;
}

void Common::AbstractLeptoProduction::setProjection(Projection proj) {
    this->ker->proj = proj;
}

void Common::AbstractLeptoProduction::setPartonicS(cdbl s) {
    checkPartonicS(s)
    this->ker->s = s;
}

void Common::AbstractLeptoProduction::setPartonicEta(cdbl eta) {
    if (eta < 0.)
        throw domain_error("partonic eta has to be larger than 0!");
    this->setPartonicS(4.*this->ker->m2*(1.+eta));
}

void Common::AbstractLeptoProduction::setXBjorken(cdbl xBj) {
    checkXBjorken(xBj)
    this->ker->xBj = xBj;
}

void Common::AbstractLeptoProduction::setHadronicS(cdbl Sh) {
    if (Sh < 4.*this->ker->m2)
        throw domain_error("hadronic cm-energy has to be larger than threshold 4m^2!");
    checkQ2(this->ker->Q2)
    this->setXBjorken(this->ker->Q2/(Sh + this->ker->Q2));
}

void Common::AbstractLeptoProduction::setPdf(const str& name, const int member) {
    if (0 != this->ker->pdf)
        delete (this->ker->pdf);
    this->ker->pdf = new PdfWrapper(name,member);
}
    
void Common::AbstractLeptoProduction::setMuR2(const DynamicScaleFactors& muR2) {
    this->ker->muR2 = muR2;
}
    
void Common::AbstractLeptoProduction::setMuF2(const DynamicScaleFactors& muF2) {
    this->ker->muF2 = muF2;
}

void Common::AbstractLeptoProduction::setMu2(const DynamicScaleFactors& mu2) {
    this->setMuF2(mu2);
    this->setMuR2(mu2);
}

void Common::AbstractLeptoProduction::setLambdaQCD(cdbl lambdaQCD) {
    if (!isfinite(lambdaQCD) || lambdaQCD <= 0.)
        throw domain_error("Lambda_QCD has to be set, finite and strict positive!");
    if (0 != this->ker->aS)
        delete (this->ker->aS);
    LHAPDF::AlphaS_Analytic *p = new LHAPDF::AlphaS_Analytic();
    p->setLambda(this->ker->nlf + 1, lambdaQCD);
    this->ker->aS = p;
}

void Common::AbstractLeptoProduction::setAlphaSByLHAPDF(const str& name, const int member) {
    if (0 != this->ker->aS)
        delete (this->ker->aS);
    this->ker->aS = LHAPDF::mkAlphaS(name,member);
}

void Common::AbstractLeptoProduction::setAlphaEM(cdbl alphaEM) {
    checkAlphaEM(alphaEM)
    this->ker->alphaEM = alphaEM;
}

void Common::AbstractLeptoProduction::setPolarizeBeams(const bool polarizeBeams) {
    this->ker->polarizeBeams = polarizeBeams;
}

void Common::AbstractLeptoProduction::setLeptonicS(cdbl Sl) {
    checkLeptonicS(Sl)
    this->ker->Sl = Sl;
}

void Common::AbstractLeptoProduction::setQ2min(cdbl Q2min) {
    if (!isfinite(Q2min) || Q2min <= 0.)
        throw domain_error("Q2min has to be finite and strict positive!");
    this->ker->Q2min = Q2min;
    this->ker->q2minHVQDIS = dblNaN;
}

void Common::AbstractLeptoProduction::setQ2minByHVQDIS(cdbl q2minHVQDIS)  {
    if (!isfinite(q2minHVQDIS) || q2minHVQDIS <= 0.)
        throw domain_error("q2minHVQDIS has to be finite and strict positive!");
    this->ker->q2minHVQDIS = q2minHVQDIS;
    this->ker->Q2min = dblNaN;
}