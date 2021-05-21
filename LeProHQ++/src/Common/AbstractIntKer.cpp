#include "Common/AbstractIntKer.h"

#include "./ME/BQED.h"
#include "./AltarelliParisi.hpp"

Common::AbstractIntKer::AbstractIntKer() : flags(true, true, true, true, true, true, true),
    muR2(4.,1.,0.,0.), muF2(4.,1.,0.,0.) {
    this->aS = new LHAPDF::AlphaS_Analytic();
}

Common::AbstractIntKer::~AbstractIntKer() {
    if (0 != this->pdf)
        delete this->pdf;
    if (0 != this->aS)
        delete this->aS;
}

cdbl Common::AbstractIntKer::getElectricCharge(cint PDGId) const {
    // see https://en.wikipedia.org/wiki/Mathematical_formulation_of_the_Standard_Model
    // + Leader + PDG
    switch(PDGId) {
        case 11: case -11: return -1.; // e
        case 1: case -1: case 3: case -3: case 5: case -5: return -1./3.; // d,s,b
        case 2: case -2: case 4: case -4: case 6: case -6: return  2./3.; // u,c,t
        case 21: return 0.; // g
        default:
            throw domain_error((boost::format("unkown particle id: %d")%PDGId).str());
    }
}

cdbl Common::AbstractIntKer::getWeakIsospin(cint PDGId) const {
    // see this->getElectricCharge
    switch (PDGId) {
        case  11: case -11: case 1: case -1: case 3: case -3: case 5: case -5: return -.5; // e,d,s,b
                            case 2: case -2: case 4: case -4: case 6: case -6: return +.5; //   u,s,c
        case 21: return 0.; // g
        default:
            throw domain_error((boost::format("unkown particle id: %d")%PDGId).str());
    }
}

cdbl Common::AbstractIntKer::getVectorialCoupling(cint PDGId) const {
    return this->getWeakIsospin(PDGId) - 2.*this->getElectricCharge(PDGId)*this->sin2ThetaWeak;
}

cdbl Common::AbstractIntKer::getAxialCoupling(cint PDGId) const {
    return this->getWeakIsospin(PDGId);
}

cdbl Common::AbstractIntKer::getScale(const DynamicScaleFactors& factors, cdbl HAQTransverseMomentum, cdbl HQPairTransverseMomentum) const {
    cdbl mu2 = factors.cM2 * this->m2 + factors.cQ2 * this->Q2 + factors.cHAQTransverseMomentum * HAQTransverseMomentum*HAQTransverseMomentum + factors.cHQPairTransverseMomentum * HQPairTransverseMomentum*HQPairTransverseMomentum;
    if (!isfinite(mu2) || mu2 <= 0.)
        throw domain_error((boost::format("all scales have to be finite and positive! (%g = %g*%g + %g*%g + %g*%g^2 + %g*%g^2)")%mu2%factors.cM2%m2%factors.cQ2%Q2%factors.cHAQTransverseMomentum%HAQTransverseMomentum%factors.cHQPairTransverseMomentum%HQPairTransverseMomentum).str());
    return mu2;
}

cdbl Common::AbstractIntKer::getAlphaS(cdbl HAQTransverseMomentum, cdbl HQPairTransverseMomentum) const {
    cdbl muR2 = this->getScale(this->muR2, HAQTransverseMomentum, HQPairTransverseMomentum);
    this->aS->setOrderQCD(1 + (this->flags.useNextToLeadingOrder ? 1 : 0));
    return this->aS->alphasQ2(muR2);
}

void Common::AbstractIntKer::getBQED(fPtr4dbl &fVV, fPtr4dbl &fVA, fPtr4dbl &fAA) const { getME(Common::ME::BQED) }

Common::AbstractIntKer::fPtr1dbl Common::AbstractIntKer::getPgq0() const {
    switch(this->proj) {
        case F2:   case FL: case xF3: return &AltarelliParisi::Pgq0;
        case x2g1: case g4: case gL:  return &AltarelliParisi::DeltaPgq0;
    }
    throw domain_error("unkonwn projection: "+strProj(this->proj));
}

Common::AbstractIntKer::fPtr1dbl Common::AbstractIntKer::getPgq1() const {
    switch(this->proj) {
        case F2:   case FL: case xF3: return &AltarelliParisi::Pgq1;
        case x2g1: case g4: case gL:  return &AltarelliParisi::DeltaPgq1;
    }
    throw domain_error("unkonwn projection: "+strProj(this->proj));
}

Common::AbstractIntKer::fPtr1dbl Common::AbstractIntKer::getPggH0() const {
    switch(this->proj) {
        case F2:   case FL: case xF3: return &AltarelliParisi::PggH0;
        case x2g1: case g4: case gL:  return &AltarelliParisi::DeltaPggH0;
    }
    throw domain_error("unkonwn projection: "+strProj(this->proj));
}

Common::AbstractIntKer::fPtr1dbl Common::AbstractIntKer::getPggH1() const {
    switch(this->proj) {
        case F2:   case FL: case xF3: return &AltarelliParisi::PggH1;
        case x2g1: case g4: case gL:  return &AltarelliParisi::DeltaPggH1;
    }
    throw domain_error("unkonwn projection: "+strProj(this->proj));
}

cdbl Common::AbstractIntKer::getNormPhZ() const {
    cdbl etaPhZ = 1./(4. * this->sin2ThetaWeak * (1. - this->sin2ThetaWeak)) * this->Q2 / (this->Q2 + this->MZ2);
    cdbl sgn = this->incomingLeptonHasNegativeCharge ? -1. : 1.;
    cdbl lam = this->incomingLeptonHasPositiveHelicity ? 1. : -1.;
    cdbl gVl = this->getVectorialCoupling(11); /// @todo initial lepton defaults to electron = 11; allow change?
    cdbl gAl = this->getAxialCoupling(11);
    if (isParityConservingProj(this->proj))
        return (gVl + sgn*lam*gAl) * etaPhZ;
    return (gAl + sgn*lam*gVl) * etaPhZ;
}

cdbl Common::AbstractIntKer::getNormZ() const {
    cdbl etaPhZ = 1./(4. * this->sin2ThetaWeak * (1. - this->sin2ThetaWeak)) * this->Q2 / (this->Q2 + this->MZ2);
    cdbl sgn = this->incomingLeptonHasNegativeCharge ? -1. : 1.;
    cdbl lam = this->incomingLeptonHasPositiveHelicity ? 1. : -1.;
    cdbl gVl = this->getVectorialCoupling(11);
    cdbl gAl = this->getAxialCoupling(11);
    if (isParityConservingProj(this->proj))
        return (gVl*gVl + gAl*gAl + 2.*sgn*lam*gVl*gAl) * etaPhZ*etaPhZ;
    return (2.*gVl*gAl + sgn*lam*(gVl*gVl + gAl*gAl)) * etaPhZ*etaPhZ;
}