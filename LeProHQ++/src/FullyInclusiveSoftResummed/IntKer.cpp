#include "FullyInclusiveSoftResummed/IntKer.h"
#include "FullyInclusiveSoftResummed/MellinFuncs.h"
#include "cg0t.h"
#include "FullyInclusiveSoftResummed/ResExp.h"
#include "FullyInclusive/PartonicVars.hpp"

FullyInclusiveSoftResummed::IntKer::IntKer() : AbstractIntKer() {
};

FullyInclusiveSoftResummed::IntKer::~IntKer() {
    if (this->interpolator)
        delete this->interpolator;
};

void FullyInclusiveSoftResummed::IntKer::setPath(cdbl a) {
    // for now used edged path
    cdbl ang = M_PI_2 + this->path_angle;
    this->jac_path = this->path_length*dcmplx(cos(ang), sin(ang));
    this->N = dcmplx(this->path_offset) + a * this->jac_path;
}

cdbl FullyInclusiveSoftResummed::IntKer::runPartonic(cdbl a1) {
    if (this->mode != Mode_cg0t_VV and this->mode != Mode_cg0t_VA and this->mode != Mode_cg0t_AA)
        return 0.;
    this->setPath(a1);
    partonic_rho(m2,s)
    dcmplx res = 1./(dcmplx(0.,M_PI))*this->jac_path*exp(-this->N * log(rho));
    fPtr2dblN fVV, fVA, fAA;
    getME(FullyInclusiveSoftResummed::coeffs::cg0t)
    switch(this->mode) {
        case Mode_cg0t_VV: res *= fVV(this->m2,-this->Q2,this->N); break;
        case Mode_cg0t_VA: res *= fVA(this->m2,-this->Q2,this->N); break;
        case Mode_cg0t_AA: res *= fAA(this->m2,-this->Q2,this->N); break;
    }
    return res.real();
}

FullyInclusiveSoftResummed::cdcmplx FullyInclusiveSoftResummed::IntKer::cg0t() const {
    fPtr2dblN fVV, fVA, fAA;
    getME(FullyInclusiveSoftResummed::coeffs::cg0t)
    cdbl eH = this->getElectricCharge(this->nlf+1);
    cdbl gVQ = this->getVectorialCoupling(this->nlf+1);
    cdbl gAQ = this->getAxialCoupling(this->nlf+1);
    dcmplx lo = 0.;
    if (isParityConservingProj(this->proj)) {
        cdcmplx eVV = fVV(this->m2,-this->Q2,this->N);
        if (this->flags.usePhoton) lo += eH*eH * eVV;
        if (this->flags.usePhotonZ) lo -= this->getNormPhZ() * eH*gVQ * eVV;
        if (this->flags.useZ) {
            cdcmplx eAA = fAA(this->m2,-this->Q2,this->N);;
            lo += this->getNormZ()*(gVQ*gVQ*eVV + gAQ*gAQ*eAA);
        }
    } else {
        cdcmplx eVA = fVA(this->m2,-this->Q2,this->N);
        if (this->flags.usePhotonZ) lo -= this->getNormPhZ() * eH*gAQ * eVA;
        if (this->flags.useZ) lo += this->getNormZ() * 2.*gVQ*gAQ * eVA;
    }
    return lo;
}

cdbl FullyInclusiveSoftResummed::IntKer::runHadronic(cdbl a1) {
    this->setPath(a1);
    dcmplx res = 1./(dcmplx(0.,M_PI))*this->jac_path;
    // N is conjugated to x/z_max
    res *= this->interpolator->evaluateN(this->N, log(this->xBj/this->getZMax()));
    cdbl alphaS = this->getAlphaS(0.);
    cdbl nLO = alphaS/(4.*M_PI*M_PI) * (this->Q2/this->m2);
    dcmplx lnDelta = 0.;
    cdbl beta0 = this->beta0lf();
    cdbl Ag1 = Color::CA;
    cdcmplx lnN = log(N);
    cdcmplx lambda = alphaS * beta0 * lnN;
    lnDelta += lnN * gLL(Ag1, beta0, lambda);
    res *= nLO * this->cg0t() * exp(lnDelta);
    return res.real();
}

cdbl FullyInclusiveSoftResummed::IntKer::operator()(cdbl a1, cdbl a2, cdbl a3) {      
    // void mode
    if (0 == this->mode) return 0.;
    // partonic mode?
    if (this->mode < Mode_F) {
        cdbl r = this->runPartonic(a1);
        return isfinite(r) ? r : 0.;
    }
    // hadronic mode?
    if (this->mode < Mode_sigma) {
        cdbl r = this->runHadronic(a1);
        return isfinite(r) ? r : 0.;
    }
    /*// leptonic mode = Mode_sigma
    cdbl r = this->runLeptonic(a1, a2, a3);
    return isfinite(r) ? r : 0.;*/
    return 0.;
}

void FullyInclusiveSoftResummed::IntKer::operator()(const double x[], const int k[], const double& weight, const double aux[], double f[]) {
    // proide simple redirect as weight isn't needed
    cdbl a0 = this->dim >= 1 ? x[0] : 0.;
    cdbl a1 = this->dim >= 2 ? x[1] : 0.;
    cdbl a2 = this->dim >= 3 ? x[2] : 0.;
    f[0] = this->operator()(a0,a1,a2);
}