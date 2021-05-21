#include "FullyInclusive/IntKer.h"

#include "coeffs/cg0.h"
#include "coeffs/cg0t.h"
#include "coeffs/dq1.h"
#include "coeffs/cq1.h"
#include "coeffs/cqBarF1.h"
#include "coeffs/cg1.h"
#include "coeffs/cgBar1.h"

FullyInclusive::IntKer::IntKer() : AbstractIntKer() {
};

FullyInclusive::IntKer::~IntKer() {
};

#define getCoeffs(cf)\
  fPtr3dbl fVV, fVA, fAA;\
  switch(this->proj) {\
    case F2:   fVV = &FullyInclusive::coeffs::cf##_F2_VV;   fAA = &FullyInclusive::coeffs::cf##_F2_AA;   fVA = 0; break;\
    case FL:   fVV = &FullyInclusive::coeffs::cf##_FL_VV;   fAA = &FullyInclusive::coeffs::cf##_FL_AA;   fVA = 0; break;\
    case x2g1: fVV = &FullyInclusive::coeffs::cf##_x2g1_VV; fAA = &FullyInclusive::coeffs::cf##_x2g1_AA; fVA = 0; break;\
    case xF3:  fVA = &FullyInclusive::coeffs::cf##_xF3_VA; fVV = 0; fAA = 0;              break;\
    case g4:   fVA = &FullyInclusive::coeffs::cf##_g4_VA;  fVV = 0; fAA = 0;              break;\
    case gL:   fVA = &FullyInclusive::coeffs::cf##_gL_VA;  fVV = 0; fAA = 0;              break;\
    default: throw domain_error("unkonwn projection: "+strProj(this->proj));\
  }

#define getCoeff(cf) {\
  getCoeffs(cf)\
  if (Mode_##cf##_VV == this->mode) return fVV(this->m2,-this->Q2,this->s);\
  if (Mode_##cf##_VA == this->mode) return fVA(this->m2,-this->Q2,this->s);\
  if (Mode_##cf##_AA == this->mode) return fAA(this->m2,-this->Q2,this->s);\
}


/**
 * @brief implement partonic coefficient as all-currents-combination
 * @param n name
 */
#define implementFIPartonicCoeff(n)\
cdbl FullyInclusive::IntKer::n() const {\
    getCoeffs(n)\
    combineCurs(fVV(this->m2,-this->Q2,this->s),fVA(this->m2,-this->Q2,this->s),fAA(this->m2,-this->Q2,this->s))\
}

implementFIPartonicCoeff(cg0)
implementFIPartonicCoeff(cg0t)
implementFIPartonicCoeff(cg1)
implementFIPartonicCoeff(cgBar1)
implementFIPartonicCoeff(cq1)
implementFIPartonicCoeff(cqBarF1)

cdbl FullyInclusive::IntKer::cgBarF1() const {
    return this->cgBar1() - this->cgBarR1();
}

cdbl FullyInclusive::IntKer::cgBarR1() const {
    cdbl b = this->beta0lf()/(16.*M_PI*M_PI);
    return b * this->cg0();
}

cdbl FullyInclusive::IntKer::runPartonic() {
    // cg0 will be used recursively so it has to come first
    getCoeff(cg0)
    getCoeff(cg0t)
    // cgBarR1 ~ cg0
    if (Mode_cgBarR1_VV == this->mode || Mode_cgBarR1_VA == this->mode || Mode_cgBarR1_AA == this->mode) {
        cuint oldMode = this->mode;
        if (Mode_cgBarR1_VV == this->mode) this->mode = Mode_cg0_VV;
        else if (Mode_cgBarR1_VA == this->mode) this->mode = Mode_cg0_VA;
        else if (Mode_cgBarR1_AA == this->mode) this->mode = Mode_cg0_AA;
        cdbl b = this->beta0lf()/(16.*M_PI*M_PI);
        cdbl cg0 = this->runPartonic();
        this->mode = oldMode;
        return b * cg0;
    }
    getCoeff(dq1)
    // cq1 and friends
    getCoeff(cq1t)
    getCoeff(cq1tp)
    getCoeff(cq1b)
    getCoeff(cq1)
    getCoeff(cqBarF1)
    // cg1 and friends
    getCoeff(cg1t)
    getCoeff(cg1tp)
    getCoeff(cg1b)
    getCoeff(cg1)
    getCoeff(cgBar1)
    // cgBarF1
    if (Mode_cgBarF1_VV == this->mode || Mode_cgBarF1_VA == this->mode || Mode_cgBarF1_AA == this->mode) {
        cuint oldMode = this->mode;
        uint jmode = 0;
        uint rmode = 0;
        if (Mode_cgBarF1_VV == this->mode) { jmode = Mode_cgBar1_VV; rmode = Mode_cgBarR1_VV; }
        else if (Mode_cgBarF1_VA == this->mode) { jmode = Mode_cgBar1_VA; rmode = Mode_cgBarR1_VA; }
        else if (Mode_cgBarF1_AA == this->mode) { jmode = Mode_cgBar1_AA; rmode = Mode_cgBarR1_VA; }
        this->mode = jmode;
        cdbl cgBar1 = this->runPartonic();
        this->mode = rmode;
        cdbl cgBarR1 = this->runPartonic();
        this->mode = oldMode;
        return cgBar1 - cgBarR1;
    }
    return 0.;
}

void FullyInclusive::IntKer::setPartonicVars(cdbl a) {
    cdbl zmax = this->getZMax();
    this->jac_z = (zmax - this->xBj);
    cdbl z = this->xBj + a * this->jac_z;
    cdbl sp = this->Q2 / z;
    this->xi = this->xBj / z;
    this->s = sp - this->Q2;
}

cdbl FullyInclusive::IntKer::runHadronic(cdbl a1) {
    this->setPartonicVars(a1);
    dbl r = 0.;
    cdbl jac = this->xi/this->xBj * this->jac_z; // = jac/z
    // LO
    if (this->flags.useLeadingOrder && this->flags.useGluonicChannel) {
        r += jac*this->Fg0();
    }
    // NLO
    if (this->flags.useNextToLeadingOrder) {
        if (this->flags.useQuarkChannel) {
            r += jac*this->Fq1();
        }
        if (this->flags.useGluonicChannel) {
            r += jac*this->Fg1();
        }
    }
    return r;
}

cdbl FullyInclusive::IntKer::Fg0() const {
    cdbl alphaS = this->getAlphaS(0.);
    cdbl curMuF2 = this->getScale(this->muF2,0.);
    cdbl g = this->pdf->xfxQ2(21,this->xi,curMuF2);
    cdbl nLO = alphaS/(4.*M_PI*M_PI) * (this->Q2/this->m2);
    cdbl r = nLO * g * this->cg0();
    return r;
}

cdbl FullyInclusive::IntKer::Fg1() const {
    cdbl alphaS = this->getAlphaS(0.);
    cdbl curMuF2 = this->getScale(this->muF2,0.);
    cdbl curMuR2 = this->getScale(this->muR2,0.);
    cdbl g = this->pdf->xfxQ2(21,this->xi,curMuF2);
    cdbl nNLO = alphaS*alphaS/(M_PI) * (this->Q2/this->m2);
    return nNLO * g * (this->cg1() + log(curMuF2/m2)*this->cgBarF1() + log(curMuR2/m2)*this->cgBarR1());
}

cdbl FullyInclusive::IntKer::Fq1() const {
    // compute matrix elements for dq1
    // as they are combined with the properties of the light quarks,
    // we have to keep the raw-current elements
    dbl dq1_VV = 0.,dq1_VA = 0.,dq1_AA = 0.;
    {
        getCoeffs(dq1)
        if (isParityConservingProj(this->proj)) {
            dq1_VV = fVV(this->m2,-this->Q2,this->s);
            dq1_AA = fAA(this->m2,-this->Q2,this->s);
        } else {
            dq1_VA = fVA(this->m2,-this->Q2,this->s);
        }
    }
    // common stuff
    cdbl alphaS = this->getAlphaS(0.);
    cdbl curMuF2 = this->getScale(this->muF2,0.);
    cdbl nNLO = alphaS*alphaS/(M_PI) * (this->Q2/this->m2);
    // join elems
    cdbl cq1j = this->cq1() + log(curMuF2/this->m2) * this->cqBarF1();
    dbl fqs = 0.;
    for (uint q = 1; q < this->nlf + 1; ++q) {
        dbl dq1 = 0.;
        cdbl eL = getElectricCharge(q);
        cdbl gVq = this->getVectorialCoupling(q);
        cdbl gAq = this->getAxialCoupling(q);
        if (isParityConservingProj(this->proj)) {
            if (this->flags.usePhoton) {
                dq1 += eL*eL * dq1_VV; 
            }
            if (this->flags.usePhotonZ) {
                dq1 -= this->getNormPhZ() *  eL*gVq           * dq1_VV;
            }
            if (this->flags.useZ) {
                dq1 += this->getNormZ()*(gVq*gVq*dq1_VV + gAq*gAq*dq1_AA);
            }
        } else {
            if (this->flags.usePhotonZ) {
                dq1 -= this->getNormPhZ() *  eL*gAq           * dq1_VA;
            }
            if (this->flags.useZ) {
                dq1 += this->getNormZ() * 2.*gVq*gAq          * dq1_VA;
            }
        }
        // combine cq1+cqBarF+dq1
        fqs += (this->pdf->xfxQ2((int)q,this->xi,curMuF2) + this->pdf->xfxQ2(-((int)q),this->xi,curMuF2))*(cq1j + dq1);
    }
    return nNLO*fqs;
}

cdbl FullyInclusive::IntKer::operator()(cdbl a1, cdbl a2, cdbl a3) {      
    // void mode
    if (0 == this->mode) return 0.;
    // partonic mode?
    if (this->mode < Mode_F) {
        cdbl r = this->runPartonic();
        return isfinite(r) ? r : 0.;
    }
    // hadronic mode?
    if (this->mode < Mode_sigma) {
        cdbl r = this->runHadronic(a1);
        return isfinite(r) ? r : 0.;
    }
    /*// leptonic mode = Mode_sigma
    cdbl r = this->runLeptonic(a1, a2, a3, a4, a5);
    return isfinite(r) ? r : 0.;*/
    return 0.;
}

void FullyInclusive::IntKer::operator()(const double x[], const int k[], const double& weight, const double aux[], double f[]) {
    // proide simple redirect as weight isn't needed
    cdbl a0 = this->dim >= 1 ? x[0] : 0.;
    cdbl a1 = this->dim >= 2 ? x[1] : 0.;
    cdbl a2 = this->dim >= 3 ? x[2] : 0.;
    f[0] = this->operator()(a0,a1,a2);
}