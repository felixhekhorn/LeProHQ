#include "Inclusive/IntKer.h"

#include "Common/Color.hpp"
#include "Common/ME/BQED.h"
#include "ME/IntA3.h"
#include "ME/IntA2.h"
#include "ME/IntA1.h"
#include "ME/IntR.h"
#include "ME/SV.h"

Inclusive::IntKer::IntKer() : AbstractIntKer() {
};

Inclusive::IntKer::~IntKer() {
};

void Inclusive::IntKer::setT1(cdbl a) {
    _sp
    cdbl beta = this->beta();
    cdbl t1min = -sp/2. * (1. + beta);
    cdbl t1max = -sp/2. * (1. - beta);
    this->jac_t1 = t1max - t1min;
    this->t1 = t1min + this->jac_t1*a;
}

void Inclusive::IntKer::setHAQTransverseMomentum(cdbl a) {
    dbl mt2Max = dblNaN;
    if(isfinite(this->HAQRapidity)) {
        cdbl Sh = this->getHadronicS();
        cdbl coshy = cosh(this->HAQRapidity);
        mt2Max = Sh/(4.*coshy*coshy);
    } else {
        cdbl ptmax = this->getHAQTransverseMomentumMax();
        mt2Max = this->m2 + ptmax*ptmax;
    }
    mt2Max = max(this->m2, mt2Max);
    this->V_HAQTransverseMass2 = mt2Max - this->m2;
    cdbl mt2 = this->m2 + this->V_HAQTransverseMass2*a;
    this->HAQTransverseMomentum = sqrt(mt2 - this->m2);
}

void Inclusive::IntKer::setHAQRapidity(cdbl a) {
    dbl y0 = dblNaN;
    if (isfinite(this->HAQTransverseMomentum)) {
        cdbl mt = sqrt(this->m2 + this->HAQTransverseMomentum*this->HAQTransverseMomentum);
        y0 = acosh(sqrt(this->getHadronicS())/(2.*mt));
    } else {
        y0 = this->getHAQRapidityMax();
    }
    this->V_HAQRapidity = 2.*y0;
    this->HAQRapidity = -y0 + this->V_HAQRapidity*a;
}

void Inclusive::IntKer::setPartonicVars() {
    // set partonic vars
    cdbl Shp = this->getHadronicS() + this->Q2;
    cdbl T1 = this->getHadronicT1();
    cdbl U1 = this->getHadronicU1();
    this->xi = min((this->s4 - U1)/(Shp + T1),1.);
    this->s = Shp*this->xi - this->Q2;
    this->t1 = this->xi*T1;
}

// setup cg0
#define init_cg0 _sp\
    cdbl n = 4.*M_PI*m2/(sp*sp) * Color::Kgph * cdbl(Color::NC) * Color::CF;\
    fPtr4dbl fVV = 0;\
    fPtr4dbl fVA = 0;\
    fPtr4dbl fAA = 0;\
    this->getBQED(fVV, fVA, fAA);
#define cg0_VV n * fVV(m2,-Q2,sp,t1)
#define cg0_VA n * fVA(m2,-Q2,sp,t1)
#define cg0_AA n * fAA(m2,-Q2,sp,t1)
implementPartonicCoeff(Inclusive,cg0)

// setup dq1
void Inclusive::IntKer::getIntA2(fPtr5dbl &fVV, fPtr5dbl &fVA, fPtr5dbl &fAA) const { getME(Inclusive::ME::IntA2) }
#define init_dq1 _sp\
    cdbl n = m2/(4.*M_PI*sp*sp) * Color::Kqph * cdbl(Color::NC) * Color::CF;\
    fPtr5dbl fVV = 0;fPtr5dbl fVA = 0;fPtr5dbl fAA = 0;\
    this->getIntA2(fVV, fVA, fAA);
cdbl Inclusive::IntKer::dq1_cur() const {
    init_dq1
    if (Mode_dq1_VV == this->mode) return n * fVV(m2,-Q2,sp,t1,s4);
    if (Mode_dq1_VA == this->mode) return n * fVA(m2,-Q2,sp,t1,s4);
    if (Mode_dq1_AA == this->mode) return n * fAA(m2,-Q2,sp,t1,s4);
    return 0.;
}

// setup oq1
void Inclusive::IntKer::getIntA3(fPtr5dbl &fVV, fPtr5dbl &fVA, fPtr5dbl &fAA) const { getME(Inclusive::ME::IntA3) }
#define init_oq1 _sp\
    cdbl n = m2/(4.*M_PI*sp*sp) * Color::Kqph * cdbl(Color::NC) * Color::CF;\
    fPtr5dbl fVV = 0;fPtr5dbl fVA = 0;fPtr5dbl fAA = 0;\
    this->getIntA3(fVV, fVA, fAA);

// setup cq1
void Inclusive::IntKer::getIntA1(fPtr5dbl &fVV, fPtr5dbl &fVA, fPtr5dbl &fAA) const { getME(Inclusive::ME::IntA1) }
#define init_cq1 _sp\
    cdbl n = m2/(4.*M_PI*sp*sp) * Color::Kqph * cdbl(Color::NC) * Color::CF;\
    fPtr5dbl fVV = 0;fPtr5dbl fVA = 0;fPtr5dbl fAA = 0;\
    this->getIntA1(fVV, fVA, fAA);
#define cq1_VV n * fVV(m2,-Q2,sp,t1,s4)
#define cq1_VA n * fVA(m2,-Q2,sp,t1,s4)
#define cq1_AA n * fAA(m2,-Q2,sp,t1,s4)
implementPartonicCoeff(Inclusive,cq1)

// setup cqBarF1
#define init_cqBarF1 _sp\
    fPtr4dbl fVV = 0;fPtr4dbl fVA = 0;fPtr4dbl fAA = 0;\
    this->getBQED(fVV, fVA, fAA);\
    cdbl u1 = s4 - sp - t1;\
    cdbl x1 = - u1/(sp + t1);\
    fPtr1dbl Pgq0 = this->getPgq0();\
    cdbl n = m2/(4.*M_PI*sp*sp) * Color::Kqph * cdbl(Color::NC) * Pgq0(x1)/u1;
#define cqBarF1_VV n * fVV(m2,-Q2,x1*sp,x1*t1)
#define cqBarF1_VA n * fVA(m2,-Q2,x1*sp,x1*t1)
#define cqBarF1_AA n * fAA(m2,-Q2,x1*sp,x1*t1)
implementPartonicCoeff(Inclusive,cqBarF1)

// soft helper coefficients
cdbl Inclusive::IntKer::getDeltaCoeffA0() const {
    return 1./(this->s4max - this->Delta);
}
cdbl Inclusive::IntKer::getDeltaCoeffA1() const {
    return log(this->s4max/this->m2)/(this->s4max - this->Delta) - 1./this->s4;
}
cdbl Inclusive::IntKer::getDeltaCoeffA2() const {
    return pow(log(this->s4max/this->m2),2)/(this->s4max - this->Delta) - 2.*log(this->s4/m2)/this->s4;
}

cdbl Inclusive::IntKer::cgBarF1(cint parts) const {
    _sp
    fPtr4dbl fVV = 0;fPtr4dbl fVA = 0;fPtr4dbl fAA = 0;
    this->getBQED(fVV, fVA, fAA);
    cdbl n = m2/(4.*M_PI*sp*sp) * Color::Kgph * cdbl(Color::NC) * Color::CF;
    dbl r = 0.;
    if (parts >= 0) { // H
        cdbl u1 = s4 - sp - t1;
        cdbl x1 = - u1/(sp + t1);
        const fPtr1dbl PggH0 = this->getPggH0();
        cdbl AP = -2./u1*PggH0(x1)*(-1.);
        combineModesAndCurs(meH,cgBarF1,fVV(m2,-Q2,x1*sp,x1*t1),fVA(m2,-Q2,x1*sp,x1*t1),fAA(m2,-Q2,x1*sp,x1*t1));
        r += n*AP*meH;
    }
    if (parts <= 0) { // S
        cdbl A0 = this->getDeltaCoeffA0();
        cdbl A1 = this->getDeltaCoeffA1();
        cdbl SoftDelta1 = -4.*Color::CA;
        cdbl u1S = -sp - t1;
        cdbl SoftDelta0 = 4.*Color::CA * log(-u1S/this->m2) - this->beta0lf();
        cdbl AS = A0*SoftDelta0 + A1*SoftDelta1;
        combineModesAndCurs(meS,cgBarF1,fVV(m2,-Q2,sp,t1),fVA(m2,-Q2,sp,t1),fAA(m2,-Q2,sp,t1));
        r += n*AS*meS;
    }
    return r;
}

cdbl Inclusive::IntKer::cgBarR1(cint parts) const {
    _sp
    cdbl n = 4.*M_PI*m2/(sp*sp) * Color::Kgph * cdbl(Color::NC) * Color::CF;
    fPtr4dbl fVV = 0;fPtr4dbl fVA = 0;fPtr4dbl fAA = 0;
    this->getBQED(fVV, fVA, fAA);
    cdbl b = this->beta0lf()/(16.*M_PI*M_PI);
    combineModesAndCurs(me,cgBarR1,fVV(m2,-Q2,sp,t1),fVA(m2,-Q2,sp,t1),fAA(m2,-Q2,sp,t1));
    cdbl r = n * b * me;
    if (0 == parts) // partonic
        return r;
    else if (parts < 0) { // S
        cdbl A0 = this->getDeltaCoeffA0();
        return A0*r;
    } else // H
        return 0.;
}

// setup cg1
void Inclusive::IntKer::getIntROK(fPtr5dbl &fVV, fPtr5dbl &fVA, fPtr5dbl &fAA)  const { getME(Inclusive::ME::IntROK) }
void Inclusive::IntKer::getIntRQED(fPtr5dbl &fVV, fPtr5dbl &fVA, fPtr5dbl &fAA) const { getME(Inclusive::ME::IntRQED) }
void Inclusive::IntKer::getSVOK(fPtr4dbl &fVV, fPtr4dbl &fVA, fPtr4dbl &fAA)    const { getME(Inclusive::ME::SVOK) }
void Inclusive::IntKer::getSVQED(fPtr4dbl &fVV, fPtr4dbl &fVA, fPtr4dbl &fAA)   const { getME(Inclusive::ME::SVQED) }
cdbl Inclusive::IntKer::cg1(cint parts) const {
    _sp
    cdbl n = m2/(4.*M_PI*sp*sp) * Color::Kgph * cdbl(Color::NC) * Color::CF;
    
    dbl r = 0.;
    if (parts >= 0) { // H
        {
            fPtr5dbl fVV = 0;fPtr5dbl fVA = 0;fPtr5dbl fAA = 0;
            this->getIntROK(fVV, fVA, fAA);
            combineModesAndCurs(meROK,cg1,fVV(m2,-Q2,sp,t1,s4),fVA(m2,-Q2,sp,t1,s4),fAA(m2,-Q2,sp,t1,s4));
            r += n * Color::CA * meROK;
        } {
            fPtr5dbl fVV = 0;fPtr5dbl fVA = 0;fPtr5dbl fAA = 0;
            this->getIntRQED(fVV, fVA, fAA);
            combineModesAndCurs(meRQED,cg1,fVV(m2,-Q2,sp,t1,s4),fVA(m2,-Q2,sp,t1,s4),fAA(m2,-Q2,sp,t1,s4));
            r += n * 2.*Color::CF * meRQED;
        }
        /*if(!isfinite(IntRQED_VV(m2,-Q2,sp,t1,s4)))  printf(">%.10e\t%.10e IntRQED_VV\n",t1,s4);\
        if(!isfinite(IntROK_VV(m2,-Q2,sp,t1,s4))) printf(">%.10e\t%.10e IntROK_VV\n",t1,s4);*/\
    }
    if (parts <= 0) { // S
        {
            cdbl A1 = this->getDeltaCoeffA1();
            cdbl A2 = this->getDeltaCoeffA2();
            cdbl SoftDelta2 = 4.*Color::CA;
            cdbl beta = this->beta();
            cdbl chi = (1. - beta)/(1. + beta);
            cdbl u1S = -sp - t1;
            cdbl SoftDelta1 = 4.*(Color::CA*(log(-t1/m2)-log(-u1S/m2))-(Color::CA - 2.*Color::CF)*(2.*m2-s)/(s*beta)*log(chi) - 2.*Color::CF);
            
            fPtr4dbl fVV = 0;fPtr4dbl fVA = 0;fPtr4dbl fAA = 0;
            this->getBQED(fVV, fVA, fAA);
            combineModesAndCurs(meB,cg1,fVV(m2,-Q2,sp,t1),fVA(m2,-Q2,sp,t1),fAA(m2,-Q2,sp,t1));
            r += n * (A2*SoftDelta2 + A1*SoftDelta1) * meB;
        }
        cdbl A0 = this->getDeltaCoeffA0();
        {
            fPtr4dbl fVV = 0;fPtr4dbl fVA = 0;fPtr4dbl fAA = 0;\
            this->getSVOK(fVV, fVA, fAA);\
            combineModesAndCurs(meSVOK,cg1,fVV(m2,-Q2,sp,t1),fVA(m2,-Q2,sp,t1),fAA(m2,-Q2,sp,t1));
            r += n * A0 * Color::CA * meSVOK;
        } {
            fPtr4dbl fVV = 0;fPtr4dbl fVA = 0;fPtr4dbl fAA = 0;\
            this->getSVQED(fVV, fVA, fAA);\
            combineModesAndCurs(meSVQED,cg1,fVV(m2,-Q2,sp,t1),fVA(m2,-Q2,sp,t1),fAA(m2,-Q2,sp,t1));
            r += n * A0 * 2.*Color::CF * meSVQED;
        }
    }
    return r;
}

cdbl Inclusive::IntKer::Fg0() const {
    cdbl alphaS = this->getAlphaS(this->HAQTransverseMomentum);
    cdbl curMuF2 = this->getScale(this->muF2,this->HAQTransverseMomentum);
    cdbl g = this->pdf->xfxQ2(21,this->xi,curMuF2);
    cdbl nLO = alphaS/(4.*M_PI*M_PI) * (this->Q2/this->m2);
    cdbl r = nLO * g * this->cg0();
    return r;
}

cdbl Inclusive::IntKer::Fg1(cint parts) const {
    cdbl alphaS = this->getAlphaS(this->HAQTransverseMomentum);
    cdbl curMuF2 = this->getScale(this->muF2,this->HAQTransverseMomentum);
    cdbl curMuR2 = this->getScale(this->muR2,this->HAQTransverseMomentum);
    cdbl g = this->pdf->xfxQ2(21,this->xi,curMuF2);
    cdbl nNLO = alphaS*alphaS/(M_PI) * (this->Q2/this->m2);
    return nNLO * g * (this->cg1(parts) + log(curMuF2/m2)*this->cgBarF1(parts) + log(curMuR2/m2)*this->cgBarR1(parts));
}

cdbl Inclusive::IntKer::Fq1() const {
    // compute matrix elements for dq1+oq1
    // as they are combined with the properties of the light quarks,
    // we have to keep the raw-current elements
    dbl dq1_VV = 0.,dq1_VA = 0.,dq1_AA = 0.;
    {
        init_dq1
        if (isParityConservingProj(this->proj)) {
            dq1_VV = n*fVV(this->m2,-this->Q2,sp,t1,s4);
            dq1_AA = n*fAA(this->m2,-this->Q2,sp,t1,s4);
        } else {
            dq1_VA = n*fVA(this->m2,-this->Q2,sp,t1,s4);
        }
    }
    dbl oq1_VV = 0.,oq1_VA = 0.,oq1_AA = 0.;
    if (this->mode > Mode_F) {
        init_oq1
        if (isParityConservingProj(this->proj)) {
            oq1_VV = n*fVV(this->m2,-this->Q2,sp,t1,s4);
            oq1_AA = n*fAA(this->m2,-this->Q2,sp,t1,s4);
        } else {
            oq1_VA = n*fVA(this->m2,-this->Q2,sp,t1,s4);
        }
    }
    // common stuff
    cdbl alphaS = this->getAlphaS(this->HAQTransverseMomentum);
    cdbl curMuF2 = this->getScale(this->muF2,this->HAQTransverseMomentum);
    cdbl nNLO = alphaS*alphaS/(M_PI) * (this->Q2/this->m2);
    cdbl eH = getElectricCharge(this->nlf+1);
    cdbl gVQ = this->getVectorialCoupling(this->nlf+1);
    cdbl gAQ = this->getAxialCoupling(this->nlf+1);
    // join elems
    cdbl cq1j = this->cq1() + log(curMuF2/this->m2) * this->cqBarF1();
    dbl fqs = 0.;
    for (uint q = 1; q < this->nlf + 1; ++q) {
        dbl dq1 = 0., oq1 = 0.;
        cdbl eL = getElectricCharge(q);
        cdbl gVq = this->getVectorialCoupling(q);
        cdbl gAq = this->getAxialCoupling(q);
        if (isParityConservingProj(this->proj)) {
            if (this->flags.usePhoton) {
                dq1 += eL*eL * dq1_VV; 
                oq1 += eH*eL * oq1_VV;
            }
            if (this->flags.usePhotonZ) {
                dq1 -= this->getNormPhZ() *  eL*gVq           * dq1_VV;
                oq1 -= this->getNormPhZ() * (eH*gVq + eL*gVQ) * oq1_VV;
            }
            if (this->flags.useZ) {
                dq1 += this->getNormZ()*(gVq*gVq*dq1_VV + gAq*gAq*dq1_AA);
                oq1 += this->getNormZ()*(gVQ*gVq*oq1_VV + gAQ*gAq*oq1_AA);
            }
        } else {
            if (this->flags.usePhotonZ) {
                dq1 -= this->getNormPhZ() *  eL*gAq           * dq1_VA;
                oq1 -= this->getNormPhZ() * (eH*gAq + eL*gAQ) * oq1_VA;
            }
            if (this->flags.useZ) {
                dq1 += this->getNormZ() * 2.*gVq*gAq          * dq1_VA;
                oq1 += this->getNormZ() * (gVQ*gAq + gVq*gAQ) * oq1_VA;
            }
        }
        // combine cq1+cqBarF+dq1+oq1
        fqs += (this->pdf->xfxQ2((int)q,this->xi,curMuF2) + this->pdf->xfxQ2(-((int)q),this->xi,curMuF2))*(cq1j + dq1 + oq1);
    }
    return nNLO*fqs;
}

cdbl Inclusive::IntKer::runPartonic(cdbl a1, cdbl a2) {
    // 1D integrations
    this->setT1(a1);
    // cg0
    if (Mode_cg0_VV == this->mode || Mode_cg0_VA == this->mode || Mode_cg0_AA == this->mode)
        return this->jac_t1*this->cg0_cur();
    // cgBarR1
    if (Mode_cgBarR1_VV == this->mode || Mode_cgBarR1_VA == this->mode || Mode_cgBarR1_AA == this->mode)
        return this->jac_t1*this->cgBarR1(0);
    // setS4 in NLOq
    _sp
    //cdbl beta = this->beta();
    this->s4max = this->s + this->s*this->t1/sp + sp*this->m2/this->t1; //(this->s*((sp*(1. - beta))/2. + this->t1)*((sp*(1. + beta))/2. + this->t1))/(sp*this->t1);
    if (s4max <= 0.)
        return 0.;
    this->s4 = a2*s4max;
    // dq1
    if (Mode_dq1_VV == this->mode || Mode_dq1_VA == this->mode || Mode_dq1_AA == this->mode)
        return this->jac_t1*this->s4max*this->dq1_cur();
    // cq1
    if (Mode_cq1_VV == this->mode || Mode_cq1_VA == this->mode || Mode_cq1_AA == this->mode)
        return this->jac_t1*this->s4max*this->cq1_cur();
    // cqBarF1
    if (Mode_cqBarF1_VV == this->mode || Mode_cqBarF1_VA == this->mode || Mode_cqBarF1_AA == this->mode)
        return this->jac_t1*this->s4max*this->cqBarF1_cur();
    // setS4 in NLOg
    cdbl V_s4 = s4max - this->Delta;
    if (V_s4 <= 0.)
        return 0.;
    this->s4 = this->Delta + a2*V_s4;
    if(this->s4 <= this->Delta || this->s4 > s4max)
        return 0.;
    // cg1
    if (Mode_cg1_VV == this->mode || Mode_cg1_VA == this->mode || Mode_cg1_AA == this->mode)
        return this->jac_t1*V_s4*this->cg1(0);
    // cgBarF1
    if (Mode_cgBarF1_VV == this->mode || Mode_cgBarF1_VA == this->mode || Mode_cgBarF1_AA == this->mode)
        return this->jac_t1*V_s4*this->cgBarF1(0);
    // cgBar1
    if (Mode_cgBar1_VV == this->mode || Mode_cgBar1_VA == this->mode || Mode_cgBar1_AA == this->mode) {
        cuint old = this->mode;
        uint F = 0, R = 0;
        if (Mode_cgBar1_VV == this->mode) { F = this->Mode_cgBarF1_VV; R = this->Mode_cgBarR1_VV; }
        else if (Mode_cgBar1_VA == this->mode) { F = this->Mode_cgBarF1_VA; R = this->Mode_cgBarR1_VA; }
        else if (Mode_cgBar1_AA == this->mode) { F = this->Mode_cgBarF1_AA; R = this->Mode_cgBarR1_AA; }
        dbl r = 0.;
        this->mode = F;
        r += this->jac_t1*V_s4*this->cgBarF1(0);
        this->mode = R;
        r += this->jac_t1*this->cgBarR1(0);
        this->mode = old;
        return r;
    }
    return 0.;
}

cdbl Inclusive::IntKer::runHadronic(cdbl a1, cdbl a2, cdbl a3) {
    dbl jac = 1.;
    if (Mode_F == this->mode) {
        this->HAQRapidity = dblNaN;
        this->HAQTransverseMomentum = dblNaN;
        this->setHAQRapidity(a1);
        this->setHAQTransverseMomentum(a2);
        jac = this->V_HAQTransverseMass2*this->V_HAQRapidity;
    } else if (Mode_dF_dHAQTransverseMomentum == this->mode) {
        this->HAQRapidity = dblNaN;
        this->setHAQRapidity(a1);
        cdbl jac2 = 2. * this->HAQTransverseMomentum;
        jac = this->V_HAQRapidity * jac2;
    } else if (Mode_dF_dHAQRapidity == this->mode) {
        this->HAQTransverseMomentum = dblNaN;
        this->setHAQTransverseMomentum(a1);
        jac = this->V_HAQTransverseMass2;
    } else if (Mode_dF_dHAQFeynmanX == this->mode) {
        cdbl Sh = this->getHadronicS();
        cdbl betah = sqrt(1. - 4.*this->m2/Sh);
        cdbl mt2Max = Sh/4.*(1. - betah*betah * this->HAQFeynmanX*this->HAQFeynmanX);
        cdbl jac_mt2 = mt2Max - this->m2;
        cdbl mt2 = this->m2 + jac_mt2*a1;
        this->HAQTransverseMomentum = sqrt(mt2 - this->m2);
        cdbl mt = sqrt(mt2);
        cdbl pLmax = sqrt(Sh)/2. * betah;
        cdbl chi = sqrt(this->HAQFeynmanX*this->HAQFeynmanX + (mt/pLmax)*(mt/pLmax));
        cdbl ey = pLmax/mt*(this->HAQFeynmanX + chi);
        this->HAQRapidity = log(ey);
        jac = jac_mt2/chi;
    }
    cdbl Shp = this->getHadronicS() + this->Q2;
    cdbl T1 = this->getHadronicT1();
    jac *= Shp/(Shp + T1);
    dbl r = 0.;
    // LO
    if (this->flags.useLeadingOrder && this->flags.useGluonicChannel) {
        this->s4 = 0;
        this->setPartonicVars();
        r += jac * this->Fg0();
    }
    // NLO
    if (this->flags.useNextToLeadingOrder) {
        // setS4
        cdbl U1 = this->getHadronicU1();
        this->s4max = Shp + T1 + U1;
        cdbl as4 = (Mode_F == this->mode) ? a3 : a2;
        if (this->flags.useQuarkChannel) {
            if (s4max > 0.) {
                this->s4 = s4max*as4;
                this->setPartonicVars();
                r += jac*s4max * this->Fq1();
            }
        }
        if (this->flags.useGluonicChannel) {
            cdbl V_s4 = s4max - this->Delta;
            if (V_s4 > 0.) {
                // split H and S, because xi=xi(s4)
                // S
                this->s4 = 0.;
                this->setPartonicVars();
                this->s4 = this->Delta + V_s4*as4;
                r += jac*V_s4 * this->Fg1(-1);
                // H
                this->setPartonicVars();
                r += jac*V_s4 * this->Fg1(1);
            }
        }
    }
    return r;
}

cdbl Inclusive::IntKer::runLeptonic(cdbl a1, cdbl a2, cdbl a3, cdbl a4, cdbl a5) {
    // setYBj
    cdbl yBjMax = 1.;
    cdbl yBjMin = 4.*this->m2 / this->Sl;
    cdbl V_yBj = yBjMax - yBjMin;
    if (V_yBj <= 0.)
        return 0.;
    cdbl yBj = yBjMin + V_yBj*a1;
    cdbl OmYBj = 1. - yBj;
    // set Q2
    dbl curQ2min = this->Q2min;
    if (!isfinite(this->Q2min)) {
        curQ2min = this->q2minHVQDIS * yBj*yBj/OmYBj/OmYBj;
    }
    cdbl Q2max = yBj * this->Sl - 4.*this->m2;
    cdbl V_Q2 = Q2max - curQ2min;
    if (V_Q2 <= 0.)
        return 0.;
    this->Q2 = curQ2min + V_Q2*a2;
    this->xBj = this->Q2/yBj/this->Sl;
    if(this->xBj >= this->getZMax())
        return 0.;
    // normalisation
    cdbl n = 2.*M_PI*this->alphaEM*this->alphaEM/yBj/this->Q2/this->Q2;
    cdbl polN = this->polarizeBeams ? 2. : 1.;
    // collect structure functions
    dbl F2_mg4 = 0.;
    dbl FL_mgL = 0.;
    dbl xF3_x2g1 = 0.;
    this->mode = Mode_F;
    if (this->polarizeBeams) {
        this->proj = x2g1;
        xF3_x2g1 = this->runHadronic(a3,a4,a5);
        this->proj = g4;
        F2_mg4 = -this->runHadronic(a3,a4,a5);
        this->proj = gL;
        FL_mgL = -this->runHadronic(a3,a4,a5);
    } else {
        this->proj = F2;
        F2_mg4 = this->runHadronic(a3,a4,a5);
        this->proj = FL;
        FL_mgL = this->runHadronic(a3,a4,a5);
        this->proj = xF3;
        xF3_x2g1 = this->runHadronic(a3,a4,a5);
    }
    this->mode = Mode_sigma;
    // combine
    cdbl Yp = 1. + OmYBj*OmYBj;
    cdbl Ym = 1. - OmYBj*OmYBj;
    cdbl YmSign = this->incomingLeptonHasNegativeCharge ? 1. : -1. ;
    cdbl sig = n*polN*(Yp*F2_mg4 + YmSign*Ym*xF3_x2g1 - yBj*yBj*FL_mgL);
    return V_yBj*V_Q2*sig;
}

cdbl Inclusive::IntKer::operator()(cdbl a1, cdbl a2, cdbl a3, cdbl a4, cdbl a5) {      
    // void mode
    if (0 == this->mode) return 0.;
    // partonic mode?
    if (this->mode < Mode_F) {
        cdbl r = this->runPartonic(a1, a2);
        //if(!isfinite(r)) printf("<%.10e\t%.10e ()\n",a1,a2);
        return isfinite(r) ? r : 0.;
    }
    // hadronic mode?
    if (this->mode < Mode_sigma) {
        cdbl r = this->runHadronic(a1, a2, a3);
        return isfinite(r) ? r : 0.;
    }
    // leptonic mode = Mode_sigma
    cdbl r = this->runLeptonic(a1, a2, a3, a4, a5);
    return isfinite(r) ? r : 0.;
}

void Inclusive::IntKer::operator()(const double x[], const int k[], const double& weight, const double aux[], double f[]) {
    // proide simple redirect as weight isn't needed
    cdbl a0 = this->dim >= 1 ? x[0] : 0.;
    cdbl a1 = this->dim >= 2 ? x[1] : 0.;
    cdbl a2 = this->dim >= 3 ? x[2] : 0.;
    cdbl a3 = this->dim >= 4 ? x[3] : 0.;
    cdbl a4 = this->dim >= 5 ? x[4] : 0.;
    f[0] = this->operator()(a0,a1,a2,a3,a4);
}