#include "FullyDiff/IntKer.h"

#include "Common/Color.hpp"
#include "Common/ME/BQED.h"
#include "FullyDiff/KinematicVars.hpp"

#include "ME/A1.h"
#include "ME/A1Counter.h"
#include "ME/A2.h"
#include "ME/A3.h"

#include "ME/R.h"
#include "ME/RCounterX.h"
#include "ME/RCounterY.h"
#include "ME/RCounterXY.h"
#include "ME/SV.h"

FullyDiff::IntKer::IntKer() : AbstractIntKer() {
};

FullyDiff::IntKer::~IntKer() {
    // delete all histograms
    for(FullyDiff::histMapT::const_iterator it = this->histMap.cbegin(); it != this->histMap.cend(); ++it) {
        delete (it->second);
    }
};

void FullyDiff::IntKer::setTheta1(cdbl a) {
    this->Theta1 = this->V_Theta1*a;
}

void FullyDiff::IntKer::setX(cdbl a) {
    // setSpRaw
    cdbl eta = this->s/(4.*m2) - 1.;
    cdbl cut = this->deltax*pow(eta,2./3.);
    cdbl xmax = 1. - cut;
    // reached numeric limit?
    if (0. != this->deltax && 1. == xmax)
        throw domain_error((boost::format("x_max = 1-deltax*eta^(2/3) = 1. - %e = %e has to be smaller than 1!")%cut%xmax).str());
    // make hard cut instead of throwing an error
    _sp
    cdbl rhoStar = min((4.*m2 + this->Q2)/sp,xmax);
    this->V_xE = xmax - rhoStar;
    this->rhoTilde = min(1. - this->xTilde*(1. - rhoStar),xmax);
    this->V_xC = xmax - this->rhoTilde;
    this->xE = rhoStar  + this->V_xE * a;
    this->xC = this->rhoTilde + this->V_xC * a;
}

void FullyDiff::IntKer::setTheta2(cdbl a) {
    this->Theta2 = this->V_Theta2*a;
}

void FullyDiff::IntKer::setY(cdbl a) {
    cdbl ymin = -1.+this->deltay;
    this->V_yE = 1. - ymin;
    this->V_yC = (-1. + this->omega) - ymin;
    this->yE = ymin + this->V_yE*a;
    this->yC = ymin + this->V_yC*a;
}

void FullyDiff::IntKer::setZ(cdbl a) {
    cdbl zMax = this->getZMax();
    this->V_z = zMax - this->xBj;
    this->z = this->xBj + this->V_z*a;
    cdbl sp = this->Q2/this->z;
    this->s = sp - this->Q2;
}

cdbl FullyDiff::IntKer::cg0() const {
    _sp
    cdbl beta = this->beta();
    cdbl t1 = -.5*sp*(1. - beta*cos(Theta1));
    cdbl n = 8.*M_PI * Color::Kgph*cdbl(Color::NC)*Color::CF * 1./(4.*sp) * m2  * beta*sin(Theta1) * this->V_Theta1;
    
    fPtr4dbl fVV = 0;fPtr4dbl fVA = 0;fPtr4dbl fAA = 0;
    this->getBQED(fVV, fVA, fAA);
    combineModesAndCurs(me,cg0,fVV(m2,-Q2,sp,t1),fVA(m2,-Q2,sp,t1),fAA(m2,-Q2,sp,t1));
    return n*me;
}

// setup cg1
void FullyDiff::IntKer::getR(fPtr7dbl &fVV, fPtr7dbl &fVA, fPtr7dbl &fAA) const { getME(FullyDiff::ME::R) }
void FullyDiff::IntKer::getRCounterX(fPtr6dbl &fVV, fPtr6dbl &fVA, fPtr6dbl &fAA) const { getME(FullyDiff::ME::RCounterX) }
void FullyDiff::IntKer::getRCounterY(fPtr6dbl &fVV, fPtr6dbl &fVA, fPtr6dbl &fAA) const { getME(FullyDiff::ME::RCounterY) }
void FullyDiff::IntKer::getRCounterXY(fPtr4dbl &fVV, fPtr4dbl &fVA, fPtr4dbl &fAA) const { getME(FullyDiff::ME::RCounterXY) }
void FullyDiff::IntKer::getSV(fPtr5dbl &fVV, fPtr5dbl &fVA, fPtr5dbl &fAA) const { getME(FullyDiff::ME::SV) }
const FullyDiff::PhasespaceValues FullyDiff::IntKer::cg1() const {
    _sp
    // norm to cg1
    cdbl ncg1 = (m2/(4.*M_PI));
    PhasespaceValues r;
    
    cdbl beta = sqrt(1. - 4.*m2/s);
    cdbl t1sc = -.5*sp*(1. - beta*cos(Theta1));
    
    { // S+V contributions
        cdbl f = Color::Kgph*cdbl(Color::NC)*Color::CF * 1./(4.*sp);
        cdbl betaTilde = sqrt(1. - this->rhoTilde);
        cdbl g = ncg1 * V_Theta1 * f * beta*sin(Theta1)/(16.*M_PI);
        
        /// @todo speed up SV? grid or parametrization? "complicated" part only depends on m2,Q2,s,cos(theta1)
        fPtr5dbl fVV = 0;fPtr5dbl fVA = 0;fPtr5dbl fAA = 0;
        this->getSV(fVV, fVA, fAA);
        combineModesAndCurs(meSV,cg1,fVV(m2,-Q2,sp,t1sc,betaTilde),fVA(m2,-Q2,sp,t1sc,betaTilde),fAA(m2,-Q2,sp,t1sc,betaTilde));
        
        r.xCyE += g * meSV;
    } { // collinear contributions
        cdbl s5E = -Q2 + xE*sp;
        cdbl beta5E = sqrt(1. - 4.*m2/s5E);
        cdbl t1c  = -.5*sp*(1. - beta5E*cos(Theta1));
        cdbl f = Color::Kgph*cdbl(Color::NC)*Color::CF * 1./sp * sin(Theta1);
        cdbl l = log(sp/m2*sp/s*omega/2.);
        const fPtr1dbl PggH0 = this->getPggH0();
        const fPtr1dbl PggH1 = this->getPggH1();
        //     (1-x)P_gg^{H,0} -> 2CA for x->1 for all projections
        // and (1-x)P_gg^{H,1} -> 0   for x->1 for all projections
        cdbl gE = ncg1 * V_xE*V_Theta1 * f/xE*beta5E * (PggH0(xE)   *(           l + 2.*log(1.-xE)        ) + 2.*PggH1(xE));
        cdbl gC = ncg1 * V_xC*V_Theta1 * f   *beta   * (2.*Color::CA*(1./(1.-xC)*l + 2.*log(1.-xC)/(1.-xC))               );
        
        fPtr4dbl fVV = 0;fPtr4dbl fVA = 0;fPtr4dbl fAA = 0;
        this->getBQED(fVV, fVA, fAA);
        combineModesAndCurs(meE,cg1,fVV(m2,-Q2,xE*sp,xE*t1c),fVA(m2,-Q2,xE*sp,xE*t1c),fAA(m2,-Q2,xE*sp,xE*t1c));
        combineModesAndCurs(meC,cg1,fVV(m2,-Q2,sp,t1sc),fVA(m2,-Q2,sp,t1sc),fAA(m2,-Q2,sp,t1sc));
        
        r.xEyC += gE*meE;
        r.xCyC -= gC*meC;
    } { // hard contributions
        cdbl s5E = -Q2 + xE*sp;
        cdbl beta5E = sqrt(1. - 4.*m2/s5E);
        cdbl f = Color::Kgph*cdbl(Color::NC)*Color::CF * s/(M_PI*pow(sp,3))*sin(Theta1);
        cdbl gEE = ncg1 * V_xE*V_yE*V_Theta1*V_Theta2 * f*beta5E/(1.-yE) * 1./(1.-xE)/(1.+yE);
        cdbl gCE = ncg1 * V_xC*V_yE*V_Theta1*V_Theta2 * f*beta  /(1.-yE) * 1./(1.-xC)/(1.+yE);
        cdbl gEC = ncg1 * V_xE*V_yC*V_Theta1*V_Theta2 * f*beta5E/2.      * 1./(1.-xE)/(1.+yC);
        cdbl gCC = ncg1 * V_xC*V_yC*V_Theta1*V_Theta2 * f*beta  /2.      * 1./(1.-xC)/(1.+yC);
        
        {
            fPtr7dbl fVV = 0;fPtr7dbl fVA = 0;fPtr7dbl fAA = 0;
            this->getR(fVV, fVA, fAA);
            combineModesAndCurs(meEE,cg1,fVV(m2,-Q2,sp,xE,yE,Theta1,Theta2),fVA(m2,-Q2,sp,xE,yE,Theta1,Theta2),fAA(m2,-Q2,sp,xE,yE,Theta1,Theta2));
            r.xEyE += gEE * meEE;
        } {
            fPtr6dbl fVV = 0;fPtr6dbl fVA = 0;fPtr6dbl fAA = 0;
            this->getRCounterX(fVV, fVA, fAA);
            combineModesAndCurs(meCE,cg1,fVV(m2,-Q2,sp,yE,Theta1,Theta2),fVA(m2,-Q2,sp,yE,Theta1,Theta2),fAA(m2,-Q2,sp,yE,Theta1,Theta2));
            r.xCyE -= gCE * meCE;
        } {
            fPtr6dbl fVV = 0;fPtr6dbl fVA = 0;fPtr6dbl fAA = 0;
            this->getRCounterY(fVV, fVA, fAA);
            combineModesAndCurs(meEC,cg1,fVV(m2,-Q2,sp,xE,Theta1,Theta2),fVA(m2,-Q2,sp,xE,Theta1,Theta2),fAA(m2,-Q2,sp,xE,Theta1,Theta2));
            r.xEyC -= gEC * meEC;
        } {
            fPtr4dbl fVV = 0;fPtr4dbl fVA = 0;fPtr4dbl fAA = 0;
            this->getRCounterXY(fVV, fVA, fAA);
            combineModesAndCurs(meCC,cg1,fVV(m2,-Q2,sp,t1sc),fVA(m2,-Q2,sp,t1sc),fAA(m2,-Q2,sp,t1sc));
            r.xCyC += gCC * meCC;
        }
    }
    
    return r;
}

const FullyDiff::PhasespaceValues FullyDiff::IntKer::cgBarF1() const {
    // collinear contributions
    _sp
    cdbl beta = this->beta();
    // norm to cg1
    cdbl ncg1 = (m2/(4.*M_PI));
    PhasespaceValues r;
    
    cdbl t1 = -.5*sp*(1. - beta*cos(Theta1));
    fPtr4dbl fVV = 0;fPtr4dbl fVA = 0;fPtr4dbl fAA = 0;
    this->getBQED(fVV, fVA, fAA);
    combineModesAndCurs(meC,cgBarF1,fVV(m2,-Q2,sp,t1),fVA(m2,-Q2,sp,t1),fAA(m2,-Q2,sp,t1));
    
    { // S+V contributions
        cdbl f = 2. * Color::Kgph*cdbl(Color::NC)*Color::CF * 1./(4.*sp);
        cdbl b0 = this->beta0lf();
        // PggS0 = b0/2 + 4CA ln(betaTilde) for all projections
        cdbl g = ncg1 * V_Theta1 * f * (b0 + 4.*Color::CA*log(1.-rhoTilde)) * beta*sin(Theta1);
        r.xCyE -= g * meC;
    } { // collinear contributions
        cdbl s5E = -Q2 + xE*sp;
        cdbl beta5E = sqrt(1. - 4.*m2/s5E);
        cdbl t1c = -.5*sp*(1. - beta5E*cos(Theta1));
        cdbl f = Color::Kgph*cdbl(Color::NC)*Color::CF * 1./sp * sin(Theta1);
        cdbl l = -1.;
        const fPtr1dbl PggH0 = this->getPggH0();
        // (1-x)P_gg^{H,0} -> 2CA for x->1 for all projections
        cdbl gE = ncg1 * V_xE*V_Theta1 * f*beta5E/xE*(PggH0(xE)*l);
        cdbl gC = ncg1 * V_xC*V_Theta1 * f*beta     *(2.*Color::CA*(1./(1.-xC)*l));
        
        combineModesAndCurs(meE,cgBarF1,fVV(m2,-Q2,xE*sp,xE*t1c),fVA(m2,-Q2,xE*sp,xE*t1c),fAA(m2,-Q2,xE*sp,xE*t1c));
        
        r.xEyC += gE*meE;
        r.xCyC -= gC*meC;
    }
    return r;
}

const FullyDiff::PhasespaceValues FullyDiff::IntKer::cgBarR1() const {
    PhasespaceValues r;
    _sp
    cdbl beta = this->beta();
    cdbl t1 = -.5*sp*(1. - beta*cos(Theta1));
    cdbl n = 8.*M_PI * Color::Kgph*cdbl(Color::NC)*Color::CF * 1./(4.*sp) * m2  * beta*sin(Theta1) * this->V_Theta1;
    
    fPtr4dbl fVV = 0;fPtr4dbl fVA = 0;fPtr4dbl fAA = 0;
    this->getBQED(fVV, fVA, fAA);
    combineModesAndCurs(me,cgBarR1,fVV(m2,-Q2,sp,t1),fVA(m2,-Q2,sp,t1),fAA(m2,-Q2,sp,t1));
    
    cdbl b = this->beta0lf()/(16.*M_PI*M_PI);
    r.xCyE += n*b*me;
    return r;
}

// setup cq1
void FullyDiff::IntKer::getA1(fPtr7dbl &fVV, fPtr7dbl &fVA, fPtr7dbl &fAA) const { getME(FullyDiff::ME::A1) }
void FullyDiff::IntKer::getA1Counter(fPtr6dbl &fVV, fPtr6dbl &fVA, fPtr6dbl &fAA) const { getME(FullyDiff::ME::A1Counter) }
const FullyDiff::PhasespaceValues FullyDiff::IntKer::cq1() const {
    _sp
    PhasespaceValues r;
    cdbl ncq1 = m2/(4.*M_PI);
    { // collinear contributions
        cdbl s5B = -Q2 + sp*xE;
        cdbl beta5B = sqrt(1. - 4.*m2/s5B);
        cdbl t1c = -.5*sp*(1.-beta5B*cos(Theta1));
        cdbl jacB = V_xE*V_Theta1;
        cdbl g = ncq1 * Color::Kqph*cdbl(Color::NC) * 1./(xE*sp)*1./(2.) * beta5B*sin(Theta1);
        cdbl l = log(sp/m2*sp/s*omega/2.*(1.-xE)*(1.-xE));
        const fPtr1dbl Pgq0 = this->getPgq0();
        cdbl vPgq0 = Pgq0(xE);
        const fPtr1dbl Pgq1 = this->getPgq1();
        cdbl vPgq1 = Pgq1(xE);
        
        fPtr4dbl fVV = 0;fPtr4dbl fVA = 0;fPtr4dbl fAA = 0;
        this->getBQED(fVV, fVA, fAA);
        combineModesAndCurs(meB,cq1,fVV(m2,-Q2,xE*sp,xE*t1c),fVA(m2,-Q2,xE*sp,xE*t1c),fAA(m2,-Q2,xE*sp,xE*t1c));
        
        r.xEyC += g*jacB*meB*(vPgq0*l + 2.*vPgq1);
    } { // hard contributions
        const KinematicVars vs(m2,-Q2,sp,xE,yE,Theta1,Theta2);
        cdbl jacE = V_xE*V_yE*V_Theta1*V_Theta2;
        cdbl jacC = V_xE*V_yC*V_Theta1*V_Theta2;
        cdbl f = ncq1 * (-1.)/(4.*M_PI)*1./sp * Color::Kqph*cdbl(Color::NC)*Color::CF * vs.beta5*sin(Theta1);
        
        {
            fPtr7dbl fVV = 0;fPtr7dbl fVA = 0;fPtr7dbl fAA = 0;
            this->getA1(fVV, fVA, fAA);
            combineModesAndCurs(meE,cq1,fVV(m2,-Q2,sp,vs.t1,vs.u1,vs.tp,vs.up),fVA(m2,-Q2,sp,vs.t1,vs.u1,vs.tp,vs.up),fAA(m2,-Q2,sp,vs.t1,vs.u1,vs.tp,vs.up));
            r.xEyE += f * jacE*meE/(1.+yE);
        } {
            fPtr6dbl gVV = 0;fPtr6dbl gVA = 0;fPtr6dbl gAA = 0;
            this->getA1Counter(gVV, gVA, gAA);
            combineModesAndCurs(meC,cq1,gVV(m2,-Q2,sp,xE,Theta1,Theta2),gVA(m2,-Q2,sp,xE,Theta1,Theta2),gAA(m2,-Q2,sp,xE,Theta1,Theta2));
            r.xEyC -= f * jacC*meC/(1.+yC);
        }
    }
    return r;
}

const FullyDiff::PhasespaceValues FullyDiff::IntKer::cqBarF1() const {
    // collinear contributions
    _sp
    PhasespaceValues r;
    cdbl s5B = -Q2 + sp*xE;
    cdbl beta5B = sqrt(1. - 4.*m2/s5B);
    cdbl t1c = -.5*sp*(1. - beta5B*cos(Theta1));
    cdbl jacB = V_xE*V_Theta1;
    cdbl g = Color::Kqph*cdbl(Color::NC) * m2/(xE*sp)*1./(8.*M_PI) * beta5B*sin(Theta1);
    cdbl l = -1.;
    fPtr1dbl Pgq0 = this->getPgq0();
    cdbl vPqg0 = Pgq0(xE);
    
    fPtr4dbl fVV = 0;fPtr4dbl fVA = 0;fPtr4dbl fAA = 0;
    this->getBQED(fVV, fVA, fAA);
    combineModesAndCurs(meB,cqBarF1,fVV(m2,-Q2,xE*sp,xE*t1c),fVA(m2,-Q2,xE*sp,xE*t1c),fAA(m2,-Q2,xE*sp,xE*t1c));
    
    r.xEyC += jacB * g*meB*vPqg0*l;
    return r;
}

// setup dq1
void FullyDiff::IntKer::getA2(fPtr7dbl &fVV, fPtr7dbl &fVA, fPtr7dbl &fAA) const { getME(FullyDiff::ME::A2) }
const FullyDiff::PhasespaceValues FullyDiff::IntKer::dq1() const {
    _sp
    PhasespaceValues r;
    cdbl jac = V_Theta1*V_Theta2*V_xE*V_yE;
    const KinematicVars vs(m2,-Q2,sp,xE,yE,Theta1,Theta2);
    cdbl f = -1./(16.*M_PI*M_PI)*m2/sp * Color::Kqph*cdbl(Color::NC)*Color::CF * vs.beta5/(1.+yE)*sin(Theta1);
    
    fPtr7dbl fVV = 0;fPtr7dbl fVA = 0;fPtr7dbl fAA = 0;
    this->getA2(fVV, fVA, fAA);
    // can't combine modes here, as is proportional to the light quark
    dbl meA2 = 0.;
    if      (Mode_dq1_VV == this->mode) meA2 = fVV(m2,-Q2,sp,vs.t1,vs.u1,vs.tp,vs.up);
    else if (Mode_dq1_VA == this->mode) meA2 = fVA(m2,-Q2,sp,vs.t1,vs.u1,vs.tp,vs.up);
    else if (Mode_dq1_AA == this->mode) meA2 = fAA(m2,-Q2,sp,vs.t1,vs.u1,vs.tp,vs.up);
    
    r.xEyE += jac * f * meA2;
    return r;
}

// setup oq1
void FullyDiff::IntKer::getA3(fPtr7dbl &fVV, fPtr7dbl &fVA, fPtr7dbl &fAA) const { getME(FullyDiff::ME::A3) }
const FullyDiff::PhasespaceValues FullyDiff::IntKer::oq1() const {
    _sp
    PhasespaceValues r;
    cdbl jac = V_Theta1*V_Theta2*V_xE*V_yE;
    const KinematicVars vs(m2,-Q2,sp,xE,yE,Theta1,Theta2);
    cdbl f = -1./(16.*M_PI*M_PI)*m2/sp * Color::Kqph*cdbl(Color::NC)*Color::CF * vs.beta5/(1.+yE)*sin(Theta1);
    
    fPtr7dbl fVV = 0;fPtr7dbl fVA = 0;fPtr7dbl fAA = 0;
    this->getA3(fVV, fVA, fAA);
    // can't combine modes here, as is proportional to the light quark
    dbl meA3 = 0.;
    if      (Mode_oq1_VV == this->mode) meA3 = fVV(m2,-Q2,sp,vs.t1,vs.u1,vs.tp,vs.up);
    else if (Mode_oq1_VA == this->mode) meA3 = fVA(m2,-Q2,sp,vs.t1,vs.u1,vs.tp,vs.up);
    else if (Mode_oq1_AA == this->mode) meA3 = fAA(m2,-Q2,sp,vs.t1,vs.u1,vs.tp,vs.up);
    
    r.xEyE += jac * f * meA3;
    return r;
}

cdbl FullyDiff::IntKer::runPartonic(cdbl a1, cdbl a2, cdbl a3, cdbl a4) {
    // 1D integrations
    this->setTheta1(a1);
    // cg0
    if (Mode_cg0_VV == this->mode || Mode_cg0_VA == this->mode || Mode_cg0_AA == this->mode)
        return this->cg0();
    // cgBarR1
    if (Mode_cgBarR1_VV == this->mode || Mode_cgBarR1_VA == this->mode || Mode_cgBarR1_AA == this->mode)
        return this->cgBarR1().tot();
    // 2D integrations
    this->setX(a2);
    // cgBarF1
    if (Mode_cgBarF1_VV == this->mode || Mode_cgBarF1_VA == this->mode || Mode_cgBarF1_AA == this->mode)
        return this->cgBarF1().tot();
    // cgBar1
    if (Mode_cgBar1_VV == this->mode || Mode_cgBar1_VA == this->mode || Mode_cgBar1_AA == this->mode) {
        cuint old = this->mode;
        uint F = 0,R = 0;
        if (Mode_cgBar1_VV == this->mode)      { F = this->Mode_cgBarF1_VV; R = this->Mode_cgBarR1_VV; }
        else if (Mode_cgBar1_VA == this->mode) { F = this->Mode_cgBarF1_VA; R = this->Mode_cgBarR1_VA; }
        else if (Mode_cgBar1_AA == this->mode) { F = this->Mode_cgBarF1_AA; R = this->Mode_cgBarR1_AA; }
        dbl r = 0.;
        this->mode = F;
        r += this->cgBarF1().tot();
        this->mode = R;
        r += this->cgBarR1().tot();
        this->mode = old;
        return r;
    }
    // cqBarF1
    if (Mode_cqBarF1_VV == this->mode || Mode_cqBarF1_VA == this->mode || Mode_cqBarF1_AA == this->mode)
        return this->cqBarF1().tot();
    // 4D integrations
    this->setTheta2(a3);
    this->setY(a4);
    // cg1
    if (Mode_cg1_VV == this->mode || Mode_cg1_VA == this->mode || Mode_cg1_AA == this->mode)
        return this->cg1().tot();
    // cq1
    if (Mode_cq1_VV == this->mode || Mode_cq1_VA == this->mode || Mode_cq1_AA == this->mode)
        return this->cq1().tot();
    // dq1
    if (Mode_dq1_VV == this->mode || Mode_dq1_VA == this->mode || Mode_dq1_AA == this->mode)
        return this->dq1().tot();
    return 0.;
}

cdbl FullyDiff::IntKer::runHadronic(cdbl a1, cdbl a2, cdbl a3, cdbl a4, cdbl a5) {
    dbl r = 0.;
    // 2D integrations
    this->setZ(a1);
    this->setTheta1(a2);
    // LO
    if (this->flags.useLeadingOrder && this->flags.useGluonicChannel) {
        PhasespacePoint p(this->m2, this->Q2, this->xBj, this->muR2, this->muF2);
        p.setupLO(this->z, this->Theta1);
        cdbl curMuR2 = p.getMuR2();
        cdbl curMuF2 = p.getMuF2();
        cdbl alphaS = this->aS->alphasQ2(curMuR2);
        cdbl nLO = alphaS/(4.*M_PI*M_PI) * Q2/m2;
        cdbl xi = this->xBj/this->z;
        cdbl nLOg = this->V_z * 1./this->z * this->pdf->xfxQ2(21,xi,curMuF2);
        cdbl cg0 = this->cg0();
        cdbl fLOg = nLO * nLOg * cg0;
        this->fillAllOrderHistograms(p, fLOg);
        r += fLOg;
    }
    // NLO
    if (this->flags.useNextToLeadingOrder) {
        // 5D integrations
        this->setX(a3);
        this->setY(a4);
        this->setTheta2(a5);
        if (this->flags.useGluonicChannel) { // NLOg
            const PhasespaceValues cg1 = this->cg1();
            const PhasespaceValues cgBarR1 = this->cgBarR1();
            const PhasespaceValues cgBarF1 = this->cgBarF1();
            r += this->combineNLOg(this->xE, this->yE, cg1.xEyE, cgBarR1.xEyE, cgBarF1.xEyE);
            r += this->combineNLOg(1.,       this->yE, cg1.xCyE, cgBarR1.xCyE, cgBarF1.xCyE);
            r += this->combineNLOg(this->xE, -1.,      cg1.xEyC, cgBarR1.xEyC, cgBarF1.xEyC);
            r += this->combineNLOg(1.,       -1.,      cg1.xCyC, cgBarR1.xCyC, cgBarF1.xCyC);
        }
        if (this->flags.useQuarkChannel) { // NLOq
            const PhasespaceValues cq1 = this->cq1();
            const PhasespaceValues cqBarF1 = this->cqBarF1();
            // collect dq1 and oq1 in advance
            cuint oldMode = this->mode;
            PhasespaceValues dq1_VV;PhasespaceValues dq1_VA;PhasespaceValues dq1_AA;
            PhasespaceValues oq1_VV;PhasespaceValues oq1_VA;PhasespaceValues oq1_AA;
            if (isParityConservingProj(this->proj)) {
                this->mode = Mode_dq1_VV;
                dq1_VV = this->dq1();
                this->mode = Mode_oq1_VV;
                oq1_VV = this->oq1();
                if (this->flags.useZ) {
                    this->mode = Mode_dq1_AA;
                    dq1_AA = this->dq1();
                    this->mode = Mode_oq1_AA;
                    oq1_AA = this->oq1();
                }
            } else {
                this->mode = Mode_dq1_VA;
                dq1_VA = this->dq1();
                this->mode = Mode_oq1_VA;
                oq1_VA = this->oq1();
            }
            this->mode = oldMode;
            r += this->combineNLOq(this->xE, this->yE, cq1.xEyE, cqBarF1.xEyE, dq1_VV.xEyE, dq1_VA.xEyE, dq1_AA.xEyE, oq1_VV.xEyE, oq1_VA.xEyE, oq1_AA.xEyE);
            r += this->combineNLOq(this->xE, -1.,      cq1.xEyC, cqBarF1.xEyC, dq1_VV.xEyC, dq1_VA.xEyC, dq1_AA.xEyC, oq1_VV.xEyC, oq1_VA.xEyC, oq1_AA.xEyC);
        }
    }
    return r;
}

cdbl FullyDiff::IntKer::operator()(cdbl a1, cdbl a2, cdbl a3, cdbl a4, cdbl a5) {      
    // void mode
    if (0 == this->mode) return 0.;
    // partonic mode?
    if (this->mode < Mode_F) {
        cdbl r = this->runPartonic(a1,a2,a3,a4);
///if(!isfinite(r)) printf("<%.10e\t%.10e ()\n",a1,a2);
        return isfinite(r) ? r : 0.;
    }
    // hadronic mode?
    if (this->mode < Mode_sigma) {
        cdbl r = this->runHadronic(a1,a2,a3,a4,a5);
///if(!isfinite(r)) printf("<%.10e\t%.10e ()\n",a1,a2);
        return isfinite(r) ? r : 0.;
    }
    return 0.;
}

void FullyDiff::IntKer::operator()(const double x[], const int k[], const double& weight, const double aux[], double f[]) {
    this->vegasWeight = &weight;
    cdbl a0 = this->dim >= 1 ? x[0] : 0.;
    cdbl a1 = this->dim >= 2 ? x[1] : 0.;
    cdbl a2 = this->dim >= 3 ? x[2] : 0.;
    cdbl a3 = this->dim >= 4 ? x[3] : 0.;
    cdbl a4 = this->dim >= 5 ? x[4] : 0.;
    cdbl i = this->operator()(a0,a1,a2,a3,a4);
    f[0] = i;
}

#define combineNLOInit \
    PhasespacePoint p(this->m2, this->Q2, this->xBj, this->muR2, this->muF2);\
    p.setupNLO(this->z,x,y,this->Theta1,this->Theta2);\
    cdbl muR2 = p.getMuR2();\
    cdbl muF2 = p.getMuF2();\
    cdbl alphaS = this->aS->alphasQ2(muR2);\
    cdbl nNLO = alphaS*alphaS/(M_PI) * Q2/m2;\
    cdbl xi = this->xBj/z;

cdbl FullyDiff::IntKer::combineNLOg(cdbl x, cdbl y, cdbl cg1, cdbl cgBarR1, cdbl cgBarF1) {
    combineNLOInit
    // combine
    cdbl nNLOg = this->V_z * 1./this->z * this->pdf->xfxQ2(21,xi,muF2);
    cdbl fNLOg = nNLO * nNLOg * (cg1 + log(muR2/this->m2)*cgBarR1 + log(muF2/this->m2)*cgBarF1);
    if (!isfinite(fNLOg) || 0. == fNLOg)
        return 0.;
    // fill histograms
    this->fillAllOrderHistograms(p, fNLOg);
    this->fillNLOHistograms(p, fNLOg);
    return fNLOg;
}

cdbl FullyDiff::IntKer::combineNLOq(cdbl x, cdbl y, cdbl cq1, cdbl cqBarF1, cdbl dq1_VV, cdbl dq1_VA, cdbl dq1_AA, cdbl oq1_VV, cdbl oq1_VA, cdbl oq1_AA) {
    combineNLOInit
    cdbl eH = this->getElectricCharge(this->nlf+1);
    cdbl gVQ = this->getVectorialCoupling(this->nlf+1);
    cdbl gAQ = this->getAxialCoupling(this->nlf+1);
    // join elems
    cdbl cq1j = cq1 + log(muF2/m2)*cqBarF1;
    dbl fqs = 0.;
    for (uint q = 1; q < this->nlf + 1; ++q) {
        dbl dq1 = 0., oq1 = 0.;
        cdbl eL = this->getElectricCharge(q);
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
        fqs += (this->pdf->xfxQ2((int)q,xi,muF2) + this->pdf->xfxQ2(-((int)q),xi,muF2)) * (cq1j + dq1 + oq1);
    }
    cdbl fNLOq = nNLO * (this->V_z * 1./this->z) * fqs ;
    if (!isfinite(fNLOq) || 0. == fNLOq)
        return 0.;
    // fill histograms
    this->fillAllOrderHistograms(p, fNLOq);
    this->fillNLOHistograms(p, fNLOq);
    return fNLOq;
}

void FullyDiff::IntKer::scaleHistograms(cdbl s) const {
    for (histMapT::const_iterator it = this->histMap.cbegin(); it != this->histMap.cend(); ++it)
        it->second->scale(s);
}

void FullyDiff::IntKer::Dvegas_init() const {
    this->scaleHistograms(0.);
}
    
void FullyDiff::IntKer::Dvegas_final(cuint iterations) const {
    this->scaleHistograms(1./cdbl(iterations));
}

#define fillTemplate(cases)\
/* something active?*/\
if (this->histMap.empty()) return;\
if (0 == this->vegasWeight) return;\
cdbl value = i*(*this->vegasWeight);\
if (!isfinite(value) || 0. == value) return;\
for (histMapT::const_iterator it = this->histMap.cbegin(); it != this->histMap.cend(); ++it) {\
    dbl var = dblNaN;\
    switch (it->first) {\
        cases\
        default: continue;\
    }\
    it->second->accumulate(var,value);\
}

void FullyDiff::IntKer::fillAllOrderHistograms(const PhasespacePoint& p, cdbl i) const {
    fillTemplate(
        case histT::log10z:         var = p.getZ();            break;
        case histT::log10xi:        var = this->xBj/p.getZ();  break;
        case histT::Theta1:         var = p.getTheta1();       break;
        
        case histT::HQPairInvMass:            var = (p.getP1() + p.getP2()).m();               break;
        case histT::HQPairDeltaPhi:           var = abs(geom3::deltaPhi(p.getP1(),p.getP2())); break;
        case histT::HQPairTransverseMomentum: var = (p.getP1() + p.getP2()).pt();              break; 
        case histT::HQPairConeSizeVariable:   var = geom3::deltaR(p.getP1(),p.getP2());        break;                      
            
        case histT::HAQRapidity:                  var = p.getP2().rapidity();     break;
        case histT::HAQTransverseMomentum:        var = p.getP2().pt();           break;
        case histT::HAQTransverseMomentumScaling:
            {cdbl ptmax = sqrt(p.getSh()/4.-m2);  var = p.getP2().pt()/ptmax;}    break;
        case histT::HAQFeynmanX:
            {cdbl plmax = sqrt(p.getSh()/4.-m2);  var = p.getP2().pz()/plmax;}    break;
    )
}

void FullyDiff::IntKer::fillNLOHistograms(const PhasespacePoint& p, cdbl i) const {
    if (!p.isNLO())
        return;
    fillTemplate(
        case histT::x:      var = p.getX();         break;
        case histT::y:      var = p.getY();         break;
        case histT::Theta2: var = p.getTheta2();    break;
    )
}
