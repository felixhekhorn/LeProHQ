#include "FullyDiff/PhasespacePoint.h"

#include <boost/format.hpp>

#include <rk/geom3.hh>
#include <gsl/gsl_rng.h>

FullyDiff::PhasespacePoint::PhasespacePoint(cdbl m2, cdbl Q2, cdbl bjorkenX, const DynamicScaleFactors& muR2Factors, const DynamicScaleFactors& muF2Factors) :
    m2(m2), Q2(Q2), bjorkenX(bjorkenX), muR2Factors(muR2Factors), muF2Factors(muF2Factors){
}

void FullyDiff::PhasespacePoint::setupLO(cdbl z, cdbl Theta1) {
    this->order = 0;
    this->z = z;
    this->sp = Q2/z;
    this->Theta1 = Theta1;
    
    // use virtual photon-parton c.m.s.
    using rk::P4;
    using geom3::UnitVector3;
    cdbl s = this->sp - this->Q2;
    cdbl sqrts = sqrt(s);
    cdbl beta = sqrt(1. - 4.*m2/s);
    this->k1 = P4(this->sp/(2.*sqrts) * UnitVector3::zAxis(),0.);
    this->q = P4((s - this->Q2)/(2.*sqrts),-k1.momentum()); // space-like!
    const UnitVector3 u(0., sin(this->Theta1),cos(this->Theta1));
    this->p1 = P4(-.5*sqrts*beta*u,sqrt(m2));
    this->p2 = P4(.5*sqrts*beta*u,sqrt(m2));
    
    // move to final frame
    this->applyLTsToFinalFrame();
}

void FullyDiff::PhasespacePoint::setupNLO(cdbl z, cdbl x, cdbl y, cdbl Theta1, cdbl Theta2) {
    this->order = 1;
    this->z = z;
    this->sp = Q2/z;
    this->Theta1 = Theta1;
    this->x = x;
    this->y = y;
    this->Theta2 = Theta2;
    const KinematicVars vs(this->m2, -this->Q2, this->sp, x, y, Theta1, Theta2);
    using rk::P4;
    using geom3::UnitVector3;
    
    this->k1 = P4(vs.k10*UnitVector3(0.,vs.sinPsi,vs.cosPsi),0.);
    this->q = P4(vs.q0,vs.absq*UnitVector3::zAxis());
    const UnitVector3 u(sin(Theta2)*sin(Theta1),cos(Theta2)*sin(Theta1),cos(Theta1));
    this->p1 = P4(.5*sqrt(vs.s5)*vs.beta5*u,sqrt(m2));
    this->p2 = P4(-.5*sqrt(vs.s5)*vs.beta5*u,sqrt(m2));
    P4 k2 (k1.momentum() + q.momentum(),0.);
    
    { // boost to virtual photon-parton c.m.s.
        using rk::Boost;
        const P4 ps = q+k1;
        const Boost partonCMS = ps.restBoost();
        this->k1.boost(partonCMS);
        this->p1.boost(partonCMS);
        this->p2.boost(partonCMS);
        k2.boost(partonCMS);
        this->q = p1+p2+k2-k1; // trick, as space-like vectors can't be boosted
    } { // align k1 to z
        using geom3::Rotation3;
        using geom3::Vector3;
        const Vector3 k1vec = this->k1.momentum();
        const Rotation3 toZ(k1vec.cross(UnitVector3::zAxis()).direction(),k1vec.angle(UnitVector3::zAxis()));
        this->q.rotate(toZ);
        this->k1.rotate(toZ);
        this->p1.rotate(toZ);
        this->p2.rotate(toZ);
        //k2.rotate(toZ);
    }
    
    // move to final frame
    this->applyLTsToFinalFrame();
}

void FullyDiff::PhasespacePoint::applyLTsToFinalFrame() {
    { // randomize around z-axis
        using geom3::Rotation3;
        using geom3::UnitVector3;
        const gsl_rng_type *T;
        gsl_rng *r;
        gsl_rng_env_setup();
        T = gsl_rng_default;
        r = gsl_rng_alloc(T);
        const Rotation3 randZ (UnitVector3::zAxis(),2.*M_PI*gsl_rng_uniform(r));
        gsl_rng_free(r);
        this->q.rotate(randZ);
        this->k1.rotate(randZ);
        this->p1.rotate(randZ);
        this->p2.rotate(randZ);
    } { // boost to virtual photon-hadron c.m.s.
        using rk::P4;
        using rk::Boost;
        const P4 ph = (this->z/this->bjorkenX)*k1;
        const P4 pCMS = ph + q;
        const Boost hCMS = pCMS.restBoost();
        this->p1.boost(hCMS);
        this->p2.boost(hCMS);
    }
}

cdbl FullyDiff::PhasespacePoint::getDynamicScale(const DynamicScaleFactors& factors) const {
    cdbl HQPairTransverseMomentum = (this->p1 + this->p2).pt();
    cdbl HAQTransverseMomentum = this->p2.pt();
    cdbl mu2 = factors.cM2 * this->m2 +
                factors.cQ2 * this->Q2 + 
                factors.cHQPairTransverseMomentum * HQPairTransverseMomentum*HQPairTransverseMomentum +
                factors.cHAQTransverseMomentum * HAQTransverseMomentum*HAQTransverseMomentum;
    if (!isfinite(mu2) || mu2 <= 0.)
        throw domain_error((boost::format("all scales have to be finite and positive! (%e = %e*%e + %e*%e + %e*%e**2 + %e*%e**2)")%mu2%factors.cM2%m2%factors.cQ2%Q2%factors.cHQPairTransverseMomentum%HQPairTransverseMomentum%factors.cHAQTransverseMomentum%HAQTransverseMomentum).str());
    return mu2;
}

cdbl FullyDiff::PhasespacePoint::getZ() const {
    return this->z;
}

cdbl FullyDiff::PhasespacePoint::getSh() const {
    return this->Q2*(1./this->bjorkenX - 1.);
}

cdbl FullyDiff::PhasespacePoint::getTheta1() const {
    return this->Theta1;
}

const rk::P4 FullyDiff::PhasespacePoint::getP1() const {
    return this->p1;
}
    
const rk::P4 FullyDiff::PhasespacePoint::getP2() const {
    return this->p2;
}

cdbl FullyDiff::PhasespacePoint::getMuR2() const {
    return this->getDynamicScale(this->muR2Factors);
}
    
cdbl FullyDiff::PhasespacePoint::getMuF2() const {
    return this->getDynamicScale(this->muF2Factors);
}

const bool FullyDiff::PhasespacePoint::isNLO() const {
    return 1 == this->order;
}
    
cdbl FullyDiff::PhasespacePoint::getX() const {
    if (!this->isNLO())
        throw domain_error("only NLO phasespace points have a x value");
    return this->x;
}
    
cdbl FullyDiff::PhasespacePoint::getY() const {
    if (!this->isNLO())
        throw domain_error("only NLO phasespace points have a y value");
    return this->y;
}
    
cdbl FullyDiff::PhasespacePoint::getTheta2() const {
    if (!this->isNLO())
        throw domain_error("only NLO phasespace points have a Theta2 value");
    return this->Theta2;
}