#ifndef DynamicScaleFactors_HPP_
#define DynamicScaleFactors_HPP_

#include <yaml-cpp/yaml.h>

/**
 * @brief computes dynamic scales
 */
struct DynamicScaleFactors {
    
/** @brief factor to m2 */
    dbl cM2;
    
/** @brief factor to Q2 */
    dbl cQ2;

/** @brief factor to \f$p_{\HepGenAntiParticle{Q}{}{},T}^2\f$ */
    dbl cHAQTransverseMomentum;

/** @brief factor to \f$(p_{\HepGenParticle{Q}{}{},T}+p_{\HepGenAntiParticle{Q}{}{}},T)^2\f$ */
    dbl cHQPairTransverseMomentum;
    
/**
 * @brief constructor
 * @param cM2 factor to m2
 * @param cQ2 factor to Q2
 * @param cHAQTransverseMomentum factor to \f$p_{\HepGenAntiParticle{Q}{}{},T}^2\f$
 * @param cHQPairTransverseMomentum factor to \f$(p_{\HepGenParticle{Q}{}{},T}+p_{\HepGenAntiParticle{Q}{}{}},T)^2\f$
 */
    DynamicScaleFactors(cdbl cM2, cdbl cQ2, cdbl cHAQTransverseMomentum, cdbl cHQPairTransverseMomentum) :
        cM2(cM2), cQ2(cQ2), cHAQTransverseMomentum(cHAQTransverseMomentum), cHQPairTransverseMomentum(cHQPairTransverseMomentum) {};
/**
 * @brief return YAML representation as string
 * @return string
 */
    str toYAML() const {
        YAML::Node r;
        r["cM2"] = this->cM2;
        r["cQ2"] = this->cQ2;
        r["cHAQTransverseMomentum"] = this->cHAQTransverseMomentum;
        r["cHQPairTransverseMomentum"] = this->cHQPairTransverseMomentum;
        r.SetStyle(YAML::EmitterStyle::Flow);
        YAML::Emitter e;
        e << r;
        return e.c_str(); 
    }

/**
 * @brief writes object as YAML to stream
 * @param os stream
 * @param c object
 * @return improved stream
 */
    friend std::ostream& operator<<(std::ostream& os, const DynamicScaleFactors& c ) {
        os << c.toYAML();
        return os;
    }
    
/**
 * @brief multiplicate each element with double
 * @param m old
 * @param a factor
 * @return m*a
 */
    friend const DynamicScaleFactors operator*(const DynamicScaleFactors m, cdbl a) {
        return DynamicScaleFactors(a*m.cM2, a*m.cQ2, a*m.cHAQTransverseMomentum, a*m.cHQPairTransverseMomentum);
    }
    
/**
 * @brief multiplicate each element with double
 * @param a factor
 * @param m old
 * @return a*m = m*a
 */
    friend const DynamicScaleFactors operator*(cdbl a, const DynamicScaleFactors m) { return m*a; }
};

#endif // DynamicScaleFactors_HPP_