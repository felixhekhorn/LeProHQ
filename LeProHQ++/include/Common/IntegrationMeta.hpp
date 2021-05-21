/**
 * @brief defines meta classes for integration routines
 */

#ifndef IntegrationMeta_H_
#define IntegrationMeta_H_

#include <gsl/gsl_integration.h>
#include <yaml-cpp/yaml.h>

#include "config.h"

namespace Common {

/**
 * @class IntegrationConfig
 * @brief configurates a single integration
 */
struct IntegrationConfig {
    
/** @name common variables */
///@{
/** @brief used integration method */
    str method;
    
/** @brief level of output */
    int verbosity = 0;

/** @brief calls */
    size_t calls = 0;
///@}
    
/** @name Monte Carlo variables */
///@{
/** @brief calls for warmup */
    size_t MC_warmupCalls = 0;
    
/** @brief iterations */
    uint MC_iterations = 5;
    
/** @brief iterations during warmup */
    uint MC_warmupIterations = 5;
    
/** @brief iterate until |chi2-1| < 0.5? */
    bool MC_adaptChi2 = true;
///@}

/** @name variables for Dvegas */
///@{
/** @brief number of bins */
    uint Dvegas_bins = 250;
///@}

/** @name variables for gsl_integration_qag */
///@{
/** @see gsl_integration_qag */
    dbl GslQag_epsabs = 5e-6;
/** @see gsl_integration_qag */
    dbl GslQag_epsrel = 1e-5;
/** @see gsl_integration_qag */
    int GslQag_key = GSL_INTEG_GAUSS41;
///@}

/**
 * @brief return YAML representation as string
 * @return string
 */
    str toYAML() const {
        YAML::Node r;
        r["method"] = method;
        r["verbosity"] = verbosity;
        r["calls"] = calls;
        if ("gsl_integration_qag" == this->method) {
            r["GslQag_epsabs"] = GslQag_epsabs;
            r["GslQag_epsrel"] = GslQag_epsrel;
            r["GslQag_key"] = GslQag_key;
        } else {
            r["MC_warmupCalls"] = MC_warmupCalls;
            r["MC_iterations"] = MC_iterations;
            r["MC_warmupIterations"] = MC_warmupIterations;
            r["MC_adaptChi2"] = MC_adaptChi2;
            if ("Dvegas" == this->method) {
                r["Dvegas_bins"] = Dvegas_bins;
            }
        }
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
    friend std::ostream& operator<<(std::ostream& os, const IntegrationConfig& c ) {
        os << c.toYAML();
        return os;
    }

};

/**
 * @class IntegrationOutput
 * @brief stores meta data for a single integration
 */
struct IntegrationOutput {
    
/** @name common variables */
///@{
/** @brief result */
    dbl result = 0;
/** @brief (absolute) error */
    dbl error = 0;
///@}
    
/** @name Monte Carlo variables */
///@{
/** @brief chi^2 */
    dbl MC_chi2 = 0;
/** @brief number of iteration to converge chi^2 */
    uint MC_chi2inter = 0;
///@}
    
/**
 * @brief constructor
 * @param result result
 * @param error (absolute) error
 * @param MC_chi2 chi^2
 * @param MC_chi2inter number of iteration to converge chi^2
 */
    IntegrationOutput(cdbl result = 0, cdbl error = 0, cdbl MC_chi2 = 0, cuint MC_chi2inter = 0) : 
        result(result), error(error), MC_chi2(MC_chi2), MC_chi2inter(MC_chi2inter) {};

/**
 * @brief return YAML representation as string
 * @return string
 */
    str toYAML() const {
        YAML::Node r;
        r["result"] = result;
        r["error"] = error;
        r["MC_chi2"] = MC_chi2;
        r["MC_chi2inter"] = MC_chi2inter;
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
    friend std::ostream& operator<<(std::ostream& os, const IntegrationOutput& c ) {
        os << c.toYAML();
        return os;
    }
};

} // namespace Common

#endif // IntegrationMeta_H_