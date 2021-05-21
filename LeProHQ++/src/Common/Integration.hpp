/**
 * @brief defines integration routines in 1-5 dimensions
 */

#ifndef Integration_H_
#define Integration_H_

#include <gsl/gsl_integration.h>
#include <gsl/gsl_monte_vegas.h>
#include <dvegas.h>

#include "config.h"
#include "Common/IntegrationMeta.hpp"
#include "gslpp/gslpp.Functor.hpp"

namespace Common {

/**
 * @brief integrates the kernel in 1 dimension
 * @param K kernel
 * @param cfg config
 * @param out output of meta data
 * @return \f$\int\limits_0^1 f(a)\,da\f$
 */
template <class IntKerT> cdbl integrate1D(IntKerT* K, const IntegrationConfig& cfg, IntegrationOutput* out)  {
    // init
    if ("gsl_integration_qag" != cfg.method && "gsl_integration_cquad" != cfg.method)
        throw invalid_argument("currently only gsl_integration_qag and gsl_integration_cquad is supported for 1D-integration!");
    gsl_function f;
    f.params = K;
    f.function = gslpp::callFunctor<IntKerT>;
    out->result = dblNaN;
    out->error = dblNaN;
    // run
    dbl res,err;
    if ("gsl_integration_qag" == cfg.method) {
        size_t calls = cfg.calls;
        gsl_integration_workspace *w = gsl_integration_workspace_alloc(calls);
        gsl_integration_qag(&f, 0., 1., cfg.GslQag_epsabs, cfg.GslQag_epsrel, calls, cfg.GslQag_key, w, &res,&err);
        gsl_integration_workspace_free(w);
        // out
        if (cfg.verbosity > 0)
            printf("[INFO] int1D(gsl_qag): % e ± %e (%.3f%%)\n",res,err,fabs(err/res*1e2));
    } else if ("gsl_integration_cquad" == cfg.method) {
        size_t calls;
        gsl_integration_cquad_workspace *w = gsl_integration_cquad_workspace_alloc(cfg.calls);
        gsl_integration_cquad(&f, 0.,1., cfg.GslQag_epsabs,cfg.GslQag_epsrel, w, &res,&err, &calls);
        gsl_integration_cquad_workspace_free(w);
        // out
        if (cfg.verbosity > 0)
            printf("[INFO] int1D(gsl_cquad): % e ± %e (%.3f%%) calls: %lu\n",res,err,fabs(err/res*1e2),calls);
    }
    out->result = res;
    out->error = err;
    //size_t reNevals;reNevals = w->size;
    return res;
}

// call gsl
#define runGslVegas \
    /* init */\
    F.dim = dim;\
    const gsl_rng_type *T;\
    gsl_rng *r;\
    gsl_rng_env_setup ();\
    T = gsl_rng_default;\
    r = gsl_rng_alloc (T);\
    double res,err;\
    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (dim);\
    /* warmup */\
    gsl_monte_vegas_params p;\
    gsl_monte_vegas_params_get(s,&p);\
    p.iterations = cfg.MC_warmupIterations;\
    p.verbose = cfg.verbosity - 3;\
    gsl_monte_vegas_params_set(s,&p);\
    gsl_monte_vegas_integrate (&F, xl, xu, dim, cfg.MC_warmupCalls, r, s, &res, &err);\
    /* run */\
    p.iterations = cfg.MC_iterations;\
    p.verbose = cfg.verbosity - 2;\
    gsl_monte_vegas_params_set(s,&p);\
    uint guard = 0;\
    out->result = dblNaN;\
    out->error = dblNaN;\
    out->MC_chi2 = dblNaN;\
    out->MC_chi2inter = 0;\
    if (cfg.MC_adaptChi2) {\
        do {\
            if (!isfinite(res)) return res;\
            gsl_monte_vegas_integrate (&F, xl, xu, dim, cfg.calls, r, s, &res, &err);\
            if (cfg.verbosity > 1)\
                printf("[INFO] int%dD(gsl):    [%d] % e ± %.3e (%.3f%%) chi2/it: %.3f\n",dim,guard,res,err,fabs(err/res*1e2),gsl_monte_vegas_chisq(s));\
        } while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5 && ++guard < 15);\
    } else {\
        gsl_monte_vegas_integrate (&F, xl, xu, dim, cfg.calls, r, s, &res, &err);\
    }\
    out->MC_chi2 = gsl_monte_vegas_chisq (s);\
    /* out */\
    if (1 == cfg.verbosity)\
        printf("[INFO] int%dD(gsl):    [%d] % e ± %.3e (%.3f%%) chi2/it: %.3f\n",dim,guard,res,err,fabs(err/res*1e2),gsl_monte_vegas_chisq(s));\
    gsl_monte_vegas_free (s);\
    gsl_rng_free (r);\
    out->result = res;\
    out->error = err;\
    out->MC_chi2inter = guard;\
    return res;
    
#define runDvegas\
    HepSource::Dvegas dv(dim,cfg.Dvegas_bins,1,*K);\
    /** @todo activate correlation between z and x? -> Dvegas dv(dim,cfg.Dvegas_bins,2,{},0,1,F); */\
    double res, err;\
    /* clear histograms */\
    K->Dvegas_init();\
    /* warm-up */\
    /* catch zero kernel */\
    try{ HepSource::VEGAS(dv,cfg.MC_warmupCalls,cfg.MC_warmupIterations,0,cfg.verbosity - 3); }\
    catch(domain_error& e){ out->result = 0; out->error = 0; out->MC_chi2 = 0; out->MC_chi2inter = 0; return 0.; }\
    HepSource::IntegrandEstimate e = dv.stats(0);\
    res = e.integral();\
    uint guard = 0;\
    out->result = dblNaN;\
    out->error = dblNaN;\
    out->MC_chi2 = dblNaN;\
    out->MC_chi2inter = 0;\
    /* run */\
    if (cfg.MC_adaptChi2) { /* adapt chi */\
        do {\
            if (!isfinite(res)) return res;\
            K->Dvegas_init();\
            HepSource::VEGAS(dv,cfg.calls,cfg.MC_iterations,1,cfg.verbosity - 2);\
            e = dv.stats(0);\
            res = e.integral();\
            err = e.standardDeviation();\
            if (cfg.verbosity > 1)\
                printf("[INFO] int%dD(Dvegas): [%d] % e ± %.3e (%.3f%%) chi2/it: %.3f\n",dim,guard,res,err,fabs(err/res*1e2),e.chiSquarePerIteration());\
        } while (fabs (e.chiSquarePerIteration() - 1.0) > 0.5 && ++guard < 15);\
    } else { /* simple run */\
        K->Dvegas_init();\
        HepSource::VEGAS(dv,cfg.calls,cfg.MC_iterations,1,cfg.verbosity - 2);\
        e = dv.stats(0);\
        res = e.integral();\
        err = e.standardDeviation();\
    }\
    K->Dvegas_final(cfg.MC_iterations);\
    if (1 == cfg.verbosity)\
        printf("[INFO] int%dD(Dvegas): [%d] % e ± %.3e (%.3f%%) chi2/it: %.3f\n",dim,guard,res,err,fabs(err/res*1e2),e.chiSquarePerIteration());\
    out->result = res;\
    out->error = err;\
    out->MC_chi2 = e.chiSquarePerIteration();\
    out->MC_chi2inter = guard;\
    return res;

/**
 * @brief integrates the kernel in 2 dimension
 * @param K kernel
 * @param cfg config
 * @param out output of meta data
 * @return \f$\int\limits_0^1 f(a_1,a_2)\,da_1da_2\f$
 */
template <class IntKerT> cdbl integrate2D(IntKerT* K, const IntegrationConfig& cfg, IntegrationOutput* out) {
    cuint dim = 2;
    if ("gsl_monte_vegas_integrate" == cfg.method) {
        gsl_monte_function F;
        F.params = K;
        F.f = gslpp::callFunctor2D<IntKerT>;
        double xl[dim] = {0., 0.};
        double xu[dim] = {1., 1.};
        runGslVegas
    } else if ("Dvegas" == cfg.method) {
        runDvegas
    }
    throw invalid_argument("currently only gsl_monte_vegas_integrate or Dvegas is supported for 2D-integration!");
}

/**
 * @brief integrates the kernel in 3 dimension
 * @param K kernel
 * @param cfg config
 * @param out output of meta data
 * @return \f$\int\limits_0^1 f(a_1,a_2,a_3)\,da_1da_2da_3\f$
 */
template <class IntKerT> cdbl integrate3D(IntKerT* K, const IntegrationConfig& cfg, IntegrationOutput* out) {
    cuint dim = 3;
    if ("gsl_monte_vegas_integrate" == cfg.method) {
        gsl_monte_function F;
        F.params = K;
        F.f = gslpp::callFunctor3D<IntKerT>;
        double xl[dim] = {0., 0., 0.};
        double xu[dim] = {1., 1., 1.};
        runGslVegas
    } else if ("Dvegas" == cfg.method) {
        runDvegas
    }
    throw invalid_argument("currently only gsl_monte_vegas_integrate or Dvegas is supported for 3D-integration!");
}

/**
 * @brief integrates the kernel in 4 dimension
 * @param K kernel
 * @param cfg config
 * @param out output of meta data
 * @return \f$\int\limits_0^1 f(a_1,a_2,a_3,a_4)\,da_1da_2da_3da_4\f$
 */
template <class IntKerT> cdbl integrate4D(IntKerT* K, const IntegrationConfig& cfg, IntegrationOutput* out) {
    cuint dim = 4;
    if ("gsl_monte_vegas_integrate" == cfg.method) {
        gsl_monte_function F;
        F.params = K;
        F.f = gslpp::callFunctor4D<IntKerT>;
        double xl[dim] = {0., 0., 0., 0.};
        double xu[dim] = {1., 1., 1., 1.};
        runGslVegas
    } else if ("Dvegas" == cfg.method) {
        runDvegas
    }
    throw invalid_argument("currently only gsl_monte_vegas_integrate or Dvegas is supported for 4D-integration!");
}

/**
 * @brief integrates the kernel in 5 dimension
 * @param K kernel
 * @param cfg config
 * @param out output of meta data
 * @return \f$\int\limits_0^1 f(a_1,a_2,a_3,a_4,a_5)\,da_1da_2da_3da_4da_5\f$
 */
template <class IntKerT> cdbl integrate5D(IntKerT* K, const IntegrationConfig& cfg, IntegrationOutput* out) {
    cuint dim = 5;
    if ("gsl_monte_vegas_integrate" == cfg.method) {
        gsl_monte_function F;
        F.params = K;
        F.f = gslpp::callFunctor5D<IntKerT>;
        double xl[dim] = {0., 0., 0., 0., 0.};
        double xu[dim] = {1., 1., 1., 1., 0.};
        runGslVegas
    } else if ("Dvegas" == cfg.method) {
        runDvegas
    }
    throw invalid_argument("currently only gsl_monte_vegas_integrate or Dvegas is supported for 5D-integration!");
}

} // namespace Common

#endif // Integration_H_
