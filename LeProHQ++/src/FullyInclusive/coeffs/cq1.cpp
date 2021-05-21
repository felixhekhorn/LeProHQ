#include "cq1.h"

#include "gslpp/gslpp.Interpolation1D.hpp"
#include "gslpp/gslpp.Interpolation2D.hpp"
#include "Common/Color.hpp"
#include "Common/EnvUtils.h"
#include "cg0t.h"

#define a11_ 1.0
#define a10_ -13./12. + 3./2.*M_LN2

// true threshold
cdbl FullyInclusive::coeffs::cq1t_(cdbl m2, cdbl q2, cdbl s, fPtr3dbl cg0t, cdbl a11, cdbl a10) {
    partonic_beta(m2,s)
    partonic_rhoq(m2,q2)
    cdbl lo = cg0t(m2,q2,s);
    cdbl n = beta*beta*rhoq/(M_PI*M_PI*(rhoq-1.)) * Color::Kqph/6./Color::Kgph;
    cdbl res = a11 * log(beta) + a10;
    return lo * n * res;
}

cdbl FullyInclusive::coeffs::cq1t_F2_VV(cdbl m2, cdbl q2, cdbl s) {
    cdbl a11 = a11_;
    cdbl a10 = a10_;
    return cq1t_(m2,q2,s,&cg0t_F2_VV,a11,a10);
}

cdbl FullyInclusive::coeffs::cq1t_FL_VV(cdbl m2, cdbl q2, cdbl s) {
    cdbl a11 = a11_ - 2./5.;
    cdbl a10 = -77./100. + 9./10.*M_LN2;
    return cq1t_(m2,q2,s,&cg0t_FL_VV,a11,a10);
}

cdbl FullyInclusive::coeffs::cq1t_x2g1_VV(cdbl m2, cdbl q2, cdbl s) {
    cdbl a11 = a11_;
    cdbl a10 = a10_ - 1./4.;
    return cq1t_(m2,q2,s,&cg0t_x2g1_VV,a11,a10);
}

cdbl FullyInclusive::coeffs::cq1t_F2_AA(cdbl m2, cdbl q2, cdbl s) {
    cdbl a11 = a11_;
    cdbl a10 = a10_;
    return cq1t_(m2,q2,s,&cg0t_F2_AA,a11,a10);
}

cdbl FullyInclusive::coeffs::cq1t_FL_AA(cdbl m2, cdbl q2, cdbl s) {
    cdbl a11 = a11_;
    cdbl a10 = a10_;
    return cq1t_(m2,q2,s,&cg0t_FL_AA,a11,a10);
}

cdbl FullyInclusive::coeffs::cq1t_x2g1_AA(cdbl m2, cdbl q2, cdbl s) {
    cdbl a11 = a11_;
    cdbl a10 = a10_;
    return cq1t_(m2,q2,s,&cg0t_x2g1_AA,a11,a10);
}

// improved threshold
cdbl FullyInclusive::coeffs::cq1tp_(cdbl m2, cdbl q2, cdbl s, fPtr3dbl cq1t, str proj_cc) {
    // true threshold
    cdbl t = cq1t(m2, q2, s);
    cdbl lnxi = log(-q2/m2);
    // load additional factor
    const boost::filesystem::path parent = Common::getPathByEnv("PARTONIC_GRIDS");
    const boost::filesystem::path mine = parent / "cq1" / ("cq1-"+proj_cc+"-thres-coeff.dat");
    if (!boost::filesystem::exists(mine))
        throw runtime_error("path does not exist: "+mine.string());
    gslpp::Interpolation1D a_int(gsl_interp_steffen,mine.string(),FullyInclusive::coeffs::cq1tp_lnxi_num);
    cdbl a = a_int.eval(lnxi);
    partonic_eta(m2, s)
    return t*(1. + a*eta);
}

// bulk
cdbl FullyInclusive::coeffs::cq1b_(cdbl m2, cdbl q2, cdbl s, str proj_cc) {
    // load additional factor
    const boost::filesystem::path parent = Common::getPathByEnv("PARTONIC_GRIDS");
    const boost::filesystem::path mine = parent / "cq1" / ("cq1-"+proj_cc+"-bulk.dat");
    if (!boost::filesystem::exists(mine))
        throw runtime_error("path does not exist: "+mine.string());
    gslpp::Interpolation2D bulk_int(gsl_interp2d_bicubic,mine.string(),FullyInclusive::coeffs::cq1b_lneta_num,FullyInclusive::coeffs::cq1b_lnxi_num);
    cdbl lnxi = log(-q2/m2);
    partonic_eta(m2, s)
    cdbl res = bulk_int.eval(log(eta),lnxi);
    return res;
}

// full
cdbl FullyInclusive::coeffs::cq1_(cdbl m2, cdbl q2, cdbl s, fPtr3dbl cq1tp, fPtr3dbl cq1b) {
    partonic_eta(m2, s)
    cdbl lneta = log(eta);
    cdbl low = FullyInclusive::coeffs::cq1b_lneta_min;
    // threshold only?
    if (lneta < low)
        return cq1tp(m2, q2, s);
    cdbl high = FullyInclusive::coeffs::cq1tp_lneta_max;
    // bulk only?
    if (lneta > high)
        return cq1b(m2, q2, s);
    // otherwise apply linear interpolation between the two
    cdbl tp = cq1tp(m2, q2, s);
    cdbl b = cq1b(m2, q2, s);
    return (tp*(lneta - high) + b*(low - lneta))/(low - high);
}

// keep all parity conserving elements
#define cq1_proj_pc(proj)\
cdbl FullyInclusive::coeffs::cq1tp_##proj(cdbl m2, cdbl q2, cdbl s) {\
    return cq1tp_(m2, q2, s, &cq1t_##proj, #proj);\
}\
cdbl FullyInclusive::coeffs::cq1b_##proj(cdbl m2, cdbl q2, cdbl s) {\
    return cq1b_(m2, q2, s, #proj);\
}\
cdbl FullyInclusive::coeffs::cq1_##proj(cdbl m2, cdbl q2, cdbl s) {\
    return cq1_(m2, q2, s, &cq1tp_##proj, &cq1b_##proj);\
}
interateAllPCProj(cq1_proj_pc)

// we assume all parity violating elements to be 0 (likely, but not proven)
#define cq1_proj_pv(proj) \
cdbl FullyInclusive::coeffs::cq1t_##proj(cdbl m2, cdbl q2, cdbl s) { return 0.; } \
cdbl FullyInclusive::coeffs::cq1tp_##proj(cdbl m2, cdbl q2, cdbl s) { return 0.; } \
cdbl FullyInclusive::coeffs::cq1b_##proj(cdbl m2, cdbl q2, cdbl s) { return 0.; } \
cdbl FullyInclusive::coeffs::cq1_##proj(cdbl m2, cdbl q2, cdbl s) { return 0.; }

interateAllPVProj(cq1_proj_pv)