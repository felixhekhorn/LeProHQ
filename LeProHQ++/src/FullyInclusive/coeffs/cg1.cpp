#include "cg1.h"

#include "gslpp/gslpp.Interpolation1D.hpp"
#include "gslpp/gslpp.Interpolation2D.hpp"
#include "Common/Color.hpp"
#include "Common/EnvUtils.h"
#include "cg0t.h"
#include "cg1_a10.h"

#define a12_ 1.0
#define a11_ -5./2. + 3.*M_LN2

// true threshold
cdbl FullyInclusive::coeffs::cg1t_(cdbl m2, cdbl q2, cdbl s, fPtr3dbl cg0t, cdbl a12, cdbl a11, fPtr2dbl a10_OK, fPtr2dbl a10_QED) {
    partonic_beta(m2,s)
    cdbl lo = cg0t(m2,q2,s);
    cdbl n = 1./(M_PI*M_PI);
    cdbl lbeta = log(beta);
    cdbl res = Color::CA * (a12 * lbeta*lbeta + a11 * lbeta - M_PI*M_PI/(16.*beta) + a10_OK(m2,q2))
                + 2.*Color::CF * (M_PI*M_PI/(16.*beta) + a10_QED(m2,q2));
    return lo * n * res;
}

cdbl FullyInclusive::coeffs::cg1t_F2_VV(cdbl m2, cdbl q2, cdbl s) {
    cdbl a12 = a12_;
    cdbl a11 = a11_;
    return cg1t_(m2,q2,s,&cg0t_F2_VV,a12,a11,&cg1_a10_F2_VV_OK,&cg1_a10_F2_VV_QED);
}

cdbl FullyInclusive::coeffs::cg1t_FL_VV(cdbl m2, cdbl q2, cdbl s) {
    cdbl a12 = a12_;
    cdbl a11 = a11_ - 2./3.;
    return cg1t_(m2,q2,s,&cg0t_FL_VV,a12,a11,&cg1_a10_FL_VV_OK,&cg1_a10_FL_VV_QED);
}

cdbl FullyInclusive::coeffs::cg1t_x2g1_VV(cdbl m2, cdbl q2, cdbl s) {
    cdbl a12 = a12_;
    cdbl a11 = a11_;
    return cg1t_(m2,q2,s,&cg0t_x2g1_VV,a12,a11,&cg1_a10_x2g1_VV_OK,&cg1_a10_x2g1_VV_QED);
}

cdbl FullyInclusive::coeffs::cg1t_F2_AA(cdbl m2, cdbl q2, cdbl s) {
    cdbl a12 = a12_;
    cdbl a11 = a11_;
    return cg1t_(m2,q2,s,&cg0t_F2_AA,a12,a11,&cg1_a10_F2_AA_OK,&cg1_a10_F2_AA_QED);
}

cdbl FullyInclusive::coeffs::cg1t_FL_AA(cdbl m2, cdbl q2, cdbl s) {
    cdbl a12 = a12_;
    cdbl a11 = a11_;
    return cg1t_(m2,q2,s,&cg0t_FL_AA,a12,a11,&cg1_a10_FL_AA_OK,&cg1_a10_FL_AA_QED);
}

cdbl FullyInclusive::coeffs::cg1t_x2g1_AA(cdbl m2, cdbl q2, cdbl s) {
    cdbl a12 = a12_;
    cdbl a11 = a11_;
    return cg1t_(m2,q2,s,&cg0t_x2g1_AA,a12,a11,&cg1_a10_x2g1_AA_OK,&cg1_a10_x2g1_AA_QED);
}

// improved threshold
cdbl FullyInclusive::coeffs::cg1tp_(cdbl m2, cdbl q2, cdbl s, fPtr3dbl cg1t, str proj_cc) {
    // true threshold
    cdbl t = cg1t(m2, q2, s);
    cdbl lnxi = log(-q2/m2);
    // load additional factor
    const boost::filesystem::path parent = Common::getPathByEnv("PARTONIC_GRIDS");
    const boost::filesystem::path mine = parent / "cg1" / ("cg1-"+proj_cc+"-thres-coeff.dat");
    if (!boost::filesystem::exists(mine))
        throw runtime_error("path does not exist: "+mine.string());
    gslpp::Interpolation1D a_int(gsl_interp_steffen,mine.string(),FullyInclusive::coeffs::cg1tp_lnxi_num);
    cdbl a = a_int.eval(lnxi);
    partonic_eta(m2, s)
    return t*(1. + a*eta);
}

// bulk
cdbl FullyInclusive::coeffs::cg1b_(cdbl m2, cdbl q2, cdbl s, str proj_cc) {
    // load additional factor
    const boost::filesystem::path parent = Common::getPathByEnv("PARTONIC_GRIDS");
    const boost::filesystem::path mine = parent / "cg1" / ("cg1-"+proj_cc+"-bulk.dat");
    if (!boost::filesystem::exists(mine))
        throw runtime_error("path does not exist: "+mine.string());
    gslpp::Interpolation2D bulk_int(gsl_interp2d_bicubic,mine.string(),FullyInclusive::coeffs::cg1b_lneta_num,FullyInclusive::coeffs::cg1b_lnxi_num);
    cdbl lnxi = log(-q2/m2);
    partonic_eta(m2, s)
    cdbl res = bulk_int.eval(log(eta),lnxi);
    return res;
}

// full
cdbl FullyInclusive::coeffs::cg1_(cdbl m2, cdbl q2, cdbl s, fPtr3dbl cg1tp, fPtr3dbl cg1b) {
    partonic_eta(m2, s)
    cdbl lneta = log(eta);
    cdbl low = FullyInclusive::coeffs::cg1b_lneta_min;
    // threshold only?
    if (lneta < low)
        return cg1tp(m2, q2, s);
    cdbl high = FullyInclusive::coeffs::cg1tp_lneta_max;
    // bulk only?
    if (lneta > high)
        return cg1b(m2, q2, s);
    // otherwise apply linear interpolation between the two
    cdbl tp = cg1tp(m2, q2, s);
    cdbl b = cg1b(m2, q2, s);
    return (tp*(lneta - high) + b*(low - lneta))/(low - high);
}

// keep all parity conserving elements
#define cg1_proj_pc(proj)\
cdbl FullyInclusive::coeffs::cg1tp_##proj(cdbl m2, cdbl q2, cdbl s) {\
    return cg1tp_(m2, q2, s, &cg1t_##proj, #proj);\
}\
cdbl FullyInclusive::coeffs::cg1b_##proj(cdbl m2, cdbl q2, cdbl s) {\
    return cg1b_(m2, q2, s, #proj);\
}\
cdbl FullyInclusive::coeffs::cg1_##proj(cdbl m2, cdbl q2, cdbl s) {\
    return cg1_(m2, q2, s, &cg1tp_##proj, &cg1b_##proj);\
}
interateAllPCProj(cg1_proj_pc)

// we assume all parity violating elements to be 0 (likely, but not proven)
#define cg1_proj_pv(proj) \
cdbl FullyInclusive::coeffs::cg1t_##proj(cdbl m2, cdbl q2, cdbl s) { return 0.; } \
cdbl FullyInclusive::coeffs::cg1tp_##proj(cdbl m2, cdbl q2, cdbl s) { return 0.; } \
cdbl FullyInclusive::coeffs::cg1b_##proj(cdbl m2, cdbl q2, cdbl s) { return 0.; } \
cdbl FullyInclusive::coeffs::cg1_##proj(cdbl m2, cdbl q2, cdbl s) { return 0.; }

interateAllPVProj(cg1_proj_pv)