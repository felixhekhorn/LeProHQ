#include "FullyInclusiveSoftResummed/Interpolation.h"

using namespace FullyInclusiveSoftResummed::Interpolation;

bool FullyInclusiveSoftResummed::Interpolation::is_close(cdbl a, cdbl b) { return abs(a-b) < 1e-10; }    

const lsdbl FullyInclusiveSoftResummed::Interpolation::make_log_lin_grid(cuint n_low, cuint n_mid, cdbl x_min, dbl x_low) {
    if (n_mid == 0)
        x_low = 1.0;
    lsdbl l;
    // log
    if (n_low > 0) {
        cdbl dlogx = exp(log(x_low / x_min) / (n_low - 1));
        dbl x = x_min;
        for (uint j = 0; j < n_low; ++j, x *= dlogx)
            l.push_back(x);
    }
    // lin
    if (n_mid > 0) {
        cdbl dx = (1.0 - x_low) / (n_mid - 1);
        dbl x = x_low;
        for (uint j = 0; j < n_mid; ++j, x += dx)
            l.push_back(x);
    }
    l.push_back(1.0);
    l.sort();
    l.unique(is_close);
    return l;
}

uint FullyInclusiveSoftResummed::Interpolation::fact(uint n){ return (n==0) || (n==1) ? 1 : n* fact(n-1); }

Interpolator::lsblockT FullyInclusiveSoftResummed::Interpolation::Interpolator::setupBlocks(cuint xgridSize, cuint polynomialDegree) const {
    // Create blocks
    lsblockT blocks;
    uint po2 = polynomialDegree / 2;
    // if degree is even, we can not split the block symmetric, e.g. deg=2 -> |-|-|
    // so, in case of doubt use the block, which lays higher, i.e.
    // we're not allowed to go so deep -> make po2 smaller
    if (polynomialDegree % 2 == 0)
        po2 -= 1;
    // iterate areas: there is 1 less then number of points
    for (uint j = 0; j < xgridSize - 1; ++j) {
        uint kmin = max(0, (int)(j - po2));
        uint kmax = kmin + polynomialDegree;
        if (kmax >= xgridSize) {
            kmax = xgridSize - 1;
            kmin = kmax - polynomialDegree;
        }
        blocks.emplace_back(kmin, kmax);
    }
    return blocks;
}

Interpolator::bfT FullyInclusiveSoftResummed::Interpolation::Interpolator::setupBasisFunction(const vdbl ygrid, cuint polynomNumber, const lsblockT& blocks) const {
    bfT areas;
    uint j = 0;
    for (auto it = blocks.cbegin(); it != blocks.cend(); ++it, ++j) {
        if (it->first <= polynomNumber && polynomNumber <= it->second)
            areas.push_back(this->setupArea(ygrid, polynomNumber, j, *it));
    }
    return areas;
}

lsdbl FullyInclusiveSoftResummed::Interpolation::Interpolator::setupArea(const vdbl ygrid, cuint polynomNumber, cuint lowerIndex, const blockT& block) const {
    cdbl xmin = ygrid[lowerIndex];
    cdbl xmax = ygrid[lowerIndex + 1];
    dbl denominator = 1.0;
    lsdbl coeffs = {1.0};
    cdbl xj = ygrid[polynomNumber];
    for (uint k = block.first; k <= block.second; ++k) {
        if (k == polynomNumber)
            continue;
        cdbl xk = ygrid[k];
        denominator *= xj - xk;
        // Mxk_coeffs = -xk * coeffs
        lsdbl Mxk_coeffs(coeffs);
        for (auto& e : Mxk_coeffs)
            e *= -xk;
        coeffs.push_front(0.0);
        // coeffs[: s + +1] += Mxk_coeffs
        auto coeffs_it = coeffs.begin();
        for (auto it = Mxk_coeffs.begin(); it != Mxk_coeffs.end(); ++it, ++coeffs_it)
            *coeffs_it += *it;
    }
    for (auto& e : coeffs)
        e /= denominator;
    // add range to the front
    coeffs.push_front(xmax);
    coeffs.push_front(xmin);
    return coeffs;
}

cdbl FullyInclusiveSoftResummed::Interpolation::Interpolator::evalPolyX(const bfT bf_config, cdbl y, bool isFirst) const {
    dbl res = 0.0;
    for (const auto areas_it : bf_config) {
        auto coeffs_it = areas_it.cbegin();
        cdbl xmin = *coeffs_it;
        coeffs_it++;
        cdbl xmax = *coeffs_it;
        coeffs_it++;
        if ((xmin < y and y <= xmax) or (isFirst and (y == xmin))) {
            dbl xx = 1.;
            for (; coeffs_it != areas_it.cend(); ++coeffs_it, xx *= y)
                res += (*coeffs_it) * xx;
        }
    }
    return res;
}

FullyInclusiveSoftResummed::cdcmplx FullyInclusiveSoftResummed::Interpolation::Interpolator::evalLogPolyN(const bfT bf_config, FullyInclusiveSoftResummed::cdcmplx N, cdbl logx) const {
    dcmplx res = 0.0;
    for (const auto areas_it : bf_config) {
        auto coeffs_it = areas_it.cbegin();
        cdbl lnxmin = *coeffs_it;
        coeffs_it++;
        cdbl lnxmax = *coeffs_it;
        coeffs_it++;
        // skip area completely?
        if (logx >= lnxmax)
            continue;
        cdcmplx umax = N * lnxmax;
        cdcmplx umin = N * lnxmin;
        cdcmplx emax = exp(N * (lnxmax - logx));
        cdcmplx emin = exp(N * (lnxmin - logx));
        uint jj = 0;
        for (; coeffs_it != areas_it.cend(); ++coeffs_it, ++jj) {
            dcmplx tmp = 0.0;
            cdcmplx factjj = fact(jj) * pow(-1, jj) / pow(N, jj+1);
            for (uint k = 0; k <= jj; ++k) {
                cdbl factk = 1.0 / fact(k);
                dcmplx pmax = emax;
                if (!(is_close(umax.real(), 0.0) && is_close(umax.imag(), 0.0) && k == 0))
                    pmax *= pow(-umax, k);
                dcmplx pmin = 0.0;
                // drop factor by analytics?
                if (logx < lnxmin)
                    pmin = pow(-umin, k) * emin;
                tmp += factk * (pmax - pmin);
            }
            res += (*coeffs_it) * factjj * tmp;
        }
    }
    return res;
}
    
FullyInclusiveSoftResummed::cdcmplx FullyInclusiveSoftResummed::Interpolation::Interpolator::evalLinPolyN(const bfT bf_config, FullyInclusiveSoftResummed::cdcmplx N, cdbl logx) const {
    dcmplx res = 0.0;
    for (const auto areas_it : bf_config) {
        auto coeffs_it = areas_it.cbegin();
        cdbl xmin = *coeffs_it;
        coeffs_it++;
        cdbl xmax = *coeffs_it;
        coeffs_it++;
        cdbl lnxmax = log(xmax);
        // skip area completely?
        if (logx >= lnxmax)
            continue;
        uint k = 0;
        for (; coeffs_it != areas_it.cend(); ++coeffs_it, ++k) {
            dcmplx low = 0.0;
            if (!is_close(xmin, 0.0)) {
                cdbl lnxmin = log(xmin);
                low = exp(N * (lnxmin - logx) + k*lnxmin);
            }
            cdcmplx up = exp(N * (lnxmax - logx) + k*lnxmax);
            res += (*coeffs_it) * (up - low) / (N + dbl(k));
        }
    }
    return res;
}

FullyInclusiveSoftResummed::Interpolation::Interpolator::Interpolator(const vdbl xgrid, cuint polynomialDegree, bool isLog) :
            xgrid(xgrid), polynomialDegree(polynomialDegree), config(), isLog(isLog) {
    // Create blocks
    cuint xgridSize = xgrid.size();
    const lsblockT blocks = this->setupBlocks(xgridSize, polynomialDegree);
    // apply logarithm?
    vdbl ygrid(xgrid);
    if (isLog) {
        for (auto& y : ygrid)
            y = log(y);
    }
    // setup a basis function for each point
    for (uint j = 0; j < xgridSize; ++j)
        this->config.push_back(this->setupBasisFunction(ygrid, j, blocks));
}

cdbl FullyInclusiveSoftResummed::Interpolation::Interpolator::evaluatePolynomX(cuint j, cdbl x) const {
    dbl y = x;
    if (this->isLog)
        y = log(y);
    return this->evalPolyX(this->config[j], y, j == 0);
}

FullyInclusiveSoftResummed::cdcmplx FullyInclusiveSoftResummed::Interpolation::Interpolator::evaluatePolynomN(cuint j, FullyInclusiveSoftResummed::cdcmplx N, cdbl logx) const {
    const bfT bf_config = this->config[j];
    if (this->isLog)
        return this->evalLogPolyN(bf_config, N, logx);
    return this->evalLinPolyN(bf_config, N, logx);
}

void FullyInclusiveSoftResummed::Interpolation::Interpolator::loadPdfData(PdfWrapper* pdf, const int pid, cdbl muF2) {
    this->data.clear();
    for (const auto x : this->xgrid) {
        this->data.push_back(pdf->xfxQ2(pid, x, muF2));
    }
}

cdbl FullyInclusiveSoftResummed::Interpolation::Interpolator::evaluateX(cdbl x) const {
    if (this->data.size() != this->config.size())
        throw invalid_argument("reference data and grid do not have the same size!");
    dbl y = x;
    if (this->isLog)
        y = log(y);
    dbl res = 0.0;
    auto ref = this->data.cbegin();
    bool isFirst = true;
    for (const auto pj : this->config) {
        res += ((*ref) * this->evalPolyX(pj, y, isFirst));
        ++ref;
        isFirst = false;
    }
    return res;
}

FullyInclusiveSoftResummed::cdcmplx FullyInclusiveSoftResummed::Interpolation::Interpolator::evaluateN(FullyInclusiveSoftResummed::cdcmplx N, cdbl logx) const {
    if (this->data.size() != this->config.size())
        throw invalid_argument("reference data and grid do not have the same size!");
    dcmplx res = 0.0;
    auto ref = this->data.cbegin();
    for (const auto pj : this->config) {
        if (this->isLog)
            res += ((*ref) * this->evalLogPolyN(pj, N, logx));
        else
            res += ((*ref) * this->evalLinPolyN(pj, N, logx));
        ++ref;
    }
    return res;
}