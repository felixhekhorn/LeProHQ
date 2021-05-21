/**
 * @brief defines interpolation routines
 */

#ifndef Interpolation_H_
#define Interpolation_H_

#include <list>
#include <vector>

#include "config.h"
#include "MellinVars.hpp"
#include "Pdf/PdfWrapper.h"

namespace FullyInclusiveSoftResummed {

/**
 * @brief namespace for interpolation functions
 */
namespace Interpolation {

/** @brief list of doubles */    
typedef list<dbl> lsdbl;

/** @brief vector of doubles */
typedef vector<dbl> vdbl;


/**
 * @brief Helper function to eliminate close double numbers
 * @param a first operand
 * @param b second operand
 * @return |a-b| < 1e-10?
 */
bool is_close(cdbl a, cdbl b);

/**
 * @brief Creates a log-lin-spaced grid
 * 
 * 1.0 is always part of the grid and the final grid is unified, i.e. esp. points that might
 * appear twice at the borders of the regions are unified.
 * @param n_low points in the small-x region
 * @param n_mid points in the mediun-x region
 * @param x_min minimum x (included)
 * @param x_low seperation point between small and medium
 * @return grid
 */
const lsdbl make_log_lin_grid(cuint n_low, cuint n_mid, cdbl x_min=1e-7, dbl x_low=0.1);

/**
 * @brief factorial function
 * @param n argument
 * @return n!
 */
uint fact(uint n);

/**
 * @brief Provides an interpolation interface.
 */
class Interpolator {
    
/** @brief grid in x-space from which the interpolators are constructed */
    vdbl xgrid;
    
/** @brief degree of the interpolation polynomial */
    uint polynomialDegree;
    
/** @brief block borders type */
    typedef pair<uint,uint> blockT;
    
/** @brief list of blocks type */
    typedef list<blockT> lsblockT;
    
/** @brief basis function configuration type */
    typedef list<lsdbl> bfT;
    
/** @brief full interpolation configuration type */
    typedef vector<bfT> vconfigT;
    
/** @brief interpolation configuration */
    vconfigT config;
    
/** @brief logarithmic interpolation? */
    bool isLog;
    
/** reference PDF values */
    lsdbl data;

/**
 * @brief Setup interpolation blocks
 * @param xgridSize grid length
 * @param polynomialDegree degree of the interpolation polynomial
 * @return block configurations
 */
    lsblockT setupBlocks(cuint xgridSize, cuint polynomialDegree) const;

/**
 * @brief Setup a single basis function
 * @param ygrid evtually transformed grid
 * @param polynomNumber associated polynomial number
 * @param blocks precomputed list of blocks
 * @return area configurations
 */    
    bfT setupBasisFunction(const vdbl ygrid, cuint polynomNumber, const lsblockT& blocks) const;
    
/**
 * @brief Setup a single area for a single basis function
 * @param ygrid evtually transformed grid
 * @param polynomNumber associated polynomial number
 * @param lowerIndex index of the left border of the area
 * @param block block borders
 * @return list with {xmin, xmax, coefficients}
 */
    lsdbl setupArea(const vdbl ygrid, cuint polynomNumber, cuint lowerIndex, const blockT& block) const;
    
/**
 * @brief Evaluates a single polynomial in x-space at a given point in x
 * @param bf_config polynomial config
 * @param y x or log(x)
 * @param isFirst is first polynomial?
 * @return \f$p_j(x)\f$
 */
    cdbl evalPolyX(const bfT bf_config, cdbl y, bool isFirst = false) const;
    
/**
 * @brief Evaluates a single logarithmic polynomial and the inversion factor in N-space at a given point
 * @param bf_config polynomial config
 * @param N evaluation point
 * @param logx logarithm of inversion point
 * @return \f$exp(- N \log(x)) \tilde p_j(N)\f$
 */
    FullyInclusiveSoftResummed::cdcmplx evalLogPolyN(const bfT bf_config, FullyInclusiveSoftResummed::cdcmplx N, cdbl logx) const;
    
/**
 * @brief Evaluates a single linear polynomial and the inversion factor in N-space at a given point
 * @param bf_config polynomial config
 * @param N evaluation point
 * @param logx logarithm of inversion point
 * @return \f$exp(- N \log(x)) \tilde p_j(N)\f$
 */
    FullyInclusiveSoftResummed::cdcmplx evalLinPolyN(const bfT bf_config, FullyInclusiveSoftResummed::cdcmplx N, cdbl logx) const;
    
public:

/**
 * @brief Setup interpolation handler
 * 
 * Interpolation is performed in x-space, but the generated basis polynomials can be evaluated (exactly)
 * also in N-space.
 * 
 * @param xgrid grid in x-space from which the interpolation is constructed
 * @param polynomialDegree degree of the interpolation polynomial
 * @param isLog logarithmic interpolation?
 */
    explicit Interpolator(const vdbl xgrid, cuint polynomialDegree, bool isLog = true);
    
/**
 * @brief Evaluates a single polynomial in x-space at a given point in x
 * @param j polynomial number
 * @param x evaluation point
 * @return \f$p_j(x)\f$
 */
    cdbl evaluatePolynomX(cuint j, cdbl x) const;
    
/**
 * @brief Evaluates a single polynomial and the inversion factor in N-space at a given point
 * @param j polynomial number
 * @param N evaluation point
 * @param logx logarithm of inversion point
 * @return \f$exp(- N \log(x)) \tilde p_j(N)\f$
 */
    FullyInclusiveSoftResummed::cdcmplx evaluatePolynomN(cuint j, FullyInclusiveSoftResummed::cdcmplx N, cdbl logx) const;
    
/**
 * @brief Load reference PDF values
 * @param pdf pdf object
 * @param pid parton identifier
 * @param muF2 factorization scale
 */
    void loadPdfData(PdfWrapper* pdf, const int pid, cdbl muF2);
    
/**
 * @brief Evaluates the full function in x-space at a given point
 * @param x evaluation point
 * @return \f$f(x)\f$
 */
    cdbl evaluateX(cdbl x) const;
    
/**
 * @brief Evaluates the full function and the inversion factor in N-space at a given point
 * @param Nevaluation point
 * @param logx logarithm of inversion point
 * @return \f$exp(- N \log(x)) \tilde f(N)\f$
 */
    FullyInclusiveSoftResummed::cdcmplx evaluateN(FullyInclusiveSoftResummed::cdcmplx N, cdbl logx) const;
};

} // namespace Interpolation

} // namespace Common

#endif // Interpolation_H_