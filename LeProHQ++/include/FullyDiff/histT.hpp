#ifndef histT_H_
#define histT_H_

#include <unordered_map>

#include "../gslpp/gslpp.Histogram.hpp"

namespace FullyDiff {
    
/**
 * @brief histogram types
 */
enum histT {
// All order:
    log10z, /**< bjorkenX < z < z_max */
    log10xi, /**< bjorkenX/z_max < xi < 1. */
    Theta1, /**< 0 < Theta1 < pi */
    
    HQPairInvMass, /**< 2m < invMassHQPair < sqrt(S) */
    HQPairDeltaPhi, /**< -pi < HQPairDeltaPhi < pi */
    HQPairTransverseMomentum, /**< 0 < HQPairTransverseMomentum < 2sqrt(S/4-m2) @todo: upper limit is actually smaller  */
    HQPairConeSizeVariable, /**< 0 < HQPairConeSizeVariable @todo: determine upper limit */
    
    HAQRapidity, /**< -y0 < HAQRapidity < y0 where y0 = artanh(sqrt(1-4m2/S)) */
    HAQTransverseMomentum, /**< 0 < HAQTransverseMomentum < sqrt(S/4-m2) */
    HAQTransverseMomentumScaling, /**< 0 < HAQTransverseMomentumScaling < 1 */
    HAQFeynmanX, /**< -1 < HAQFeynmanX < 1 */
    
// NLO:
    x, /**< bjorkenX/z_max < x < 1 */
    y, /**< -1 < y < 1 */
    Theta2 /**< 0 < Theta2 < pi */
};

/**
 * @brief shorthand for map of histograms
 */
typedef unordered_map<histT,gslpp::Histogram*> histMapT;

} // namespace FullyDiff

#endif // histT_H_