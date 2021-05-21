#ifndef FullyDiff_PhasespaceValues_H
#define FullyDiff_PhasespaceValues_H

#include "config.h"

namespace FullyDiff {

/**
 * @brief holds the values of all kernels
 */
struct PhasespaceValues {
    
/**
 * @brief kernel at event x and event y
 */
    dbl xEyE = 0.;
    
/**
 * @brief kernel at counter event x and event y
 */
    dbl xCyE = 0.;
    
/**
 * @brief kernel at event x and counter event y
 */
    dbl xEyC = 0.;
    
/**
 * @brief kernel at counter event x and counter event y
 */
    dbl xCyC = 0.;
    
/**
 * @brief total weight
 * @return sum of all parts
 */
    inline cdbl tot() const {
        return this->xEyE + this->xCyE + this->xEyC + this->xCyC;
    }
};

} // namespace FullyDiff

#endif // FullyDiff_PhasespaceValues_H