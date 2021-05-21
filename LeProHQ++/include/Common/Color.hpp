/**
 * @file Color.hpp
 * @brief defines color constants
 */
 
#ifndef Color_HPP_
#define Color_HPP_

#include "../config.h"
 
/** @brief [debug] flag to unset QED-like contributions */
#define _fQED 1.

/**
 * @brief constants related to color space
 */
namespace Color {
    
/**
 * @brief number of Colors \f$N_C \f$
 */
cuint NC = 3;

/**
 * @brief Casimir constant of fundamental representation \f$C_F=\frac{N_C^2 - 1}{2N_C} \f$
 */
cdbl CF = cdbl(NC*NC-1)/cdbl(2*NC);

/**
 * @brief Casimir constant of adjoint representation \f$C_A=N_C \f$
 */
cdbl CA = cdbl(NC);

/**
 * @brief Color avarage of \f$\Pg\f$ and \f$\Pgg\f$ in initial state \f$K_{\Pg\Pgg} = \frac 1 {N_C^2 -1} \f$
 */
cdbl Kgph = 1./(cdbl(NC*NC-1));

/**
 * @brief Color avarage of \f$\Pq\f$ and \f$\Pgg\f$ in initial state \f$K_{\Pq\Pgg} = \frac 1 {N_C} \f$
 */
cdbl Kqph = 1./(cdbl(NC));

}

#endif // Color_HPP_