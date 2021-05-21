/**
 * @file Projection.hpp
 * @brief defines projections and related functions
 */
#ifndef Projection_HPP_
#define Projection_HPP_

#include <string>

/**
 * @brief on what structure function to project on
 */
enum Projection {
    F2,  /**< unpolarized, parity conserving F_2(x,Q2) */
    FL,  /**< unpolarized, parity conserving F_L(x,Q2) = F_2(x,Q2) - 2xF_1(x,Q2) */
    xF3, /**< unpolarized, parity violating x*F_3(x,Q2) */
    x2g1, /**< polarized, parity conserving 2x*g_1(x,Q2) (Note the naming limitations of c++) */
    g4,   /**< polarized, parity violating g_4(x,Q2) */
    gL    /**< polarized, parity violating g_L(x,Q2) = g_4(x,Q2) - 2xg_5(x,Q2) */
};

/**
 * @brief does projection conserve parity?
 * @param proj projection
 * @return projection conserves parity?
 */
inline bool isParityConservingProj(Projection proj) { return proj == F2 || proj == FL || proj == x2g1; };

/**
 * @brief cast projection to string
 * @param proj projection
 * @return (string) projection
 */
inline std::string strProj(Projection proj) {
    if (F2 == proj) return "F2";
    if (FL == proj) return "FL";
    if (xF3 == proj) return "xF3";
    if (x2g1 == proj) return "x2g1";
    if (g4 == proj) return "g4";
    if (gL == proj) return "gL";
    return "unknown projection!";
};

#endif // Projection_HPP_