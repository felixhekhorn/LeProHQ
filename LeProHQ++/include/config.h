/**
 * @file config.h
 * @brief includes all common header files and defines types and shortcuts
 */

#ifndef Config_H_
#define Config_H_

#include <stdlib.h>
#include <cmath>
#include <string>
#include <stdexcept>
#include <iostream>
#include <stdio.h>

using namespace std;

/** @name shorthands */
///@{
/** @brief unsigned int shorthand */
typedef unsigned int uint;
/** @brief const unsigned int shorthand */
typedef const uint cuint;
/** @brief const int shorthand */
typedef const int cint;

/** @brief string shorthand */
typedef string str;

/** @brief defines floating point precision */
typedef double dbl;
/** @brief defines floating point precision */
typedef const dbl cdbl;

/** @brief NaN (often (mis-)used by me as synonym to null-pointer/void) */
cdbl dblNaN = nan("");
///@}

/**
 * @brief iterates all projections by macro
 * @param MEHfnc macro
 */
#define interateAllProj(MEHfnc) interateAllPCProj(MEHfnc) interateAllPVProj(MEHfnc)
    
/**
 * @brief iterates all parity conserving projections by macro
 * @param MEHfnc macro
 */
#define interateAllPCProj(MEHfnc) \
    MEHfnc(F2_VV) MEHfnc(F2_AA) \
    MEHfnc(FL_VV) MEHfnc(FL_AA) \
    MEHfnc(x2g1_VV) MEHfnc(x2g1_AA)
    
/**
 * @brief iterates all parity violating projections by macro
 * @param MEHfnc macro
 */
#define interateAllPVProj(MEHfnc) MEHfnc(xF3_VA) MEHfnc(g4_VA) MEHfnc(gL_VA)
    
#endif // Config_H_