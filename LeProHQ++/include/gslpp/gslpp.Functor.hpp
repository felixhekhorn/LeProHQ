#ifndef gslpp_Functor_HPP
#define gslpp_Functor_HPP

/**
 * @brief object-oriented wrappers to GSL
 */
namespace gslpp {

/**
 * @brief wrapper for gsl
 * @param x integration variable
 * @param params functor of class F
 * @return params(x)
 */
template <class F> inline double callFunctor (double x, void * params) {
    F* o = (F *) params;
    return o->operator()(x);
}

/**
 * @brief wrapper for gsl for dim = 1
 * @param x integration variables
 * @param dim dimension = 1
 * @param params functor of class F
 * @return params(x)
 */
template <class F> inline double callFunctor1D (double x[], size_t dim, void * params) {
    F* o = (F *) params;
    return o->operator()(x[0]);
}

/**
 * @brief wrapper for gsl for dim = 2
 * @param x integration variables
 * @param dim dimension = 2
 * @param params functor of class F
 * @return params(x,y)
 */
template <class F> inline double callFunctor2D (double x[], size_t dim, void * params) {
    F* o = (F *) params;
    return o->operator()(x[0],x[1]);
}

/**
 * @brief wrapper for gsl for dim = 3
 * @param x integration variables
 * @param dim dimension = 3
 * @param params functor of class F
 * @return params(x,y,z)
 */
template <class F> inline double callFunctor3D (double x[], size_t dim, void * params) {
    F* o = (F *) params;
    return o->operator()(x[0],x[1],x[2]);
}

/**
 * @brief wrapper for gsl for dim = 4
 * @param x integration variables
 * @param dim dimension = 4
 * @param params functor of class F
 * @return params(a1,a2,a3,a4)
 */
template <class F> inline double callFunctor4D (double x[], size_t dim, void * params) {
    F* o = (F *) params;
    return o->operator()(x[0],x[1],x[2],x[3]);
}

/**
 * @brief wrapper for gsl for dim = 5
 * @param x integration variables
 * @param dim dimension = 5
 * @param params functor of class F
 * @return params(a1,a2,a3,a4,a5)
 */
template <class F> inline double callFunctor5D (double x[], size_t dim, void * params) {
    F* o = (F *) params;
    return o->operator()(x[0],x[1],x[2],x[3],x[4]);
}

} // namespace gsl

#endif // gslpp_Functor_HPP