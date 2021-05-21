#ifndef gslpp_Interpolation1D_HPP
#define gslpp_Interpolation1D_HPP

#include <stdexcept>
#include <string>
#include <cstring>
#include <unordered_map>

#include <gsl/gsl_spline.h>

#include "./gslpp.Matrix.hpp"

/**
 * @brief object-oriented wrappers to GSL
 */
namespace gslpp {

/**
 * @brief wrapper to gsl_spline
 */
class Interpolation1D {
    
/**
 * @brief spline object
 */
    gsl_spline* spline = 0;
    
/**
 * @brief number of points
 */
    size_t size;
    
/**
 * @brief accereleration object
 */
    gsl_interp_accel* acc = 0;
    
    void init_from_matrix(const gsl_interp_type * T, const gslpp::Matrix *m) {
        if (2 != m->get_n1())
            throw std::invalid_argument("Matrix has to be 2xn");
        this->size = m->get_n2();
        this->spline = gsl_spline_alloc (T, this->size);
        this->acc = gsl_interp_accel_alloc();
        gsl_spline_init(this->spline,m->row(0).vector.data,m->row(1).vector.data,this->size);
    }
    
public:

/**
 * @brief constructor with matrix
 * @param T interpolation type
 * @param m matrix with x and y points
 */
    explicit Interpolation1D(const gsl_interp_type * T, const gslpp::Matrix *m) {
        this->init_from_matrix(T, m);
    }
    
/**
 * @brief constructor with matrix from file
 * @param T interpolation type
 * @param path matrix path
 * @param n number of points
 */
    explicit Interpolation1D(const gsl_interp_type * T, const std::string& path, size_t n) {
        const gslpp::Matrix m(2,n);
        m.readFromFile(path);
        this->init_from_matrix(T, &m);
    }
    
    /*static std::unordered_map<std::string, Interpolation1D> cache;
    
    static Interpolation1D load(const gsl_interp_type * T, const std::string& path, size_t n, bool cache = true) {
        if (cache) {
//            Interpolation1D::cache;
        }
        Interpolation1D obj(T,path,n);
        if (cache)
            Interpolation1D::cache[path] = obj;
        return obj;
    }*/
    
/**
 * @brief destructor
 */
    ~Interpolation1D() {
        if(0 != this->spline)
            gsl_spline_free(this->spline);
        if (0 != this->acc)
            gsl_interp_accel_free(this->acc);
    }

/**
 * @brief evaluate interpolation
 * @param x interpolation point
 * @return interpolated value
 */
    double eval(const double x) const {
        return gsl_spline_eval(this->spline,x,this->acc);
    }
};

} // namespace gsl

#endif // gslpp_Interpolation1D_HPP