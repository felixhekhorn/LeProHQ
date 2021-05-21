#ifndef gslpp_Interpolation2D_HPP
#define gslpp_Interpolation2D_HPP

#include <stdexcept>
#include <string>
#include <cstring>
#include <unordered_map>

#include <gsl/gsl_spline2d.h>

#include "./gslpp.Matrix.hpp"

/**
 * @brief object-oriented wrappers to GSL
 */
namespace gslpp {

/**
 * @brief wrapper to gsl_spline
 */
class Interpolation2D {
    
/**
 * @brief spline object
 */
    gsl_spline2d* spline = 0;
    
/**
 * @brief number of x-points
 */
    size_t nx;
    
/**
 * @brief number of y-points
 */
    size_t ny;
    
/**
 * @brief x-accereleration object
 */
    gsl_interp_accel* xacc = 0;
    
/**
 * @brief y-accereleration object
 */
    gsl_interp_accel* yacc = 0;
    
    double *xs = 0;
    double *ys = 0;
    double *zs = 0;
    
/**
 * @brief constructor delegation
 * @param T interpolation type
 * @param m matrix with x/y points as first column/row and z values in the bulk
 */
    void init_from_matrix(const gsl_interp2d_type * T, const gslpp::Matrix *m) {
        if (m->get_n1() < 2 || m->get_n2() < 2)
            throw std::invalid_argument("Matrix has to be bigger then 2x2");
        this->nx = m->get_n1() - 1;
        this->ny = m->get_n2() - 1;
        // alloc everything
        this->spline = gsl_spline2d_alloc (T, this->nx, this->ny);
        this->xacc = gsl_interp_accel_alloc();
        this->yacc = gsl_interp_accel_alloc();
        this->xs = (double*)malloc(nx  * sizeof(double));
        this->ys = (double*)malloc(ny * sizeof(double));
        this->zs = (double*)malloc(nx * ny * sizeof(double));
        // get matrix border
        gsl_vector vxs = m->subcolumn(0,1,this->nx).vector;
        gsl_vector vys = m->subrow(0,1,this->ny).vector;
        // assign everything
        for (unsigned int j = 0; j < this->nx; ++j) {
            for (unsigned int k = 0; k < this->ny; ++k) {
                if (0 == j) {
                    this->ys[k] = gsl_vector_get(&vys,k);
                }
                gsl_spline2d_set(this->spline, this->zs, j, k, m->get(1+j,1+k));
            }
            this->xs[j] = gsl_vector_get(&vxs,j);
        }
        // finally init
        gsl_spline2d_init(this->spline,this->xs,this->ys, this->zs, this->nx, this->ny);
    }
    
public:

/**
 * @brief constructor with matrix
 * @param T interpolation type
 * @param m matrix with x/y points as first column/row and z values in the bulk
 */
    explicit Interpolation2D(const gsl_interp2d_type * T, const gslpp::Matrix *m) {
        this->init_from_matrix(T, m);
    }
    
/**
 * @brief 
 * @param T interpolation type
 * @param path matrix path
 * @param nx x size
 * @param ny y size
 */
    explicit Interpolation2D(const gsl_interp2d_type * T, const std::string& path, size_t nx, size_t ny) {
        const gslpp::Matrix m (1+nx,1+ny);
        m.readFromFile(path);
        this->init_from_matrix(T, &m);
    }
    
    
/**
 * @brief destructor
 */
    ~Interpolation2D() {
        if(0 != this->spline)
            gsl_spline2d_free(this->spline);
        if (0 != this->xacc)
            gsl_interp_accel_free(this->xacc);
        if (0 != this->yacc)
            gsl_interp_accel_free(this->yacc);
        if (0 != this->xs)
            free(this->xs);
        if (0 != this->ys)
            free(this->ys);
        if (0 != this->zs)
            free(this->zs);
    }

/**
 * @brief evaluate interpolation
 * @param x interpolation point on x-axis
 * @param y interpolation point on y-axis
 * @return interpolated value
 */
    double eval(const double x, const double y) const {
        return gsl_spline2d_eval(this->spline,x,y,this->xacc,this->yacc);
    }
};

} // namespace gsl

#endif // gslpp_Interpolation2D_HPP