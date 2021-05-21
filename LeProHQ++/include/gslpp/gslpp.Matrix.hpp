#ifndef gslpp_Matrix_HPP
#define gslpp_Matrix_HPP

#include <stdexcept>
#include <string>
#include <cstring>

#include <gsl/gsl_matrix.h>

/**
 * @brief object-oriented wrappers to GSL
 */
namespace gslpp {

/**
 * @brief wrapper to gsl_matrix
 */
class Matrix {
    
/**
 * @brief matrix object
 */
    gsl_matrix* m = 0;
    
/**
 * @brief first dimension
 */
    std::size_t n1 = 0;
    
/**
 * @brief second dimension
 */
    std::size_t n2 = 0;
    
public:

/**
 * @brief constructor - initializes dimensions
 * @param n1 first dimension
 * @param n2 second dimension
 */
    explicit Matrix(const std::size_t n1, const std::size_t n2) : n1(n1), n2(n2) {
        this->m = gsl_matrix_alloc(n1, n2);
    }
    
/**
 * @brief destructor
 */
    ~Matrix() {
        if(this->m)
            gsl_matrix_free(this->m);
    }
    
/**
 * @brief read data
 * @param stream target stream
 * @return gsl_matrix_fscanf
 */
    const int fscanf(FILE * stream) const {
        return gsl_matrix_fscanf(stream, this->m);
    }
    
/**
 * @brief read matrix from file path
 * @param path source path
 */
    void readFromFile(const std::string& path) const {
        FILE* f = fopen(path.c_str(),"r");
        if (f == NULL)
            throw ios::failure(std::strerror(errno));
        this->fscanf(f);
        fclose(f);
    }
    
/**
 * @brief access matrix element
 * @param i row
 * @param j column
 * @return element
 */
    double get(const size_t i, const size_t j) const {
        return gsl_matrix_get(this->m, i, j);
    }
    
/**
 * @brief returns first dimension
 * @return first dimension 
 */
    const size_t get_n1() const {
        return this->n1;
    }
    
/**
 * @brief returns second dimension
 * @return second dimension 
 */
    const size_t get_n2() const {
        return this->n2;
    }
    
/**
 * @brief obtain a row
 * @param i row number
 * @return gsl_vector_view row view
 */
    gsl_vector_view row(size_t i) const {
        return gsl_matrix_row(this->m, i);
    }
    
/**
 * @brief obtain a column
 * @param j column number
 * @return gsl_vector_view column view
 */
    gsl_vector_view column(size_t j) const {
        return gsl_matrix_column(this->m, j);
    }
    
/**
 * @brief obtain a subrow
 * @param i row number
 * @param offset offset
 * @param n length
 * @return gsl_vector_view row view
 */
    gsl_vector_view subrow(size_t i, size_t offset, size_t n) const {
        return gsl_matrix_subrow(this->m, i, offset, n);
    }
    
/**
 * @brief obtain a subcolumn
 * @param j column number
 * @param offset offset
 * @param n length
 * @return gsl_vector_view row view
 */
    gsl_vector_view subcolumn(size_t j, size_t offset, size_t n) const {
        return gsl_matrix_subcolumn(this->m, j, offset, n);
    }
};

} // namespace gsl

#endif // gslpp_Matrix_HPP