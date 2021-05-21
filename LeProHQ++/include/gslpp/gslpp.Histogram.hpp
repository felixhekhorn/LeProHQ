#ifndef gslpp_Histogram_HPP
#define gslpp_Histogram_HPP

#include <stdexcept>
#include <string>

#include <gsl/gsl_histogram.h>

/**
 * @brief object-oriented wrappers to GSL
 */
namespace gslpp {

/**
 * @brief wrapper to gsl_histogram
 */
class Histogram {
    
/**
 * @brief number of bins
 */
    size_t size;
    
/**
 * @brief histogram object
 */
    gsl_histogram* h = 0;
    
/**
 * @brief are ranges set?
 */
    bool initialized = false;
    
/**
 * @brief file path
 */
    std::string path;
    
public:

/**
 * @brief constructor - initializes histogram
 * @param size number of bins
 */
    explicit Histogram(const std::size_t size) : size(size) {
        this->h = gsl_histogram_alloc(size);
    }
    
/**
 * @brief destructor
 */
    ~Histogram() {
        if(0 != this->h)
            gsl_histogram_free(this->h);
    }

/**
 * @brief sets uniform ranges
 * @param min minimun value
 * @param max maximum value
 * @return gsl_histogram_set_ranges_uniform
 */
    const int setRangesUniform(const double min, const double max) {
        if (max <= min)
            throw std::invalid_argument("expects min < max!");
        const int r = gsl_histogram_set_ranges_uniform(this->h,min,max);
        this->initialized = true;
        return r;
    }

/**
 * @brief sets logarithmic ranges
 * @param min minimun value
 * @param max maximum value
 * @return gsl_histogram_set_ranges
 */
    const int setRangesLog10(const double min, const double max) {
        if (min <= 0. || max <= 0. || max <= min)
            throw std::invalid_argument("expects 0 < min < max!");
        const double mi = log10(min);
        const double ma = log10(max);
        const double d = (ma - mi)/this->size;
        double x = mi;
        size_t s = this->size + 1;
        double range[s];
        for (size_t j = 0; j < s; ++j) {
            range[j] = pow(10.,x);
            x += d;
        }
        const int ret = gsl_histogram_set_ranges(this->h,range,s);
        this->initialized = true;
        return ret;
    }
    
/**
 * @brief is histogram initialized? (i.e. are ranges set?)
 * @return initialized?
 */
    const bool isInitialized() const {
        return this->initialized;
    }
    
/**
 * @brief sets file path
 * @param path file path
 */
    void setPath(const std::string& path) {
        this->path = path;
    }
    
/**
 * @brief accumulates by 1.0
 * @param x value
 * @return gsl_histogram_accumulate
 */
    const int increment(const double x) const {
        return this->accumulate(x, 1.);
    }
    
/**
 * @brief adds w to the bin containing x
 * @param x value
 * @param w weight
 * @return gsl_histogram_accumulate
 */
    const int accumulate(const double x, const double w) const {
        return gsl_histogram_accumulate(this->h, x, w);
    }
    
/**
 * @brief scales all histogram bins
 * @param s scale
 * @return gsl_histogram_scale
 */
    const int scale(const double s) const {
        return gsl_histogram_scale(this->h,s);
    }

/**
 * @brief outputs data
 * @param stream target stream
 * @param range_format printf-format for bin ranges
 * @param bin_format printf-format for bin value
 * @return gsl_histogram_fprintf
 */
    const int fprintf(FILE * stream, const char * range_format, const char * bin_format) const {
        // empty?
        if (!this->isInitialized())
            return 0;
        return gsl_histogram_fprintf(stream, this->h, range_format, bin_format);
    }
    
/**
 * @brief writes histogram to the file path
 */
    void writeToFile() const {
        FILE* f = fopen(this->path.c_str(),"w");
        if (f == NULL)
            throw ios::failure(strerror(errno));
        this->fprintf(f, "% e", "% e");
        fclose(f);
    }
    
};

} // namespace gsl

#endif // gslpp_Histogram_HPP