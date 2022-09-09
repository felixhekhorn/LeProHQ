/**
 * @file PineAPPL.hpp
 * @brief object oriented interface to PineAPPL
 */
#ifndef PineAPPL_HPP_
#define PineAPPL_HPP_

#include <string>
#include <list>

#include <pineappl_capi.h>

/** @brief object oriented interface to PineAPPL  */
namespace PineAPPL {

/** @brief Key-value storage for passing optional information during grid creation */
struct KeyVal {

    /** @brief underlying raw object */ 
    pineappl_keyval *raw;

    /** @brief constructor */
    KeyVal() { this->raw = pineappl_keyval_new(); }

    /** @brief destructor */
    ~KeyVal() { pineappl_keyval_delete(this->raw); }

    
    /** @name setter */
    ///@{
    void set_double(const std::string& key, const double value) const { pineappl_keyval_set_double(this->raw, key.c_str(), value); }
    void set_bool(const std::string& key, const bool value) const { pineappl_keyval_set_bool(this->raw, key.c_str(), value); }
    void set_int(const std::string& key, const int value) const { pineappl_keyval_set_int(this->raw, key.c_str(), value); }
    void set_string(const std::string& key, const std::string& value) const { pineappl_keyval_set_string(this->raw, key.c_str(), value.c_str()); }
    ///@}

    /** @name getter */
    ///@{
    double get_double(const std::string& key) const { return pineappl_keyval_double(this->raw, key.c_str()); }
    bool get_bool(const std::string& key) const { return pineappl_keyval_bool(this->raw, key.c_str()); }
    int get_int(const std::string& key) const { return pineappl_keyval_int(this->raw, key.c_str()); }
    std::string get_string(const std::string& key) const { return pineappl_keyval_string(this->raw, key.c_str()); }
    ///@}

};

/** @brief Entry in luminosity function  */
struct LumiEntry {
    int32_t pid1;
    int32_t pid2;
    double weight;
};

/** @brief Luminosity function */
struct Lumi {

    /** @brief underlying raw object */
    pineappl_lumi *raw;

    /** @brief constructor */
    Lumi(){ this->raw = pineappl_lumi_new(); }

    /** @brief destructor */
    ~Lumi(){ pineappl_lumi_delete(this->raw); }

    size_t count() const { return pineappl_lumi_count(this->raw); }

    /**
     * @brief add a luminosity function
     * @param c luminosity function
     */
    void add(const list<LumiEntry>& c) const {
        auto itr = c.cbegin();
        auto end = c.cend();
        size_t n = 0;
        int32_t *pids = new int32_t;
        double *weights = new double;
        for (; itr != end; ++itr, ++n) {
            cout << (*itr).pid1 << "," << itr->pid2 << " " << itr->weight <<endl;
            pids[2*n] = itr->pid1;
            pids[2*n+1] = itr->pid2;
            weights[n] = itr->weight;
        }
        pineappl_lumi_add(this->raw, n, pids, weights);
        delete pids;
        delete weights;
    }

    /**
     * @brief Returns the number of combinations of the luminosity function `lumi` for the specified entry.
     * @param entry position in lumi
     */
    size_t combinations(size_t entry) const { return pineappl_lumi_combinations(this->raw, entry); }
};

/** @brief Coupling powers for each grid. */
struct Order {
    /** @brief Exponent of the strong coupling. */
    uint32_t alphas;
    /** @brief Exponent of the electromagnetic coupling. */
    uint32_t alpha;
    /** @brief Exponent of the logarithm of the scale factor of the renomalization scale. */
    uint32_t logxir;
    /** @brief Exponent of the logarithm of the scale factor of the factorization scale. */
    uint32_t logxif;
};

}


#endif // PineAPPL_HPP_
