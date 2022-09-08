/**
 * @file PineAPPL.hpp
 * @brief object oriented interface to PineAPPL
 */
#ifndef PineAPPL_HPP_
#define PineAPPL_HPP_

#include <string>

#include <pineappl_capi.h>

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
    void set(const std::string& key, const double value) { pineappl_keyval_set_double(this->raw, key.c_str(), value); }
    void set(const std::string& key, const bool value) { pineappl_keyval_set_bool(this->raw, key.c_str(), value); }
    void set(const std::string& key, const int value) { pineappl_keyval_set_int(this->raw, key.c_str(), value); }
    void set(const std::string& key, const std::string& value) { pineappl_keyval_set_string(this->raw, key.c_str(), value.c_str()); }
    ///@}

    /** @name getter */
    ///@{
    double get_double(const std::string& key) { return pineappl_keyval_double(this->raw, key.c_str()); }
    bool get_bool(const std::string& key) { return pineappl_keyval_bool(this->raw, key.c_str()); }
    int get_int(const std::string& key) { return pineappl_keyval_int(this->raw, key.c_str()); }
    std::string get_string(const std::string& key) { return pineappl_keyval_string(this->raw, key.c_str()); }
    ///@}

};

}


#endif // PineAPPL_HPP_
