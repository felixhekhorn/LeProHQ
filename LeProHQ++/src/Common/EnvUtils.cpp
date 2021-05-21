#include "Common/EnvUtils.h"

#include <stdexcept>
#include <algorithm>

/**
 * @brief read environment variable
 * @see http://stackoverflow.com/questions/631664/accessing-environment-variables-in-c
 * @param key name
 * @return value - empty if not set
 */
const std::string Common::getEnvVar( std::string const & key ) {
    char * val = getenv( key.c_str() );
    return val == NULL ? std::string("") : std::string(val);
}

/**
 * @brief is environment variable set
 * @param key name
 * @return set?
 */
const bool Common::hasEnvVar( std::string const & key ) {
    char * val = getenv( key.c_str() );
    return !(val == NULL);
}

/**
 * @brief obtain a path from an environment variable
 * @param key name
 * @return path
 */
const boost::filesystem::path Common::getPathByEnv(const std::string& key) {
    if (!Common::hasEnvVar(key))
        throw std::invalid_argument("The environment variale "+key+" has to be set!");
    const std::string path = Common::getEnvVar(key);
    const boost::filesystem::path p(path);
    if (!boost::filesystem::exists(p) || !boost::filesystem::is_directory(p))
        throw std::invalid_argument("The environment variale "+key+" has to be set to an existing directory!");
    return p;
}
