/**
 * @brief defines environment variables access
 */

#ifndef EnvUtils_H_
#define EnvUtils_H_

#include <boost/filesystem.hpp>

namespace Common {

/**
 * @brief read environment variable
 * @see http://stackoverflow.com/questions/631664/accessing-environment-variables-in-c
 * @param key name
 * @return value - empty if not set
 */
const std::string getEnvVar( std::string const & key );

/**
 * @brief is environment variable set
 * @param key name
 * @return set?
 */
const bool hasEnvVar( std::string const & key );

/**
 * @brief obtain a path from an environment variable
 * @param key name
 * @return path
 */
const boost::filesystem::path getPathByEnv(const std::string& key);

} // namespace Common

#endif // EnvUtils_H_