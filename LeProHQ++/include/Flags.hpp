#ifndef Flags_HPP_
#define Flags_HPP_

#include <yaml-cpp/yaml.h>

/**
 * @class Flags
 * @brief controls active channels, orders and bosons
 */
struct Flags {
    
    /** @name partonic channels: gluon, light quark */
    ///@{
    /** @brief use gluonic initial state? */
    bool useGluonicChannel = true;
    
    /** @brief use light quark initial state? */
    bool useQuarkChannel = true;
    ///@}
    
    /** @name computing orders = number of loops */
    ///@{
    /** @brief use leading order calculations? */
    bool useLeadingOrder = true;
    
    /** @brief use (pure) next-to-leading order calculations? */
    bool useNextToLeadingOrder = true;
    ///@}
    
    /** @name exchanged bosons: photon, Z */
    ///@{
    /** @brief use photon exchange? */
    bool usePhoton = true;
    
    /** @brief use interference between photon and Z exchange? */
    bool usePhotonZ = true;
    
    /** @brief use Z exchange? */
    bool useZ = true;
    ///@}
    
    /**
     * @brief constructor
     * @param useGluonicChannel use gluonic initial state?
     * @param useQuarkChannel use light quark initial state?
     * @param useLeadingOrder use leading order calculations?
     * @param useNextToLeadingOrder use (pure) next-to-leading order calculations?
     * @param usePhoton use photon exchange?
     * @param usePhotonZ use interference between photon and Z exchange?
     * @param useZ use Z exchange?
     */
    Flags(bool useGluonicChannel, bool useQuarkChannel,
            bool useLeadingOrder, bool useNextToLeadingOrder,
            bool usePhoton, bool usePhotonZ, bool useZ) :
        useGluonicChannel(useGluonicChannel), useQuarkChannel(useQuarkChannel), 
        useLeadingOrder(useLeadingOrder), useNextToLeadingOrder(useNextToLeadingOrder),
        usePhoton(usePhoton), usePhotonZ(usePhotonZ), useZ(useZ) {}

/**
 * @brief return YAML representation as string
 * @return string
 */
    str toYAML() const {
        YAML::Node r;
        r["useGluonicChannel"] = this->useGluonicChannel;
        r["useQuarkChannel"] = this->useQuarkChannel;
        r["useLeadingOrder"] = this->useLeadingOrder;
        r["useNextToLeadingOrder"] = this->useNextToLeadingOrder;
        r["usePhoton"] = this->usePhoton;
        r["usePhotonZ"] = this->usePhotonZ;
        r["useZ"] = this->useZ;
        r.SetStyle(YAML::EmitterStyle::Flow);
        YAML::Emitter e;
        e << r;
        return e.c_str(); 
    }

/**
 * @brief writes object as YAML to stream
 * @param os stream
 * @param c object
 * @return improved stream
 */
    friend std::ostream& operator<<(std::ostream& os, const Flags& c) {
        os << c.toYAML();
        return os;
    }
};

#endif // Flags_HPP_