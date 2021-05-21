#ifndef PDFWRAPPER_H
#define PDFWRAPPER_H

#include <string>

#include <boost/filesystem.hpp>

#include <LHAPDF/LHAPDF.h>

/**
 * @brief Wrapper to take different pdf sources into account
 * 
 * Defaults to LHAPDF\cite LHAPDF6 and delivers ports to
 * - MorfinTungB\cite Morfin:1990ck implemented as *MorfinTungB*
 * - CTEQ3\cite Lai:1994bb implemented as *CTEQ3M*
 * - GRV94\cite Gluck:1994uf implemented as *GRV94LO* and *GRV94NLO*
 * - GRSV96\cite Gluck:1995yr for polarized PDFs implemented as *GRSV96STDLO* and *GRSV96STDNLO* [NB: patched version (name collision)]
 * - DSSV2014\cite PhysRevLett.113.012001 for polarized PDFs implemented as *DSSV2014* with support for all members [NB: patched version (member loading)]
 * 
 * The location of the member file for DSSV2014 is read from the environment variable *DSSV2014_GRIDS*.
 * 
 * The location of the files for GRSV96 is read from the environment variable *GRSV96_GRIDS*
 * 
 * For PDG particle id see http://pdg.lbl.gov/mc_particle_id_contents.html
 */
class PdfWrapper {

/**
 * @brief port to DSSV?
 */
    bool isDSSV2014;

/**
 * @brief port to CTEQ3?
 */
    bool isCTEQ3;

/**
 * @brief port to GRSV96?
 */
    bool isGRSV96;

/**
 * @brief port to GRV94?
 */
    bool isGRV94;

/**
 * @brief port to MorfinTungB?
 */
    bool isMorfinTungB;
    
/**
 * @brief set name
 */
    std::string setname;

/**
 * @brief member id
 */
    int member;
    
/**
 * @brief parton distribution function
 */
    LHAPDF::PDF* lha;
    
/**
 * @brief path to active GRSV96 file
 */
    std::string GRSV96_path;
    
/**
 * @brief reads a path from an environment variable
 * @param pathName environment variable
 * @return path
 */
    const boost::filesystem::path getPathByEnv(const std::string& pathName);
    
public:

/**
 * @brief constructor
 * @param setname set name
 * @param member member id
 */
    PdfWrapper(const std::string &setname, const int member);
    
/**
 * @brief destructor
 */
    ~PdfWrapper();
    
/**
 * @brief return setname
 * @return setname
 */
    const std::string getSetname() const { return this->setname; };
    
/**
 * @brief return member id
 * @return member id
 */
    const int getMember() const { return this->member; };

/**
 * @brief valid range for PDF?
 * @param Q2 energy scale to check
 * @return in range?
 */
    const bool inRangeQ2(const double Q2) const;

/**
 * @brief return x*distribution for particle pid
 * @param pid PDG particle number
 * @param x
 * @param Q2
 * @return x*f(x,Q²)
 */
    const double xfxQ2(const int pid, const double x, const double Q2) const;

};

#endif // PDFWRAPPER_H
