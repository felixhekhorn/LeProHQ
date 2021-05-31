#include "Pdf/PdfWrapper.h"

#include "config-buildtime.hpp"

#include "Common/EnvUtils.h"

/**
 * @brief converst a c-style string to a fortran-style string
 * @see http://stackoverflow.com/questions/10163485/passing-char-arrays-from-c-to-fortran
 * @param fstring pointer to fortran-style string
 * @param fstring_len length of fortran-style string
 * @param cstring c-style string
 */
void ConvertToFortran(char* fstring, const std::size_t fstring_len, const char* cstring);
void ConvertToFortran(char* fstring, const std::size_t fstring_len, const char* cstring) {
    std::size_t inlen = std::strlen(cstring);
    std::size_t cpylen = std::min(inlen, fstring_len);
    // raise truncation error or warning?
    std::copy(cstring, cstring + cpylen, fstring);
    std::fill(fstring + cpylen, fstring + fstring_len, ' ');
}

/** @brief Fortran wrappers */
const unsigned int DSSV14_fp_len = 100;
const unsigned int GRSV96_fp_len = 130;
extern "C" {
// MorfinTungB
#ifdef HAVE_MORFINTUNGB
double tungb_(double* X, double* Q2, double* UV, double* DV, double* GL, double* UBAR, double* CSEA, double* BSEA, double* TSEA, int* IFLAG);
#endif
// CTEQ3
#ifdef HAVE_CTEQ3
double ctq3pd_(int* Iset, int* Iparton, double* X, double* Q, int* Irt);
#endif
// GRV94
#ifdef HAVE_GRV
double grv94lo_(double* X, double* Q2, double* UV, double* DV, double* DEL, double* UDB, double* SB, double* GL);
double grv94ho_(double* X, double* Q2, double* UV, double* DV, double* DEL, double* UDB, double* SB, double* GL);
#endif
// GRSV96
#ifdef HAVE_GRSV
extern struct {
   int iini;
} intini_;
double parpol_(char path[GRSV96_fp_len], double* X, double* Q2, double* UV, double* DV, double* QB, double* ST, double* GL);
#endif
// DSSV14
#ifdef HAVE_DSSV
void dssvini_(char rpath[DSSV14_fp_len], int* member);
void dssvgupdate_(double* X, double* Q2, double* DUV, double* DDV, double* DUBAR, double* DDBAR, double* DSTR, double* DGLU);
#endif
}
PdfWrapper::PdfWrapper(const std::string &setname, const int member) : setname(setname), member(member), lha(0) {
    this->isDSSV14 = ("DSSV14" == setname);
    this->isCTEQ3 = ("CTEQ3M" == setname);
    this->isGRSV96 = ("GRSV96STDLO" == setname) || ("GRSV96STDNLO" == setname);
    this->isGRV94 = ("GRV94LO" == setname) || ("GRV94NLO" == setname);
    this->isMorfinTungB = ("MorfinTungB" == setname);
    if ((this->isCTEQ3 || this->isGRSV96 || this->isGRV94 || this->isMorfinTungB) && 0 != member)
        throw LHAPDF::UserError("pdf "+setname+" has only a central member!");
    /** @todo add verbosity flag as member? */
    if (this->isDSSV14) {
#ifdef HAVE_DSSV
        // setup path
        int m = member;
        /** @todo bind Common::getPathByEnv again to PdfWrapper and specialize again the error */
        const boost::filesystem::path path = Common::getPathByEnv("DSSV14_GRIDS");
        char rpath[DSSV14_fp_len];
        ConvertToFortran(rpath,DSSV14_fp_len,path.c_str());
        std::cout << "[INFO] PdfWrapper loading "<<setname<<" member #"<<member<<" from "<<path<<std::endl;
        // init
        dssvini_(rpath,&m);
#else
        throw std::runtime_error("DSSV14 support not enabled");
#endif
    } else if(this->isGRSV96) {
#ifdef HAVE_GRSV
        // find path
        boost::filesystem::path p = Common::getPathByEnv("GRSV96_GRIDS");
        // add file
        if ("GRSV96STDLO" == setname) p /= "STDLO.GRID";
        else if ("GRSV96STDNLO" == setname) p /= "STDNLO.GRID";
        this->GRSV96_path = p.string();
        std::cout << "[INFO] PdfWrapper loading "<<setname<<" from "<<this->GRSV96_path<<std::endl;
        intini_.iini = 0;
#else
        throw std::runtime_error("GRSV96 support not enabled");
#endif
    } // some analytic PDFs (and caching deactivated)
    else if(this->isCTEQ3) {
#ifdef HAVE_CTEQ3
        std::cout << "[INFO] PdfWrapper loading "<<setname<<std::endl;
#else
        throw std::runtime_error("CTEQ3 support not enabled");
#endif
    } else if(this->isGRV94) {
#ifdef HAVE_GRV
        std::cout << "[INFO] PdfWrapper loading "<<setname<<std::endl;
#else
        throw std::runtime_error("GRV94 support not enabled");
#endif
    } else if(this->isMorfinTungB) {
#ifdef HAVE_MORFINTUNGB
        std::cout << "[INFO] PdfWrapper loading "<<setname<<std::endl;
#else
        throw std::runtime_error("MorfinTungB support not enabled");
#endif
    } else {
        this->lha = LHAPDF::mkPDF(setname,member);
    }
}

PdfWrapper::~PdfWrapper() {
    if (0 != this->lha)
        delete (this->lha);
}

const bool PdfWrapper::inRangeQ2(const double Q2) const {
    if (this->isDSSV14) {
        return (Q2 >= .8)      && (Q2 <= 1e6);
    } else if (this->isCTEQ3) {
        return (Q2 >= 1.6*1.6) && (Q2 <= 10e3*10e3);
    } else if (this->isGRSV96) {
        return (Q2 >= .4)      && (Q2 <= 1e4);
    } else if (this->isGRV94) {
        return true; // no limits of fit given
    } else if (this->isMorfinTungB) {
        return (Q2 >= 3.*3.) && (Q2 <= 1e4*1e4);
    } 
    return this->lha->inRangeQ2(Q2);
}

const double PdfWrapper::xfxQ2(const int pid, const double x, const double Q2) const {
    if (this->isDSSV14) {
#ifdef HAVE_DSSV
        double x_ = x, Q2_ = Q2,uv, dv, ubar, dbar, s, glu;
        dssvgupdate_(&x_, &Q2_, &uv, &dv, &ubar, &dbar, &s, &glu);
        if (21 == pid) return glu;
        if ( 1 == pid) return uv + ubar;
        if (-1 == pid) return ubar;
        if ( 2 == pid) return dv + dbar;
        if (-2 == pid) return dbar;
        if (-3 == pid || 3 == pid) return s;
        return 0.;
#else
        throw std::runtime_error("DSSV support not enabled");
#endif
    } else if (this->isCTEQ3) {
#ifdef HAVE_CTEQ3
        int set, ret;
        if ("CTEQ3M" == setname) set = 1;
        else throw LHAPDF::UserError("unknown CTEQ3 set \""+setname+"\"!");
        int parton = (pid == 21) ? 0 : pid;
        double x_ = x, Q_ = sqrt(Q2);
        double f = ctq3pd_(&set,&parton,&x_,&Q_,&ret);
        if (1 == ret)
            return 0.;
        return x*f;
#else
        throw std::runtime_error("CTEQ3 support not enabled");
#endif
    } else if (this->isGRSV96) {
#ifdef HAVE_GRSV
        char path[GRSV96_fp_len];
        ConvertToFortran(path,GRSV96_fp_len,this->GRSV96_path.c_str());
        double x_ = x, Q2_ = Q2,uv, dv, qb, st, gl;
        parpol_(path,&x_,&Q2_,&uv,&dv,&qb,&st,&gl);
        if (21 == pid) return gl;
        if ( 1 == pid) return uv + qb;
        if (-1 == pid) return qb;
        if ( 2 == pid) return dv + qb;
        if (-2 == pid) return qb;
        if (-3 == pid || 3 == pid) return st;
        return 0.;
#else
        throw std::runtime_error("GRSV support not enabled");
#endif
    } else if (this->isGRV94) {
#ifdef HAVE_GRV
        double x_ = x, Q2_ = Q2,uv, dv, del, udb, sb, gl;
        if ("GRV94LO" == setname)       grv94lo_(&x_,&Q2_,&uv,&dv,&del,&udb,&sb,&gl);
        else if ("GRV94NLO" == setname) grv94ho_(&x_,&Q2_,&uv,&dv,&del,&udb,&sb,&gl);
        else throw LHAPDF::UserError("unknown GRV94 set \""+setname+"\"!");
        const double ub = (udb - del)/2.;
        const double db = (udb + del)/2.;
        if (21 == pid) return gl;
        if ( 1 == pid) return uv + ub;
        if (-1 == pid) return ub;
        if ( 2 == pid) return dv + db;
        if (-2 == pid) return db;
        if (-3 == pid || 3 == pid) return sb;
        return 0.;
#else
        throw std::runtime_error("GRV support not enabled");
#endif
    } else if (this->isMorfinTungB) {
#ifdef HAVE_MORFINTUNGB
        int flag = 2;
        double x_ = x, Q2_ = Q2, uv, dv, gl, ubar, csea, bsea, tsea;
        tungb_(&x_, &Q2_, &uv, &dv, &gl, &ubar, &csea, &bsea, &tsea, &flag);
        if (21 == pid) return gl;
        if ( 1 == pid) return uv + ubar;
        if (-1 == pid) return ubar;
        if ( 2 == pid) return dv + ubar;
        if (-2 == pid) return ubar;
        if (-3 == pid || 3 == pid) return ubar;
        if (-4 == pid || 4 == pid) return csea;
        if (-5 == pid || 5 == pid) return bsea;
        if (-6 == pid || 6 == pid) return tsea;
        return 0.;
#else
        throw std::runtime_error("MORFINTUNGB support not enabled");
#endif
    }
    return this->lha->xfxQ2(pid,x,Q2);
}