#include "config.h"

#include <stdlib.h>
#include <boost/format.hpp>

#include "FullyInclusiveSoftResummed/Interpolation.h"

int testInterpolation();

int main(int argc, char **argv) {
    try {
    return testInterpolation();
    //return testGSL();
    //return testME();
    //return testPartonic();
    //return testHadronic();
    //return testHadronicDiff();
    //return testHadronicInclusiveDiff();
    //return testHadronicRes();
    //return testLeptonic();
    } catch(const std::exception& e) {
        cout << "Hoppala ..." << endl << e.what();
    }
    
    return EXIT_SUCCESS;
}

int testInterpolation() {
    list<dbl> lxgrid = FullyInclusiveSoftResummed::Interpolation::make_log_lin_grid(7,10);
    vector<dbl> xgrid {make_move_iterator(begin(lxgrid)),make_move_iterator(end(lxgrid))};
    for (auto& it: xgrid)
        cout << it << endl;
    cout << "---" << endl;
    //const vector<dbl> xgrid = {exp(-1), 1.0};
    FullyInclusiveSoftResummed::Interpolation::Interpolator itp(xgrid, 1);
    cout << "0,.4 " << itp.evaluatePolynomX(0, .4) << endl;
    cout << "1,.4 " << itp.evaluatePolynomX(1, .4) << endl;
    /*Common::Interpolation::dcmplx N = 1.;
    dbl logx = -1.5;
    cout << itp.evaluatePolynomN(1,N,logx) << endl;
    logx = -2.;
    cout << itp.evaluatePolynomN(1,N,logx) << endl;
    N = Common::Interpolation::dcmplx(1.,1.);
    cout << itp.evaluatePolynomN(1,N,logx) << endl;*/
    /*
    logx = -2.;
    cout << itp.evaluatePolyN(0,N,logx) << "\t" << exp(-N*logx)*(1./N - 1./(N+1.)) << endl;
    N = Common::Interpolation::dcmplx(1.,1.);
    cout << itp.evaluatePolyN(0,N,logx) << "\t" << exp(-N*logx)*(1./N - 1./(N+1.)) << endl;*/
    /*PdfWrapper pdf("MSTW2008nlo90cl",0);
    cdbl q2 = 10;
    const int pid = 21;
    cdbl x = 3e-3;
    itp.loadPdfData(pdf, pid, q2);
    cout << "g() = " << itp.evaluateN(complex<double>(1.,1.),-2) << endl;*/
    return EXIT_SUCCESS;
}