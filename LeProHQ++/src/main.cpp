#include "config.h"

#include <stdlib.h>
#include <boost/format.hpp>

#include "InclusiveLeptoProduction.h"
#include "FullyDiffLeptoProduction.h"
#include "FullyInclusiveLeptoProduction.h"
#include "FullyInclusiveSoftResummedLeptoProduction.h"

int testPartonic();

int main(int argc, char **argv) {
    try {
    return testPartonic();
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



int testPartonic() {
    cuint nlf = 4;
    cdbl m2 = pow(4.75,2);
    cdbl Q2 = 1e2;

    /*cdbl Delta = 1.e-6;
    InclusiveLeptoProduction o(nlf,m2,Delta);
    printf("[INFO] InclusiveLeptoProduction(%d,%g,%g)\n",nlf,m2,Delta);*/
    
    /*cdbl xTilde = .8;
    cdbl omega = 1.;
    cdbl deltax = 1e-6;
    cdbl deltay = 7e-6;
    FullyDiffLeptoProduction o(nlf,m2,xTilde,omega,deltax,deltay);
    printf("[INFO] FullyDiffLeptoProduction(%d,%g,%g,%g,%g,%g)\n",nlf,m2,xTilde,omega,deltax,deltay);*/

    /*FullyInclusiveLeptoProduction o(nlf,m2);
    printf("[INFO] FullyInclusiveLeptoProduction(nlf=%d,m2=%g)\n",nlf,m2);*/
    
    FullyInclusiveSoftResummedLeptoProduction o(nlf,m2);
    printf("[INFO] FullyInclusiveSoftResummedLeptoProduction(nlf=%d,m2=%g)\n",nlf,m2);

    o.setQ2(Q2);
    printf("[INFO] Q2 = %g\n",Q2);
    cuint N = 11;
    for (uint j = 0; j < N; ++j) {
        cdbl eta = pow(10.,-3.+3./(N-1)*j);
        o.setPartonicEta(eta);
        //o.getIntegrationConfig("cq1_AA")->verbosity = 1;
        //o.getIntegrationConfig("cg1_VV")->method = "gsl_monte_vegas_integrate";
        //o.getIntegrationConfig("cg1_VV")->calls = 40000;
        //o.getIntegrationConfig("cg1_VV")->Dvegas_bins = 40;
        //o.setDelta(eta/1000.);
        //cout << "l=" << 20.-10.*log(eta) << endl;
        //o.setPath(2,M_PI_4,20.-10.*log(eta));
        o.setProjection(F2);
        cdbl cF2 = o.cg0t_VV();
        o.setProjection(FL);
        cdbl cFL = o.cg0t_VV();
        o.setProjection(x2g1);
        cdbl cx2g1 = o.cg0t_VV();
        cout << boost::format("%e\t% e\t% e\t% e")%eta%(cF2)%cFL%cx2g1 << endl;
    }
    
    return EXIT_SUCCESS;
}