#include "config.h"

#include <stdlib.h>
#include <boost/format.hpp>

#include "InclusiveLeptoProduction.h"
#include "FullyDiffLeptoProduction.h"
#include "FullyInclusiveLeptoProduction.h"
#include "FullyInclusiveSoftResummedLeptoProduction.h"

//#include "pineappl_capi.h"
#include "PineAPPLpp.hpp"

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

int testPineAPPL() {
    cout << "Bla" << endl;

    PineAPPL::KeyVal k;
    k.set("d", 21.123);
    k.set("i", 2);
    k.set("b", false);
    k.set("s", std::string("äöüß"));
    cout << k.get_double("d")<< "\t" << k.get_int("i") << "\t" << k.get_bool("b")<< "\t" << k.get_string("s");

    /* // create a new luminosity function for the $\gamma\gamma$ initial state
    auto* lumi = pineappl_lumi_new();
    int32_t pdg_ids[] = { 11, 21 };
    double weights[] = { 1.0 };
    pineappl_lumi_add(lumi, 1, pdg_ids, weights);

    // only LO $\alpha_\mathrm{s}^0 \alpha^2 \log^0(\xi_\mathrm{R}) \log^0(\xi_\mathrm{F})$
    uint32_t orders[] = { 0, 1, 0, 0 };

    // we bin in rapidity from 0 to 2.4 in steps of 0.1
    double bins[] = {
        0.0,1.0
    };

    // create the PineAPPL grid with default interpolation and binning parameters
    auto* grid = pineappl_grid_new(lumi, 1, orders, 1, bins, keyval);

    // now we no longer need `keyval` and `lumi`
    pineappl_lumi_delete(lumi);
    // destroy the object
    pineappl_grid_delete(grid);*/
    return EXIT_SUCCESS;
}

int main(int argc, char **argv) {
    try {
        return testPineAPPL();
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
