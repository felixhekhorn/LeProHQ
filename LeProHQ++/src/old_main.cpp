#include "config.h"

#include <stdlib.h>
#include <boost/format.hpp>


#include "InclusiveLeptoProduction.h"
#include "FullyDiffLeptoProduction.h"
#include "FullyInclusiveLeptoProduction.h"
#include "FullyInclusiveSoftResummedLeptoProduction.h"

int testPartonic();
int testHadronic();
int testHadronicDiff();
int testHadronicInclusiveDiff();
int testHadronicRes();
int testLeptonic();

int main(int argc, char **argv) {
    try {
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

    //FullyInclusiveLeptoProduction o(nlf,m2);
    //printf("[INFO] FullyInclusiveLeptoProduction(nlf=%d,m2=%g)\n",nlf,m2);
    
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
        cout << "l=" << 20.-10.*log(eta) << endl;
        o.setPath(2,M_PI_4,20.-10.*log(eta));
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

int testHadronic() {
    cuint nlf = 3;
    cdbl m2 = pow(1.5,2);
    const str unpolPDF = "gonly"; // "MSTW2008nlo90cl";
    const str polPDF = "gonly"; // "DSSV2014";
    /*cuint nlf = 4;
    cdbl m2 = pow(4.75,2);
    const str unpolPDF = "gonly";
    const str polPDF = "gonly";*/
    cdbl Q2 = 1e1;
    const DynamicScaleFactors mu02 (4.,1., 0., 0.);
    
    /*cdbl Delta = 1.e-6;
    InclusiveLeptoProduction o2(nlf,m2,Delta);
    InclusiveLeptoProduction oL(nlf,m2,Delta);
    InclusiveLeptoProduction oP(nlf,m2,Delta);*/

    /*cdbl xTilde = .8;
    cdbl omega = 1.;
    cdbl deltax = 1e-6;
    cdbl deltay = 7e-6;
    FullyDiffLeptoProduction o2(nlf,m2,xTilde,omega,deltax,deltay);
    FullyDiffLeptoProduction oL(nlf,m2,xTilde,omega,deltax,deltay);
    FullyDiffLeptoProduction oP(nlf,m2,xTilde,omega,deltax,deltay);*/
    
    /*FullyInclusiveLeptoProduction o2(nlf,m2);
    FullyInclusiveLeptoProduction oL(nlf,m2);
    FullyInclusiveLeptoProduction oP(nlf,m2);*/
    
    FullyInclusiveSoftResummedLeptoProduction o2(nlf,m2);
    FullyInclusiveSoftResummedLeptoProduction oL(nlf,m2);
    FullyInclusiveSoftResummedLeptoProduction oP(nlf,m2);
    list<dbl> lxgrid = FullyInclusiveSoftResummed::Interpolation::make_log_lin_grid(30,30);
    vector<dbl> xgrid {make_move_iterator(begin(lxgrid)),make_move_iterator(end(lxgrid))};
    cuint polynomial_degree = 3;
    o2.setInterpolation(xgrid, polynomial_degree);oL.setInterpolation(xgrid, polynomial_degree);oP.setInterpolation(xgrid, polynomial_degree);
    

    o2.setProjection(F2);oL.setProjection(FL);oP.setProjection(x2g1);
    
    o2.setQ2(Q2);oL.setQ2(Q2);oP.setQ2(Q2);
    o2.setMu2(mu02);oL.setMu2(mu02);oP.setMu2(mu02);
    o2.setPdf(unpolPDF,0);oL.setPdf(unpolPDF,0);oP.setPdf(polPDF,0);
    o2.setAlphaSByLHAPDF(unpolPDF,0);oL.setAlphaSByLHAPDF(unpolPDF,0);oP.setAlphaSByLHAPDF(unpolPDF,0);
    //cdbl lambdaQCD = .194;
    //o2.setLambdaQCD(lambdaQCD);oL.setLambdaQCD(lambdaQCD);oP.setLambdaQCD(lambdaQCD);
    
    //o2.flags().useLeadingOrder = oL.flags().useLeadingOrder = oP.flags().useLeadingOrder = false;
    o2.flags().useNextToLeadingOrder = oL.flags().useNextToLeadingOrder = oP.flags().useNextToLeadingOrder = false;
    //o2.flags().useGluonicChannel = oL.flags().useGluonicChannel = oP.flags().useGluonicChannel = false;
    //o2.flags().useQuarkChannel = oL.flags().useQuarkChannel = oP.flags().useQuarkChannel = false;
    o2.flags().usePhotonZ = oL.flags().usePhotonZ = oP.flags().usePhotonZ = false;
    o2.flags().useZ = oL.flags().useZ = oP.flags().useZ = false;
    
    //o2.getIntegrationConfig("F")->verbosity = 1;
    //oL.getIntegrationConfig("F")->verbosity = 1;
    //oP.getIntegrationConfig("F")->verbosity = 1;
    
    cuint N = 31;
    for (uint j = 0; j < N; ++j) {
        cdbl x = pow(10,0.-1./(N-1)*(cdbl)j);
        o2.setXBjorken(x);oL.setXBjorken(x);oP.setXBjorken(x);
        cdbl c2 = o2.F();
        cdbl e2 = o2.getIntegrationOutput().error;
        cdbl cL = oL.F();
        cdbl eL = oL.getIntegrationOutput().error;
        cdbl cP = oP.F();
        cdbl eP = oP.getIntegrationOutput().error;
        printf("%e\t% e\t% e\t% e\t% e\t% e\t% e\n",x,c2,e2,cL,eL,cP,eP);
        
        /*oL.flags() = Flags(true,false,true,false,true,false,false);
        cdbl c0g = oL.F();
        oL.flags() = Flags(true,false,false,true,true,false,false);
        cdbl c1g = oL.F();
        oL.flags() = Flags(false,true,false,true,true,false,false);
        cdbl c1q = oL.F();
        oL.flags() = Flags(true,true,true,true,true,false,false);
        cdbl c = oL.F();
        printf("%e\t % e+\t% e+\t% e=\t% e==\t% e\n",x,c0g,c1g,c1q,c0g+c1g+c1q,c);*/
    }
    
    return EXIT_SUCCESS;
}

int testHadronicDiff() {
    cuint nlf = 3;
    cdbl m2 = pow(4.5,2);
    cdbl Q2 = 1e4;
    const DynamicScaleFactors mu02 (4.,1., 0., 0.);
    
#define useHI2 1
    
#ifdef useHI2
    cdbl Delta = 1.e-6;
    InclusiveLeptoProduction o2(nlf,m2,Delta);
    InclusiveLeptoProduction oL(nlf,m2,Delta);
    InclusiveLeptoProduction oP(nlf,m2,Delta);
#else // ifdef useHI2
    cdbl xTilde = .8;
    cdbl omega = 1.;
    cdbl deltax = 1e-6;
    cdbl deltay = 7e-6;
    FullyDiffLeptoProduction o2(nlf,m2,xTilde,omega,deltax,deltay);
    FullyDiffLeptoProduction oL(nlf,m2,xTilde,omega,deltax,deltay);
    FullyDiffLeptoProduction oP(nlf,m2,xTilde,omega,deltax,deltay);
#endif // ifdef useHI2

    o2.setProjection(xF3);
    oL.setProjection(g4);
    oP.setProjection(gL);
    
    o2.setQ2(Q2);oL.setQ2(Q2);oP.setQ2(Q2);
    o2.setPdf("MSTW2008nlo90cl",0);oL.setPdf("DSSV2014",0);oP.setPdf("DSSV2014",0);
    o2.setMu2(mu02);oL.setMu2(mu02);oP.setMu2(mu02);
    o2.setAlphaSByLHAPDF("MSTW2008nlo90cl",0);oL.setAlphaSByLHAPDF("MSTW2008nlo90cl",0);oP.setAlphaSByLHAPDF("MSTW2008nlo90cl",0);
    //cdbl lambdaQCD = .194;
    //o2.setLambdaQCD(lambdaQCD);oL.setLambdaQCD(lambdaQCD);oP.setLambdaQCD(lambdaQCD);
    
    //o2.flags().useNextToLeadingOrder = oL.flags().useNextToLeadingOrder = oP.flags().useNextToLeadingOrder = false;
    o2.flags().useGluonicChannel = oL.flags().useGluonicChannel = oP.flags().useGluonicChannel = false;
    //o2.flags().usePhotonZ = oL.flags().usePhotonZ = oP.flags().usePhotonZ = false;
    //o2.flags().useZ = oL.flags().useZ = oP.flags().useZ = false;
    
    cdbl xBj = 1e-3;
    o2.setXBjorken(xBj);oL.setXBjorken(xBj);oP.setXBjorken(xBj);
    cout << "[INFO] xBj=" << xBj << endl;

/*
    cdbl ptMaxSc = .1;
    cuint N = 10;
#ifdef useHI2
    for (uint j = 0; j < N; ++j) {
        cdbl ptSc = ptMaxSc/(N)*(j+.5);
        cdbl c2 = o2.dF_dHAQTransverseMomentumScaling(ptSc);
        cdbl e2 = o2.getIntegrationOutput().error;
        cdbl cL = oL.dF_dHAQTransverseMomentumScaling(ptSc);
        cdbl eL = oL.getIntegrationOutput().error;
        cdbl cP = oP.dF_dHAQTransverseMomentumScaling(ptSc);
        cdbl eP = oP.getIntegrationOutput().error;
        printf("%e\t% e\t% e\t% e\t% e\t% e\t% e\n",ptSc,c2,cL,cP,e2,eL,eP);
    }
#else // ifdef useHI2
    o2.getIntegrationConfig("F")->verbosity = 1;
    oL.getIntegrationConfig("F")->verbosity = 1;
    oP.getIntegrationConfig("F")->verbosity = 1;
    / *o2.getIntegrationConfig("F")->calls = 2000;
    oL.getIntegrationConfig("F")->calls = 2000;
    oP.getIntegrationConfig("F")->calls = 7000;* /
    o2.activateHistogram(FullyDiff::histT::HAQTransverseMomentumScaling, N, "/home/Felix/Physik/PhD/data2/debug/xt-e2.dat",0.,ptMaxSc);
    oL.activateHistogram(FullyDiff::histT::HAQTransverseMomentumScaling, N, "/home/Felix/Physik/PhD/data2/debug/xt-eL.dat",0.,ptMaxSc);
    oP.activateHistogram(FullyDiff::histT::HAQTransverseMomentumScaling, N, "/home/Felix/Physik/PhD/data2/debug/xt-eP.dat",0.,ptMaxSc);
    o2.F();oL.F();oP.F();
#endif // ifdef useHI2
*/

    cdbl yMin = -4.2,yMax = 1.;
    cuint N = 20;
#ifdef useHI2
    /*{
        cdbl c2 = o2.F();
        cdbl cL = oL.F();
        cdbl cP = oP.F();
        printf("% e\t% e\t% e\n",c2,cL,cP);
    }*/
    for (uint j = 0; j < N; ++j) {
        cdbl y = yMin + (yMax-yMin)/N*(j+.5);
        cdbl c2 = o2.dF_dHAQRapidity(y);
        cdbl e2 = o2.getIntegrationOutput().error;
        cdbl cL = oL.dF_dHAQRapidity(y);
        cdbl eL = oL.getIntegrationOutput().error;
        cdbl cP = oP.dF_dHAQRapidity(y);
        cdbl eP = oP.getIntegrationOutput().error;
        printf("% e\t% e\t% e\t% e\t% e\t% e\t% e\n",y,c2,e2,cL,eL,cP,eP);
    }
#else // ifdef useHI2
    o2.getIntegrationConfig("F")->verbosity = 1;
    oL.getIntegrationConfig("F")->verbosity = 1;
    oP.getIntegrationConfig("F")->verbosity = 1;
    o2.getIntegrationConfig("F")->calls = 1000000;
    oL.getIntegrationConfig("F")->calls = 1000000;
    oP.getIntegrationConfig("F")->calls = 2000000;
    o2.getIntegrationConfig("F")->MC_adaptChi2 = false;
    oL.getIntegrationConfig("F")->MC_adaptChi2 = false;
    oP.getIntegrationConfig("F")->MC_adaptChi2 = false;
    o2.activateHistogram(FullyDiff::histT::HAQRapidity, N, "/home/Felix/Physik/PhD/data2/debug/y-e2.dat",yMin,yMax);
    oL.activateHistogram(FullyDiff::histT::HAQRapidity, N, "/home/Felix/Physik/PhD/data2/debug/y-eL.dat",yMin,yMax);
    oP.activateHistogram(FullyDiff::histT::HAQRapidity, N, "/home/Felix/Physik/PhD/data2/debug/y-eP.dat",yMin,yMax);
    o2.F();oL.F();
    oP.F();
#endif // ifdef useHI2
    
    return EXIT_SUCCESS;
}

int testHadronicInclusiveDiff() {
    cuint nlf = 3;
    cdbl m2 = pow(1.5,2);
    cdbl Q2 = 1e1;
    const DynamicScaleFactors mu02 (4.,1.,0., 0.);
    
    cdbl Delta = 1.e-6;
    InclusiveLeptoProduction o2(nlf,m2,Delta);
    o2.setProjection(F2);
    
    o2.setQ2(Q2);
    o2.setPdf("MSTW2008nlo90cl",0);
    o2.setMu2(mu02);
    o2.setAlphaSByLHAPDF("MSTW2008nlo90cl",0);
    
    //o2.flags().useLeadingOrder = false;
    o2.flags().useNextToLeadingOrder = false;
    //o2.flags().useGluonicChannel = false;
    //o2.flags().useQuarkChannel = false;
    o2.flags().usePhotonZ = false;
    o2.flags().useZ = false;
    
    cdbl xBj = 1e-3;
    o2.setXBjorken(xBj);
    cout << "[INFO] xBj=" << xBj << ", Q2="<<Q2 << endl;

    //o2.getIntegrationConfig("dF_dHAQ")->method = "gsl_integration_cquad";
    //o2.getIntegrationConfig("dF_dHAQ")->GslQag_epsabs = 1e-22;
    //o2.getIntegrationConfig("dF_dHAQ")->GslQag_epsrel = 1e-6;
    
    /*cuint N = 10;
    for (uint j = 0; j < N; ++j) {
        const DynamicScaleFactors mu02a (4.,1.,0.,0.);
        const DynamicScaleFactors mu02b (4.,1.,4.,0.);
    
        cdbl x = pow(10.,-3.+3./cdbl(N)*j);
        o2.setXBjorken(x);
        o2.setMu2(mu02a);
        cdbl c2a = o2.F();
        cdbl e2a = o2.getIntegrationOutput().error;
        o2.setMu2(mu02b);
        cdbl c2b = o2.F();
        cdbl e2b = o2.getIntegrationOutput().error;
        printf("% e\t% e\t% e\t% e\t% e\t% e\n",x,c2a/c2b-1.,c2a,e2a,c2b,e2b);
    }*/
    
    /*cdbl ptMaxSc = 1;
    cuint N = 5;
    for (uint j = 0; j < N; ++j) {
        cdbl ptSc = ptMaxSc/(N)*(j+.5);
        cdbl c2 = o2.dF_dHAQTransverseMomentumScaling(ptSc);
        cdbl e2 = o2.getIntegrationOutput().error;
        printf("%e\t% e\t% e\n",ptSc,c2,e2);
    }*/
    
    /*o2.setPdf("DSSV2014",5);
    o2.setProjection(x2g1);
    o2.setMu2(mu02);
    cuint N = 5;
    cdbl pt = .5;
    cdbl c0 = o2.dF_dHAQTransverseMomentum(pt);
    printf("c0 = % e\n",c0);
    for (uint j = 0; j < N; ++j) {
        cdbl a = pow(10.,-1.+2./cdbl(N-1)*j);
        const DynamicScaleFactors mu2a = mu02*a;
        o2.setMu2(mu2a);
        cdbl c2 = o2.dF_dHAQTransverseMomentum(pt);
        cdbl e2 = o2.getIntegrationOutput().error;
        cout << boost::format("%e\t% e\t% e\t% e\t")%a%c2%e2%(c2/c0) << endl;
    }*/
    
    cdbl y = -4.075700;
    cdbl c2 = o2.dF_dHAQRapidity(y);
    cdbl e2 = o2.getIntegrationOutput().error;
    printf("%e\t% e\t% e\n",y,c2,e2);
        
    /*cdbl y0 = 3.;
    cuint N = 10;
    for (uint j = 0; j < N; ++j) {
        cdbl y = y0*(-1.+2./N*(j+.5));
        cdbl c2 = o2.dF_dHAQRapidity(y);
        cdbl e2 = o2.getIntegrationOutput().error;
        printf("%e\t% e\t% e\n",y,c2,e2);
    }*/
    
    /*cdbl ymin = -4.2, ymax=-4.1;
    cuint N = 10;
    for (uint j = 0; j < N; j+=2) {
        cdbl y = ymin + (ymax-ymin)/cdbl(N-1)*j;
        o2.setPdf("MSTW2008nlo90cl",0);
        o2.setProjection(F2);
        cdbl c2 = o2.dF_dHAQRapidity(y);
        cdbl e2 = o2.getIntegrationOutput().error;
        //printf("proj = FL\n");
        //o2.setProjection(FL);
        cdbl cL = 0.;//o2.dF_dHAQRapidity(y);
        cdbl eL = 0.;//o2.getIntegrationOutput().error;
        //o2.setPdf("DSSV2014",0);
        //o2.setProjection(x2g1);
        cdbl cg = 0.;//o2.dF_dHAQRapidity(y);
        cdbl eg = 0.;//o2.getIntegrationOutput().error;
        printf("%e\t% e\t% e\t% e\t% e\n",y,c2-cL,sqrt(e2*e2+eL*eL),cg,eg);
    }*/
    
    return EXIT_SUCCESS;
}

int testHadronicRes() {
    cuint nlf = 4;
    cdbl m2 = pow(4.5,2);
    cdbl Q2 = 1e1;
    const DynamicScaleFactors mu02 (4.,1.,0., 0.);
    
    FullyInclusiveSoftResummedLeptoProduction o(nlf,m2);
    cout << "[INFO] FullyInclusiveSoftResummedLeptoProduction(" << nlf << ", "<<m2 << ")" << endl;
    o.setProjection(F2);
    
    o.setQ2(Q2);
    o.setPdf("MSTW2008nlo90cl",0);
    o.setMu2(mu02);
    o.setAlphaSByLHAPDF("MSTW2008nlo90cl",0);
    
    list<dbl> lxgrid = FullyInclusiveSoftResummed::Interpolation::make_log_lin_grid(7,10);
    vector<dbl> xgrid {make_move_iterator(begin(lxgrid)),make_move_iterator(end(lxgrid))};
    o.setInterpolation(xgrid, 3);
    dbl zmax = Q2 / (4.*m2 + Q2);
    cout << "[INFO] zmax = " << zmax << endl;
    
    for (auto xBj : {6.309573e-02,3.981072e-02,2.511886e-02}){
        o.setXBjorken(xBj);
        cout << "[INFO] xBj=" << xBj << ", Q2="<<Q2 << endl;
        cout << xBj << "\t" << o.F() << endl;
    }
    return EXIT_SUCCESS;
}

int testLeptonic() {
    cuint nlf = 3;
    cdbl m2 = pow(1.5,2);
    cdbl Delta = 1.e-6;
    DynamicScaleFactors mu02 (4.,1., 0., 0.);
    cdbl lambdaQCD = .239;
    
    InclusiveLeptoProduction oF(nlf,m2,Delta);
    oF.setPolarizeBeams(false);
    InclusiveLeptoProduction og(nlf,m2,Delta);
    og.setPolarizeBeams(true);
    
    oF.setPdf("MSTW2008nlo90cl",0);og.setPdf("DSSV2014",0);
    oF.setMu2(mu02);og.setMu2(mu02);
    oF.setLambdaQCD(lambdaQCD);og.setLambdaQCD(lambdaQCD);
    
    oF.flags().useNextToLeadingOrder = og.flags().useNextToLeadingOrder = false;
    /*oF.flags().useGluonicChannel = og.flags().useGluonicChannel = false;*/
    oF.flags().usePhotonZ = og.flags().usePhotonZ = false;
    oF.flags().useZ = og.flags().useZ = false;
    
    cdbl Q2min = 2.;
    oF.setQ2min(Q2min);og.setQ2min(Q2min);
    
    uint N = 11;
    for (uint j = 0; j < N; ++j) {
        cdbl Sl = 200. + j*50.;
        oF.setLeptonicS(Sl);og.setLeptonicS(Sl);
        cdbl sig = oF.sigma();
        cdbl DeltaSig = og.sigma();
        printf("%e\t%e\t%e\n",Sl,sig,DeltaSig);
    }
    
    return EXIT_SUCCESS;
}