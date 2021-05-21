#include "config.h"

#include <stdlib.h>
#include <boost/format.hpp>

int testME();

int main(int argc, char **argv) {
    try {
    return testME();
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

/*#include "../src/FullyDiff/ME/R.h"
#include "../src/FullyDiff/ME/RCounterX.h"
#include "../src/FullyDiff/ME/RCounterY.h"
#include "../src/FullyDiff/ME/RCounterXY.h"
#include "../src/FullyDiff/ME/SV.h"*/
#include "../src/FullyDiff/ME/A3.h"
int testME() {
    cdbl m2 = 1.;
    cdbl q2 = -10.;
    cdbl sp = 5.;
    cdbl t1 = -3.;
    cdbl u1 = -1.;
    cdbl tp = -1.;
    cdbl up = -1.;
    cout << FullyDiff::ME::A3_x2g1_VV(m2,q2,sp,t1,u1,tp,up) << endl;
    
    /*cdbl m2 = 22.5625;
    cdbl q2 = -1e3;
    cdbl sp = 1090.34025;
    cdbl xE = 0.999958608836369;
    cdbl yE = 3.4999999999341114e-06;
    cdbl Theta1 = 1.5707963267948966;
    cdbl Theta2 = 1.5707963267948966;
    cout << FullyDiff::ME::R_FL_VV(m2,q2,sp,xE,yE,Theta1,Theta2) << endl;*/
    
    /*cdbl m2 = 1.;
    cdbl q2 = -10.;
    cdbl sp = 5.;
    cdbl x = .3;
    cdbl y = 0.;
    cdbl Theta1 = .5;
    cdbl Theta2 = 1.;
    cdbl t1 = -sp/2.*(1.-sqrt(1. - 4.*m2/(sp+q2))*cos(Theta1));
    cout << (FullyDiff::ME::RCounterY_FL_VV(m2,q2,sp,x,Theta1,Theta2)) << endl;*/
    
    /*cdbl m2 = 1.;
    cdbl q2 = -10.;
    cdbl s = 4.*m2*(1.+1.);
    cdbl sp = s-q2;
    cdbl Theta1 = 1.;
    cdbl t1 = -sp/2.*(1.-sqrt(1. - 4.*m2/(sp+q2))*cos(Theta1));
    cdbl betaTilde = .9;
    cout << (FullyDiff::ME::SVQED_FL_VV(m2,q2,sp,t1,betaTilde)) << endl;*/
    return EXIT_SUCCESS;
}
