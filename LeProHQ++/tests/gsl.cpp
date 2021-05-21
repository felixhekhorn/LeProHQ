#include "config.h"

#include <stdlib.h>
#include <boost/format.hpp>

int testGSL();

int main(int argc, char **argv) {
    try {
    //return testInterpolation();
    return testGSL();
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

#include "gslpp/gslpp.Matrix.hpp"
#include "gslpp/gslpp.Interpolation2D.hpp"
int testGSL() {
    uint nx = 4;
    uint ny = 5;
    gslpp::Matrix m (nx,ny);
    m.readFromFile("test.dat");
    for (uint i = 0; i < nx-1; i++) {
        for (uint j = 0; j < ny-1; j++) {
            printf ("%d\t%g\t%d\t%g\t%g\n", i, m.get(1+i,0), j, m.get(0,1+j), m.get(1+i,1+j));
        }
    }
    printf("----\n");
    gslpp::Interpolation2D spl(gsl_interp2d_bilinear,&m);
    double x = 1.1;
    double y = 3.1;
    printf("test: %g %g %g\n",x,y,spl.eval(x,y));
    /*gslpp::Matrix m (2,n);
    m.readFromFile("test.dat");
    */
    /*gsl_vector_view r = m.row(0);
    for (uint i = 0; i < n; i++) {
        printf ("%d\t%g\n", i, r.vector.data[i]);
    }*/
    /*gslpp::Interpolation1D spl (gsl_interp_polynomial,"test.dat",n);
    for (uint i = 0; i <= 10; i++) {
        cdbl x = .2 * i;
        printf ("%d\t%g\t%g\n",i, x, spl.eval(x));
    }*/
    return EXIT_SUCCESS;
}
