import numpy as np
from scipy.integrate import quad
from LeProHQ import dq1, Adler

# taken from yadism
def dq(z, proj, cc, xi):
    if xi * (1. - z) / z <= 4 :
        return 0.0
    eta = xi / 4.0 * (1.0 / z - 1.0) - 1.0
    eta = min(eta, 1e8)
    r = (
        xi / np.pi
        / z
        * (4.0 * np.pi) ** 2
        * dq1(proj,cc, xi, eta)
    )
    return r

def test_adler():
    for proj in ("F2", "FL", "x2g1"):
        for cc in ("VV", "AA"):
            for xi in (1e-2, 1, 1e4, 1e4):
                # check quad integration against MMa integration, which should be correct up to 25 digits (before truncation)
                from_quad = quad(dq,0.,1.,args=(proj, cc, xi))
                from_mma = Adler(proj, cc, xi)
                np.testing.assert_allclose(from_quad[0],from_mma,rtol=1e-4,err_msg=f"{proj=},{cc=},{xi=}")
    
