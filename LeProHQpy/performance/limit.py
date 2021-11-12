import timeit

import numpy as np

import LeProHQ
from LeProHQ import bmsn

for xi in [1.,1e1,1e2,1e3,1e4]:
    print(f"xi = {xi}")
    a = []
    b = []
    zs = []
    for eta in [1e-1,1.,1e1]:
        # eta(z) = xi/4*(1/z-1) -1
        # (1+eta)(4/xi) = 1/z-1
        z = 1./(1. + (1. + eta)*4./xi)
        zs.append(z)

        norm = xi/np.pi * z # * (4.*np.pi)**2

        a.append(norm * LeProHQ.cg0("FL","VV", xi, eta))
        b.append(bmsn.clg1am0_a0(z))
    a = np.array(a)
    b = np.array(b)
    print("zs:", zs)
    print(a)
    print(b)
    print(a/b)
