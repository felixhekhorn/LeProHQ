import timeit

import numpy as np
from scipy import integrate

import LeProHQ


def fake_int(xi):
    return integrate.quad(lambda z: LeProHQ.cg1("F2", "VV", xi, z), 0, 10)


[fake_int(xi) for xi in np.linspace(1, 10, 25)]
