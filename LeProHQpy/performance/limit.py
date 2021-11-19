import numpy as np
import matplotlib.pyplot as plt

import LeProHQ
from LeProHQ import bmsn, bfkl
from LeProHQ.cg0 import cg0t
from LeProHQ.cg1 import cg1t

xis = np.geomspace(1e1, 1e2 - 1, 5)
etas = np.geomspace(1e2, 1e5, 10)
us = np.geomspace(1e-1, 1.0 - 1e-4, 5)

for xi in xis:
    zmax = 1.0 / (1.0 + 4.0 / xi)
    print(f"xi = {xi}, zmax = {zmax}")
    lxi = np.log(xi)
    exas = []
    hscs = []
    thrs = []
    ns = []
    f1s = []
    # zs = []
    # for eta in etas:
    #     # eta(z) = xi/4*(1/z-1) -1
    #     # (1+eta)(4/xi) = 1/z-1
    #     z = 1.0 / (1.0 + (1.0 + eta) * 4.0 / xi)
    #     norm = xi / np.pi / z  # * (4.0 * np.pi) ** 2
    #     etap = xi / 4 * (1 / z - 1) - 1
    #     # print(eta,etap)
    #     zs.append(z)
    #     exa = LeProHQ.cg1("FL", "VV", xi, eta)
    #     he = bfkl.cg1_LL_FL_VV(z, xi) / norm
    #     print(eta, he, bfkl.cg1_asy_LL_FL_VV(z,xi)/norm)
    #     exas.append(exa)
    #     ns.append(he)

    zs = us * zmax
    etas = []
    for z in zs:
        eta = xi / 4.0 * (1.0 / z - 1.0) - 1.0
        etas.append(eta)

        norm = xi / np.pi / z  # * (4.0 * np.pi) ** 2
        # fmt: off
        exa,thr,hsc = LeProHQ.cg0("FL", "VV", xi, eta),cg0t("FL", "VV", xi, eta),(bmsn.clg1am0_a0(z)) / norm
        # exa,thr,hsc = LeProHQ.cg1("FL", "VV", xi, eta),cg1t("FL", "VV", xi, eta),(bmsn.clg2am0_a0(z) + bmsn.clg2am0_aq(z) * lxi + bfkl.cg1_LL_FL_VV(z,xi) - bfkl.cg1_asy_LL_FL_VV(z,xi)) / norm
        # exa,thr,hsc = LeProHQ.cg0("F2", "VV", xi, eta),cg0t("F2", "VV", xi, eta),(bmsn.c2g1am0_a0(z) + bmsn.c2g1am0_aq(z)*lxi) / norm
        # exa,thr,hsc = LeProHQ.cg1("F2", "VV", xi, eta),cg1t("F2", "VV", xi, eta),(bmsn.c2g2am0_a0(z) + bmsn.c2g2am0_aq(z) * lxi + bmsn.c2g2am0_aq2(z) * lxi**2 + bfkl.cg1_LL_F2_VV(z,xi) - bfkl.cg1_asy_LL_F2_VV(z,xi) ) / norm
        # fmt: on
        A, B, C, D = 1.7, 2.5, 2.5, 1.2
        a = c = 2.5
        b = d = 5.0
        h = A + (B - A) / (1.0 + np.exp(a * (lxi - b)))
        k = C + (D - C) / (1.0 + np.exp(c * (lxi - d)))
        f1 = 1.0 / (1.0 + (eta / h) ** k)
        n = thr * f1 + hsc * (1.0 - f1)
        # print(exa,thr,f1, hsc, (1-f1),n)
        exas.append(exa)
        thrs.append(thr)
        hscs.append(hsc)
        f1s.append(f1)
        ns.append(n)

    exas = np.array(exas)
    # thrs = np.array(thrs)
    # hscs = np.array(hscs)
    # f1s = np.array(f1s)
    ns = np.array(ns)
    print("etas:", etas)
    # print("zs:", zs)
    print("ls:", exas)
    # print("ts:", thrs)
    # print("bs:", hscs)
    # print("f1s:", f1s)
    print("ns:", ns)
    # # print(a)
    # # print(b)
    print("all/N", exas / ns)
    # # print("all/thr",l/t)
    # # plt.plot(etas, a/b)
    # # plt.show()
    print()
