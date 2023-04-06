import numpy as np
import matplotlib.pyplot as plt

from LeProHQ.cg1_ import cg1hv
from LeProHQ.cq1_ import cq1hv
from LeProHQ import cg1, cq1

etas = np.geomspace(1e-1, 1e4, 50)

def plot_cg1(proj, cc, xi = 1e3):
    full = np.array([cg1(proj,cc,xi,eta) for eta in etas])
    hv = np.array([cg1hv(proj,cc,xi,eta) for eta in etas])
    plot(proj, cc, xi, full, hv, "g")

def plot_cq1(proj, cc, xi = 1e3):
    full = np.array([cq1(proj,cc,xi,eta) for eta in etas])
    hv = np.array([cq1hv(proj,cc,xi,eta) for eta in etas])
    plot(proj, cc, xi, full, hv, "q")

def plot(proj, cc, xi, full, hv, parton):
    fig, ax = plt.subplots(2,1,sharex=True)
    fig.suptitle(f"$c_{parton}^1(\\eta,\\xi={xi})$ for {proj}_{cc}")
    ax[0].plot(etas, full, label="exact")
    ax[0].plot(etas, hv, label="limit")
    ax[0].set_xscale("log")
    ax[1].plot(etas, np.abs(hv/full - 1.))
    ax[1].set_yscale("log")
    ax[1].set_xlabel("$\\eta$")
    fig.legend()
    fig.savefig(f"c{parton}1-{proj}_{cc}-{xi}.pdf")

for proj in ["F2","FL"]:
    for cc in ["VV", "AA"]:
        plot_cg1(proj, cc)
        plot_cq1(proj, cc)
