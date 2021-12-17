import pathlib

import matplotlib.pyplot as plt
import numpy as np
from LeProHQ.utils import datadir, load_2d_interpolation, raw_cb

path = pathlib.Path(__file__).parent / "bulk_border"


def check(parton, proj, cc):
    """Plot the bulk borders"""
    # load data
    cf = f"c{parton}1"
    grid_bulk, bulk_int = load_2d_interpolation(
        str(datadir) + f"/{cf}/{cf}-{proj}_{cc}-bulk.dat"
    )
    xis = np.exp(grid_bulk[0][1:])
    etas = np.exp(grid_bulk[1:, 0])
    # xi_min
    xi_min = [raw_cb(proj, cc, xis[0], eta, grid_bulk, bulk_int) for eta in etas]
    fig, ax = plt.subplots()
    fig.suptitle(f"$\\xi_{{min}}$ = {xis[0]}")
    ax.plot(np.log(etas), xi_min)
    ax.set_xlabel("$ln(\eta)$")
    fig.savefig(path / f"{parton}-{proj}-{cc}-xi_min.png")
    plt.close(fig)
    # xi_max
    xi_max = [raw_cb(proj, cc, xis[-1], eta, grid_bulk, bulk_int) for eta in etas]
    fig, ax = plt.subplots()
    fig.suptitle(f"$\\xi_{{max}}$ = {xis[-1]}")
    ax.set_xlabel("$ln(\eta)$")
    ax.plot(np.log(etas), xi_max)
    fig.savefig(path / f"{parton}-{proj}-{cc}-xi_max.png")
    plt.close(fig)
    # eta_min
    eta_min = [raw_cb(proj, cc, xi, etas[0], grid_bulk, bulk_int) for xi in xis]
    fig, ax = plt.subplots()
    fig.suptitle(f"$\\eta_{{min}}$ = {etas[0]}")
    ax.plot(np.log(xis), eta_min)
    ax.set_xlabel("$ln(\\xi)$")
    fig.savefig(path / f"{parton}-{proj}-{cc}-eta_min.png")
    plt.close(fig)
    # eta_max
    eta_max = [raw_cb(proj, cc, xi, etas[-1], grid_bulk, bulk_int) for xi in xis]
    fig, ax = plt.subplots()
    fig.suptitle(f"$\\eta_{{max}}$ = {etas[-1]}")
    ax.plot(np.log(xis), eta_max)
    ax.set_xlabel("$ln(\\xi)$")
    fig.savefig(path / f"{parton}-{proj}-{cc}-eta_max.png")
    plt.close(fig)


def check_all(parton):
    """Iterate all projections"""
    for c in ["V", "A"]:
        for p in ["F2", "FL", "x2g1"]:
            check(parton, p, f"{c}{c}")


check_all("g")
