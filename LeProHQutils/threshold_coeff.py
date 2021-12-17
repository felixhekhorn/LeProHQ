import pathlib

import matplotlib.pyplot as plt
import numpy as np
from LeProHQ.cg1 import cg1t
from LeProHQ.cq1 import cq1t
from scipy.optimize import minimize_scalar

path = pathlib.Path(__file__).parent / "threshold_data/"

parton_map = {
    "q": cq1t,
    "g": cg1t,
}


def cp1tp(parton, proj, cc, xi, eta, a):
    return parton_map[parton](proj, cc, xi, eta) * (1.0 + a * eta)


def get_ker(parton, proj, cc, xi, etas, fs):
    """Returns fitting kernel"""

    def ker(a):
        gs = [cp1tp(parton, proj, cc, xi, eta, a) for eta in etas]
        return np.sum((fs - gs) ** 2)

    return ker


def fit_thres(
    parton,
    proj,
    cc,
    show_plot=False,
    show_eta_plot=False,
    save_plot=False,
    save_eta_plot=False,
):
    """Run the fitting procedure"""
    # read from file
    fp_in = "c%s1-%s_%s-thres.dat" % (
        parton,
        proj,
        cc,
    )
    print(f"loading input {fp_in}")
    m = np.loadtxt(path / fp_in)
    xis = np.exp(m[0][1:])
    etas = np.exp(m[1:, 0])
    data = m[1:, 1:].T
    # run fit
    cs = []
    for j, xi in enumerate(xis):
        res = minimize_scalar(get_ker(parton, proj, cc, xi, etas, data[j]))
        if not res.success:
            print("Problem:", parton, proj, cc, j, res, "FAILED!!!")
            continue
        # show the user
        if res.success and (show_eta_plot or save_eta_plot):
            fig, ax = plt.subplots()
            fig.suptitle(f"xi = {xi}, a = {res.x}")
            ax.plot(np.log(etas), data[j], "x")
            # plot approximation and raw
            etas2 = np.geomspace(min(etas), max(etas), 50)
            ax.plot(
                np.log(etas2),
                [cp1tp(parton, proj, cc, xi, eta, 0) for eta in etas2],
                "b",
            )
            ax.plot(
                np.log(etas2),
                [cp1tp(parton, proj, cc, xi, eta, res.x) for eta in etas2],
                "r",
            )
            if show_eta_plot:
                fig.show()
            if save_eta_plot:
                fig.savefig(path / "img" / f"{parton}-{proj}-{cc}-{j}.png")
            plt.close(fig)
        # save for me
        print(parton, proj, cc, xi, res.x, res.fun)
        cs.append(res.x)
    # write back
    fp_out = "c%s1-%s_%s-thres-coeff.dat" % (parton, proj, cc)
    print(f"write output {fp_out}")
    np.savetxt(path / fp_out, np.array([np.log(xis), cs]), fmt="%.5e")
    if show_plot or save_plot:
        fig, ax = plt.subplots()
        fig.suptitle(fp_in)
        ax.plot(np.log(xis), cs)
        if show_plot:
            fig.show()
        if save_plot:
            fig.savefig(path / "img" / f"{parton}-{proj}-{cc}.png")
        plt.close(fig)


def fit_all(
    parton, show_plot=False, show_eta_plot=False, save_plot=True, save_eta_plot=False
):
    """Iterate all projections"""
    for c in ["V", "A"]:
        for p in ["F2", "FL", "x2g1"]:
            fit_thres(
                parton,
                p,
                f"{c}{c}",
                show_plot=show_plot,
                show_eta_plot=show_eta_plot,
                save_plot=save_plot,
                save_eta_plot=save_eta_plot,
            )


fit_all("g")
