#!/usr/bin/python3
# -*- coding: utf-8 -*-
import abc

import matplotlib.pyplot as plt
import numpy as np
import runner
import yaml
from LeProHQpp import Projection as proj

# global parameters
nlf = 3
m2 = 1.5 ** 2
Delta = 1e-6
Q2 = 10
lambdaQCD = 0.194
pdf = ("MorfinTungB", 0)
intCfg = {"verbosity": 1}
default_flags = {"usePhoton": True, "usePhotonZ": False, "useZ": False}
flags = {  # split by different contributions
    "LO": {"useLeadingOrder": True, "useNextToLeadingOrder": False},  # LO result
    # "NLO": {"useLeadingOrder": True , "useNextToLeadingOrder": True}, # NLO result
}

# rapidity configuration
n_rap_bins = 20
# pt configuration
n_pt_bins = 20
mu2_pt = (
    4.0,
    1.0,
    4.0,
    0.0,
)

objArgs = (nlf, m2, Delta)


class AbstractRunner(abc.ABC):
    """
    Generate the plots.

    Parameters
    ----------
        n_bins : int
            number of bins
    """

    def __init__(self):
        self.r = runner.InclusiveRunner()
        self.m = None

    def append_all(self):
        """Add all curves with all configurations"""
        for k, cfg in enumerate(self.setups):
            for l, flag_label in enumerate(flags.keys()):
                self.append(cfg, flag_label, (k, l))

    @abc.abstractmethod
    def append(self, cfg, flag_label, tag):
        """
        Add a single curve

        Parameters
        ----------
            cfg : dict
                configuration
            flag_label : str
                flags identifier
            tag : tuple
                config identifier
        """

    def run(self, n_jobs=None):
        """
        Compute the data.

        Parameters
        ----------
            n_jobs : int
                number ofjobs
        """
        self.append_all()
        l = self.r.run(n_jobs)
        self.m = self.reorder(l)

    def reorder(self, l):
        """
        Map linear output into a matrix

        Parameters
        ----------
            l : list(dict)
                linear output

        Returns
        -------
            m : np.ndarray
                nested data accoring to tag and bin
        """
        m_raw = np.full((len(self.setups), len(flags), len(self.bins)), np.nan)
        m = {"res": m_raw.copy(), "error": m_raw.copy()}
        for e in l:
            m["res"][(*e["tag"], e["bin"])] = e["res"]
            m["error"][(*e["tag"], e["bin"])] = e["IntegrationOutput"].error
        return m

    def dump(self, fn):
        """
        Write data to file.

        Parameters
        ----------
            fn : str
                file name
        """
        print(f"[INFO] Dumping to {fn} ...")
        s = []
        for e in self.setups:
            d = e.copy()
            d["proj_"] = str(d["proj_"])
            d["yrange"] = list(e["yrange"])
            s.append(d)
        out = {
            "sorting": ["setups", "flags", "bins"],
            "setups": s,
            "flags": [f for f in flags],
            "bins": self.bins.tolist(),
            "res": self.m["res"].tobytes(),
            "error": self.m["error"].tobytes(),
        }
        with open(fn, "w") as o:
            yaml.safe_dump(out, o)

    def load(self, fn):
        """
        Read data from file.

        Parameters
        ----------
            fn : str
                file name
        """
        print(f"[INFO] Loading {fn} ...")
        with open(fn, "r") as o:
            mm = yaml.safe_load(o)
        self.m = {}
        shape = (
            len(mm["setups"]),
            len(mm["flags"]),
            len(mm["bins"]),
        )
        self.m["res"] = np.frombuffer(mm["res"]).reshape(shape)
        self.m["error"] = np.frombuffer(mm["error"]).reshape(shape)


class RapidityRunner(AbstractRunner):
    """
    Generate the rapidity plots.

    Parameters
    ----------
        n_bins : int
            number of bins
    """

    setups = (
        dict(num=11, proj_=proj.F2, x=0.1, y0=2.0, yrange=(1e-6, 1e-2)),
        dict(num=12, proj_=proj.F2, x=0.01, y0=3.5, yrange=(1e-5, 1e-1)),
        dict(num=13, proj_=proj.F2, x=0.001, y0=4.5, yrange=(1e-5, 1e-1)),
        dict(num=14, proj_=proj.F2, x=0.0001, y0=5.5, yrange=(1e-5, 1e-1)),
        dict(num=16, proj_=proj.FL, x=0.1, y0=2.0, yrange=(1e-7, 1e-3)),
        dict(num=17, proj_=proj.FL, x=0.01, y0=3.0, yrange=(1e-6, 1e-2)),
        dict(num=18, proj_=proj.FL, x=0.001, y0=4.5, yrange=(1e-6, 1e-2)),
        dict(num=19, proj_=proj.FL, x=0.0001, y0=5.5, yrange=(1e-5, 1e-1)),
    )
    """list of configurations for each figure"""

    def __init__(self, n_bins):
        super().__init__()
        self.bins = np.linspace(-1.0, 1.0, n_bins)

    def append(self, cfg, flag_label, tag):
        xBj = cfg["x"]
        y0 = cfg["y0"]
        for bin_, y in enumerate(y0 * self.bins):
            cur_flags = default_flags.copy()
            cur_flags.update(flags[flag_label])
            self.r.append(
                {
                    "objArgs": objArgs,
                    "projection": cfg["proj_"],
                    "Q2": Q2,
                    "xBjorken": xBj,
                    "lambdaQCD": lambdaQCD,
                    "pdf": pdf,
                    "IntegrationConfig": intCfg,
                    "run": ("dF_dHAQRapidity", y),
                    "flags": cur_flags,
                    "bin": bin_,
                    "tag": tag,
                    "msg": f"x={xBj}, Q2={Q2}, y={y}",
                    "IntegrationOutput": True,
                }
            )

    def plot(self, fp="fig%d.pdf", show=False):
        """
        Draw all plots

        Parameters
        ----------
            fp : str
                template file name
            show : bool
                show plot on screen?
        """
        print("[INFO] Plotting ...")
        for k, cfg in enumerate(self.setups):
            fn = fp % cfg["num"]
            fig = plt.figure()
            fig.suptitle(f"x = {cfg['x']:g}")
            ax = fig.add_subplot(111)
            ax.set_xlabel("y (rapidity)")
            kind = "2" if cfg["proj_"] == proj.F2 else "L"
            ax.set_ylabel(f"$dF_{kind}(x,Q^2,m_c^2,y)/dy$")
            # y_paper = -y_me
            ax.semilogy(
                -self.bins * cfg["y0"],
                self.m["res"][k, 0],
                label="LO",
            )
            if "NLO" in flags:
                ax.semilogy(
                    -self.bins * cfg["y0"],
                    self.m["res"][k, 1],
                    label="NLO",
                )
            ax.set_ylim(cfg["yrange"])
            ax.set_xlim([-cfg["y0"], cfg["y0"]])
            ax.tick_params(bottom=True, top=True, left=True, right=True, which="both")
            fig.legend()
            fig.savefig(fn)
            if show:
                plt.show()
            plt.close(fig)


class TransverseMomentumRunner(AbstractRunner):
    """
    Generate the transverse momentum plots.

    Parameters
    ----------
        n_bins : int
            number of bins
    """

    setups = (
        dict(num=1, proj_=proj.F2, x=0.1, ptmax=5.0, yrange=(1e-5, 1e-2)),
        dict(num=2, proj_=proj.F2, x=0.01, ptmax=8.0, yrange=(1e-4, 1e-1)),
        dict(num=3, proj_=proj.F2, x=0.001, ptmax=15.0, yrange=(1e-4, 1e-1)),
        dict(num=4, proj_=proj.F2, x=0.0001, ptmax=20.0, yrange=(1e-4, 1e-1)),
        dict(num=6, proj_=proj.FL, x=0.1, ptmax=5.0, yrange=(1e-6, 1e-3)),
        dict(num=7, proj_=proj.FL, x=0.01, ptmax=10.0, yrange=(1e-5, 1e-2)),
        dict(num=8, proj_=proj.FL, x=0.001, ptmax=15.0, yrange=(1e-5, 1e-2)),
        dict(num=9, proj_=proj.FL, x=0.0001, ptmax=20.0, yrange=(1e-5, 1e-2)),
    )
    """list of configurations for each figure"""

    def __init__(self, n_bins):
        super().__init__()
        self.bins = np.linspace(0.0, 1.0, n_bins)

    def append(self, cfg, flag_label, tag):
        xBj = cfg["x"]
        ptmax = cfg["ptmax"]
        for bin_, pt in enumerate(ptmax * self.bins):
            cur_flags = default_flags.copy()
            cur_flags.update(flags[flag_label])
            self.r.append(
                {
                    "objArgs": objArgs,
                    "projection": cfg["proj_"],
                    "Q2": Q2,
                    "xBjorken": xBj,
                    "lambdaQCD": lambdaQCD,
                    "pdf": pdf,
                    "IntegrationConfig": intCfg,
                    "run": ("dF_dHAQTransverseMomentum", pt),
                    "flags": cur_flags,
                    "bin": bin_,
                    "mu2": mu2_pt,
                    "tag": tag,
                    "msg": f"x={xBj}, Q2={Q2}, pt={pt}",
                    "IntegrationOutput": True,
                }
            )

    def plot(self, fp="fig%d.pdf", show=False):
        """
        Draw all plots

        Parameters
        ----------
            fp : str
                template file name
            show : bool
                show plot on screen?
        """
        print("[INFO] Plotting ...")
        for k, cfg in enumerate(self.setups):
            fn = fp % cfg["num"]
            fig = plt.figure()
            fig.suptitle(f"x = {cfg['x']:g}")
            ax = fig.add_subplot(111)
            ax.set_xlabel("$p_t$ (GeV/c)")
            kind = "2" if cfg["proj_"] == proj.F2 else "L"
            ax.set_ylabel(f"$dF_{kind}(x,Q^2,m_c^2,p_t)/dp_t$")
            ax.semilogy(
                self.bins * cfg["ptmax"],
                self.m["res"][k, 0],
                label="LO",
            )
            if "NLO" in flags:
                ax.semilogy(
                    self.bins * cfg["ptmax"],
                    self.m["res"][k, 1],
                    label="NLO",
                )
            ax.set_ylim(cfg["yrange"])
            ax.set_xlim([0, cfg["ptmax"]])
            ax.tick_params(bottom=True, top=True, left=True, right=True, which="both")
            fig.legend()
            fig.savefig(fn)
            if show:
                plt.show()
            plt.close(fig)


# rr = RapidityRunner(n_rap_bins)
# rr.run(-2)
# # rr.dump("rap-lo.yaml")
# # rr.load("rap-lo.yaml")
# rr.plot()

pr = TransverseMomentumRunner(n_pt_bins)
pr.run(-2)
# pr.dump("pt.yaml")
# pr.load("pt-lo.yaml")
pr.plot()
