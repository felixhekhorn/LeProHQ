#!/usr/bin/env python
# -*- coding: utf-8 -*-
import abc
import pathlib

import matplotlib.pyplot as plt
import numpy as np
import yaml

import plots
import runner
from LeProHQpp import FullyDiffHistT
from LeProHQpp import Projection as proj

# global parameters
nlf = 3
m2 = 1.5 ** 2
xTilde = 0.8
omega = 1.0
deltax = 1e-6
deltay = 7e-6
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
n_rap_bins = 50
# pt configuration
n_pt_bins = 50
mu2_pt = (
    4.0,
    1.0,
    4.0,
    0.0,
)

objArgs = (nlf, m2, xTilde, omega, deltax, deltay)

dir = pathlib.Path(__file__).parent / "FullyDiff"


class AbstractRunner(abc.ABC):
    """
    Generate the plots.

    Parameters
    ----------
        n_bins : int
            number of bins
    """

    def __init__(self):
        self.r = runner.FullyDiffRunner()
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


class RapidityRunner(AbstractRunner):
    """
    Generate the rapidity plots.

    Parameters
    ----------
        path : str
            histogram path
        n_bins : int
            number of bins
    """

    setups = plots.y

    def __init__(self, path, n_bins):
        super().__init__()
        self.path = path
        self.n_bins = n_bins

    def append(self, cfg, flag_label, tag):
        xBj = cfg["x"]
        y0 = cfg["y0"]
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
                "run": "F",
                "activateHistograms": [
                    (
                        FullyDiffHistT.HAQRapidity,
                        self.n_bins,
                        str(self.path) % (cfg["num"], flag_label),
                        -y0,
                        y0,
                    )
                ],
                "flags": cur_flags,
                "tag": tag,
                "msg": f"x={xBj}, Q2={Q2}, y, {flag_label}",
                "IntegrationOutput": True,
            }
        )

    def plot(self, dir, fp, show=False):
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
            fn = str(fp) % cfg["num"]
            fig = plt.figure()
            fig.suptitle(f"x = {cfg['x']:g}")
            ax = fig.add_subplot(111)
            ax.set_xlabel("y (rapidity)")
            kind = "2" if cfg["proj_"] == proj.F2 else "L"
            ax.set_ylabel(f"$dF_{kind}(x,Q^2,m_c^2,y)/dy$")
            ax.set_yscale("log")
            for flag_label in flags:
                data = np.loadtxt(str(dir) % (cfg["num"], flag_label))
                # y_paper = -y_me
                ax.bar(
                    -data[:, 0],
                    data[:, 2] / (data[:, 1] - data[:, 0]),
                    -(data[:, 1] - data[:, 0]),
                    align="edge",
                    label=flag_label,
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
        path : str
            histogram path
        n_bins : int
            number of bins
    """

    setups = plots.pt

    def __init__(self, path, n_bins):
        super().__init__()
        self.path = path
        self.n_bins = n_bins

    def append(self, cfg, flag_label, tag):
        xBj = cfg["x"]
        ptmax = cfg["ptmax"]
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
                "run": "F",
                "activateHistograms": [
                    (
                        FullyDiffHistT.HAQTransverseMomentum,
                        self.n_bins,
                        str(self.path) % (cfg["num"], flag_label),
                        0,
                        ptmax,
                    )
                ],
                "flags": cur_flags,
                "mu2": mu2_pt,
                "msg": f"x={xBj}, Q2={Q2}, pt, {flag_label}",
                "IntegrationOutput": True,
            }
        )

    def plot(self, dir, fp, show=False):
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
            fn = str(fp) % cfg["num"]
            fig = plt.figure()
            fig.suptitle(f"x = {cfg['x']:g}")
            ax = fig.add_subplot(111)
            ax.set_xlabel("$p_t$ (GeV/c)")
            kind = "2" if cfg["proj_"] == proj.F2 else "L"
            ax.set_ylabel(f"$dF_{kind}(x,Q^2,m_c^2,p_t)/dp_t$")
            ax.set_yscale("log")
            for flag_label in flags:
                data = np.loadtxt(str(dir) % (cfg["num"], flag_label))
                ax.bar(
                    data[:, 0],
                    data[:, 2] / (data[:, 1] - data[:, 0]),
                    (data[:, 1] - data[:, 0]),
                    align="edge",
                    label=flag_label,
                )
            ax.set_ylim(cfg["yrange"])
            ax.set_xlim([0, cfg["ptmax"]])
            ax.tick_params(bottom=True, top=True, left=True, right=True, which="both")
            fig.legend()
            fig.savefig(fn)
            if show:
                plt.show()
            plt.close(fig)


datadir = dir / "fig%d-%s.dat"
# yr = RapidityRunner(datadir, n_rap_bins)
# yr.run(-2)
# yr.plot(datadir, dir / "fig%d.pdf")

ptr = TransverseMomentumRunner(datadir, n_pt_bins)
# ptr.run(-2)
ptr.plot(datadir, dir / "fig%d.pdf")
