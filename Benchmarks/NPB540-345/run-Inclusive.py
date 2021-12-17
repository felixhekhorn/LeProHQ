#!/usr/bin/env python
# -*- coding: utf-8 -*-
import abc
import pathlib

import matplotlib.pyplot as plt
import numpy as np
import yaml

import config
import runner
from LeProHQpp import Projection as proj

# global parameters
Delta = 1e-6

objArgs = (config.nlf, config.m2, Delta)

dir = pathlib.Path(__file__).parent / "Inclusive"


class AbstractRunner(abc.ABC):
    """
    Generate the plots.

    Parameters
    ----------
        setups : list(tuple)
            configurations
        bins : numpy.ndarray
            bins
    """

    def __init__(self, setups, bins):
        self.r = runner.InclusiveRunner()
        self.m = None
        self.setups = setups
        self.bins = bins

    def append_all(self):
        """Add all curves with all configurations"""
        for k, cfg in enumerate(self.setups):
            if cfg["orderFlag"] in config.flags:
                self.append(cfg, k)

    @abc.abstractmethod
    def append(self, cfg, tag):
        """
        Add a single curve

        Parameters
        ----------
            cfg : dict
                configuration
            tag : tuple
                config identifier
        """

    def run(self, n_jobs=None):
        """
        Compute the data.

        Parameters
        ----------
            n_jobs : int
                number of jobs
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
        m_raw = np.full((len(self.setups), len(self.bins)), np.nan)
        m = {"res": m_raw.copy(), "error": m_raw.copy()}
        for e in l:
            m["res"][(e["tag"], e["bin"])] = e["res"]
            m["error"][(e["tag"], e["bin"])] = e["IntegrationOutput"].error
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
            if "proj_" in d:
                d["proj_"] = str(d["proj_"])
            s.append(d)
        out = {
            "sorting": ["setups", "bins"],
            "setups": s,
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
            len(mm["bins"]),
        )
        self.m["res"] = np.frombuffer(mm["res"]).reshape(shape)
        self.m["error"] = np.frombuffer(mm["error"]).reshape(shape)


class Fig10Runner(AbstractRunner):
    """
    Generate Figure 10.

    Parameters
    ----------
        n_bins : int
            number of bins
    """

    def __init__(self, n_bins):
        super().__init__(config.setup_fig10, np.geomspace(1e-3, 1.0, n_bins))

    def append(self, cfg, tag):
        for bin_, xt in enumerate(self.bins):
            cur_flags = config.default_flags.copy()
            cur_flags.update(config.flags[cfg["orderFlag"]])
            self.r.append(
                {
                    "objArgs": objArgs,
                    "projection": proj.x2g1,
                    "Q2": config.Q2,
                    "hadronicS": cfg["hadronicS"],
                    "lambdaQCD": config.lambdaQCD,
                    "pdf": (cfg["pdf"], 0),
                    "IntegrationConfig": config.intCfg,
                    "run": ("dF_dHAQTransverseMomentumScaling", xt),
                    "muF2": (1.0 * cfg["cMuF2"], 0.0, 1.0 * cfg["cMuF2"], 0.0),
                    "muR2": (1.0 * cfg["cMuR2"], 0.0, 1.0 * cfg["cMuR2"], 0.0),
                    "flags": cur_flags,
                    "bin": bin_,
                    "tag": tag,
                    "msg": f"S_h = {cfg['hadronicS']},  xt={xt}, muf2 = {cfg['cMuF2']}, mur2 = {cfg['cMuR2']}",
                    "IntegrationOutput": True,
                }
            )

    def search(self, fnc):
        vals = list(
            filter(
                fnc,
                self.setups,
            )
        )
        res = []
        for v in vals:
            idx = self.setups.index(v)
            res.append(self.m["res"][idx])
        return res

    def plot(self, fp, show=False):
        """
        Draw all plots

        Parameters
        ----------
            fp : str
                file name
            show : bool
                show plot on screen?
        """
        print("[INFO] Plotting ...")
        fig = plt.figure()
        fig.suptitle(f"Figure 10")
        ax = fig.add_subplot(111)
        ax.set_xlabel("$x_T$")
        ax.set_ylabel(r"$d\Delta\sigma_{\gamma p}^c/dx_T$ [nb]")
        # low scale centrals
        central_lo_low = self.search(
            lambda e: np.isclose(e["hadronicS"], 10.0 ** 2)
            and np.isclose(e["cMuF2"], 1)
            and np.isclose(e["cMuR2"], 1)
            and e["orderFlag"] == "LO"
        )
        ax.semilogx(
            self.bins,
            central_lo_low[0] * config.norm_sigma * config.nbTimesGeV2,
            label=r"LO, $\sqrt{S}=10$ GeV",
        )
        central_nlo_low = self.search(
            lambda e: np.isclose(e["hadronicS"], 10.0 ** 2)
            and np.isclose(e["cMuF2"], 1)
            and np.isclose(e["cMuR2"], 1)
            and e["orderFlag"] == "NLO"
        )
        ax.semilogx(
            self.bins,
            (central_nlo_low[0] * config.norm_sigma * config.nbTimesGeV2),
            label="NLO, $\sqrt{S}=10$ GeV",
        )
        # low scale bands
        band_lo_low = self.search(
            lambda e: np.isclose(e["hadronicS"], 10.0 ** 2) and e["orderFlag"] == "LO"
        )
        ax.semilogx(
            self.bins,
            np.max(np.array(band_lo_low), axis=0)
            * config.norm_sigma
            * config.nbTimesGeV2,
        )
        ax.semilogx(
            self.bins,
            np.min(np.array(band_lo_low), axis=0)
            * config.norm_sigma
            * config.nbTimesGeV2,
        )
        # high scale centrals
        central_lo_high = self.search(
            lambda e: np.isclose(e["hadronicS"], 200.0 ** 2)
            and np.isclose(e["cMuF2"], 1)
            and np.isclose(e["cMuR2"], 1)
            and e["orderFlag"] == "LO",
        )
        ax.semilogx(
            self.bins,
            (central_lo_high[0] * config.norm_sigma * config.nbTimesGeV2 / 8.0),
            label="LO, $\sqrt{S}=200$ GeV",
            ls="dashed",
        )
        central_nlo_high = self.search(
            lambda e: np.isclose(e["hadronicS"], 200.0 ** 2)
            and np.isclose(e["cMuF2"], 1)
            and np.isclose(e["cMuR2"], 1)
            and e["orderFlag"] == "NLO"
        )
        ax.semilogx(
            self.bins,
            (central_nlo_high[0] * config.norm_sigma * config.nbTimesGeV2 / 8.0),
            label="NLO, $\sqrt{S}=200$ GeV",
            ls="dashed",
        )
        ax.set_ylim(-70, 225)
        ax.set_xlim(1e-3, 1.0)
        ax.tick_params(bottom=True, top=True, left=True, right=True, which="both")
        fig.legend()
        fig.savefig(fp)
        if show:
            plt.show()
        plt.close(fig)


class Fig11Runner(AbstractRunner):
    """
    Generate Figure 11.

    Parameters
    ----------
        n_bins : int
            number of bins
    """

    def __init__(self, n_bins):
        super().__init__(config.setup_fig11_fig12, np.linspace(-2.0, 2.0, n_bins))

    def append(self, cfg, tag):
        for bin_, y in enumerate(self.bins):
            cur_flags = config.default_flags.copy()
            cur_flags.update(config.flags[cfg["orderFlag"]])
            self.r.append(
                {
                    "objArgs": objArgs,
                    "projection": cfg["proj_"],
                    "Q2": config.Q2,
                    "hadronicS": 10 ** 2,
                    "lambdaQCD": config.lambdaQCD,
                    "pdf": (cfg["pdf"], 0),
                    "IntegrationConfig": config.intCfg,
                    "run": ("dF_dHAQRapidity", y),
                    "muF2": (2.0, 0.0, 0.0, 0.0),
                    "muR2": (2.0, 0.0, 0.0, 0.0),
                    "flags": cur_flags,
                    "bin": bin_,
                    "tag": tag,
                    "msg": f"proj = {cfg['proj_']}, pdf = {cfg['pdf']}, y = {y}",
                    "IntegrationOutput": True,
                }
            )

    def plot(self, fp, show=False):
        """
        Draw all plots

        Parameters
        ----------
            fp : str
                file name
            show : bool
                show plot on screen?
        """
        print("[INFO] Plotting ...")
        fig = plt.figure()
        fig.suptitle(f"Figure 11")
        ax = fig.add_subplot(111)
        ax.set_xlabel("y")
        ax.set_ylabel(r"$d\Delta\sigma_{\gamma p}^c/dy$ [nb]")
        ax.set_ylim(-2.5, 27.5)
        ax.set_xlim(-2.0, 2.0)
        ax.tick_params(bottom=True, top=True, left=True, right=True, which="both")
        ax.axhline(y=0, color="black", lw=0.2)
        ax.plot(
            self.bins,
            (self.m["res"][0] * config.norm_sigma * config.nbTimesGeV2),
            label="LO",
        )
        ax.plot(
            self.bins,
            (self.m["res"][1] * config.norm_sigma * config.nbTimesGeV2),
            label="NLO",
        )
        ax.plot(
            self.bins,
            (self.m["res"][2] * config.norm_sigma * config.nbTimesGeV2 / 8.0),
            label=r"$d\sigma_{\gamma p}^c/dy/8$",
        )
        fig.legend()
        fig.savefig(fp)
        if show:
            plt.show()
        plt.close(fig)


class Fig12Runner(AbstractRunner):
    """
    Generate Figure 12.

    Parameters
    ----------
        n_bins : int
            number of bins
    """

    def __init__(self, n_bins):
        super().__init__(config.setup_fig11_fig12, np.linspace(-0.99, 0.99, n_bins))

    def append(self, cfg, tag):
        for bin_, xF in enumerate(self.bins):
            cur_flags = config.default_flags.copy()
            cur_flags.update(config.flags[cfg["orderFlag"]])
            self.r.append(
                {
                    "objArgs": objArgs,
                    "projection": cfg["proj_"],
                    "Q2": config.Q2,
                    "hadronicS": 10 ** 2,
                    "lambdaQCD": config.lambdaQCD,
                    "pdf": (cfg["pdf"], 0),
                    "IntegrationConfig": config.intCfg,
                    "run": ("dF_dHAQFeynmanX", xF),
                    "mu2": (2.0, 0.0, 0.0, 0.0),
                    "flags": cur_flags,
                    "bin": bin_,
                    "tag": tag,
                    "msg": f"proj = {cfg['proj_']}, pdf = {cfg['pdf']}, xF = {xF}",
                    "IntegrationOutput": True,
                }
            )

    def plot(self, fp, show=False):
        """
        Draw all plots

        Parameters
        ----------
            fp : str
                file name
            show : bool
                show plot on screen?
        """
        print("[INFO] Plotting ...")
        fig = plt.figure()
        fig.suptitle(f"Figure 12")
        ax = fig.add_subplot(111)
        ax.set_xlabel("$x_F$")
        ax.set_ylabel(r"$d\Delta\sigma_{\gamma p}^c/dx_F$ [nb]")
        ax.set_ylim(-3.0, 23.0)
        ax.set_xlim(-1.0, 1.0)
        ax.tick_params(bottom=True, top=True, left=True, right=True, which="both")
        ax.axhline(y=0, color="black", lw=0.2)
        ax.plot(
            self.bins,
            (self.m["res"][0] * config.norm_sigma * config.nbTimesGeV2),
            label="LO",
        )
        ax.plot(
            self.bins,
            (self.m["res"][1] * config.norm_sigma * config.nbTimesGeV2),
            label="NLO",
        )
        ax.plot(
            self.bins,
            (self.m["res"][2] * config.norm_sigma * config.nbTimesGeV2 / 11.0),
            label=r"$d\sigma_{\gamma p}^c/dx_F/11$",
        )
        fig.legend()
        fig.savefig(fp)
        if show:
            plt.show()
        plt.close(fig)


# Figure 10
# n_xt_bins = 20
# fig10 = Fig10Runner(n_xt_bins)
# fig10.run(1)
# fig10.dump(dir / "fig10.yaml")
# #fig10.load(dir / "fig10.yaml")
# fig10.plot(dir / "fig10.pdf")

# Figure 11
# n_y_bins = 20
# fig11 = Fig11Runner(n_y_bins)
# #fig11.run(1)
# #fig11.dump(dir / "fig11.yaml")
# fig11.load(dir / "fig11.yaml")
# fig11.plot(dir / "fig11.pdf")

# Figure 12
n_xF_bins = 20
fig12 = Fig12Runner(n_xF_bins)
fig12.run(1)
# fig12.dump(dir / "fig12.yaml")
# fig12.load(dir / "fig12.yaml")
fig12.plot(dir / "fig12.pdf")
