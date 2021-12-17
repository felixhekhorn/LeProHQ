#!/usr/bin/env python
# -*- coding: utf-8 -*-
import abc
import pathlib

import matplotlib.pyplot as plt
import numpy as np

import config
import runner
from LeProHQpp import FullyDiffHistT
from LeProHQpp import Projection as proj

# global parameters
xTilde = 0.8
omega = 1.0
deltax = 1e-6
deltay = 7e-6

objArgs = (config.nlf, config.m2, xTilde, omega, deltax, deltay)

dir = pathlib.Path(__file__).parent / "FullyDiff"


class AbstractRunner(abc.ABC):
    """
    Generate the plots.

    Parameters
    ----------
        setups : list(tuple)
            configurations
        n_bins : int
            number of bins
    """

    def __init__(self, path, setups):
        self.r = runner.FullyDiffRunner()
        self.m = None
        self.path = path
        self.setups = setups

    def append_all(self):
        """Add all curves with all configurations"""
        for k, cfg in enumerate(self.setups):
            if cfg["orderFlag"] in config.flags:
                self.append(cfg, k)

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
                number of jobs
        """
        self.append_all()
        l = self.r.run(n_jobs)

    def search(self, fnc):
        vals = list(
            filter(
                fnc,
                self.setups,
            )
        )
        idc = []
        for v in vals:
            idx = self.setups.index(v)
            idc.append(idx)
        return idc


class Fig10Runner(AbstractRunner):
    """
    Generate Figure 10.

    Parameters
    ----------
        path : str
            histogram path
        n_bins : int
            number of bins
    """

    def __init__(self, path, n_bins):
        super().__init__(path, config.setup_fig10)
        self.n_bins = n_bins

    def append(self, cfg, tag):
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
                "run": "F",
                "activateHistograms": [
                    (
                        FullyDiffHistT.HAQTransverseMomentumScaling,
                        self.n_bins,
                        str(self.path) % (10, tag, "xt"),
                        # 1e-3,
                        # 1,
                    )
                ],
                "muF2": (1.0 * cfg["cMuF2"], 0.0, 1.0 * cfg["cMuF2"], 0.0),
                "muR2": (1.0 * cfg["cMuR2"], 0.0, 1.0 * cfg["cMuR2"], 0.0),
                "flags": cur_flags,
                "tag": tag,
                "msg": f"S_h = {cfg['hadronicS']}, muf2 = {cfg['cMuF2']}, mur2 = {cfg['cMuR2']}",
                "IntegrationOutput": True,
            }
        )

    def plot(self, dir, fp, show=False):
        """
        Draw figure.

        Parameters
        ----------
            fp : str
                template file name
            show : bool
                show plot on screen?
        """
        print("[INFO] Plotting ...")
        fig = plt.figure()
        fig.suptitle(f"Figure 10")
        ax = fig.add_subplot(111)
        ax.set_xlabel("$x_T$")
        ax.set_ylabel(r"$d\Delta\sigma_{\gamma p}^c/dx_T$ [nb]")
        ax.set_xscale("log")
        # low scale centrals
        central_lo_low = self.search(
            lambda e: np.isclose(e["hadronicS"], 10.0 ** 2)
            and np.isclose(e["cMuF2"], 1)
            and np.isclose(e["cMuR2"], 1)
            and e["orderFlag"] == "LO"
        )
        fn = str(dir) % (10, central_lo_low[0], "xt")
        if pathlib.Path(fn).exists():
            data = np.loadtxt(fn)
            ax.plot(
                0.5 * (data[:, 1] + data[:, 0]),
                data[:, 2]
                / (data[:, 1] - data[:, 0])
                * config.norm_sigma
                * config.nbTimesGeV2,
                # xerror=(data[:, 1] - data[:, 0]),
                # align="edge",
                label=r"LO, $\sqrt{S}=10$ GeV",
            )
        central_nlo_low = self.search(
            lambda e: np.isclose(e["hadronicS"], 10.0 ** 2)
            and np.isclose(e["cMuF2"], 1)
            and np.isclose(e["cMuR2"], 1)
            and e["orderFlag"] == "NLO"
        )
        fn = str(dir) % (10, central_nlo_low[0], "xt")
        if pathlib.Path(fn).exists():
            data = np.loadtxt(fn)
            ax.bar(
                data[:, 0],
                data[:, 2]
                / (data[:, 1] - data[:, 0])
                * config.norm_sigma
                * config.nbTimesGeV2,
                (data[:, 1] - data[:, 0]),
                align="edge",
                label=r"NLO, $\sqrt{S}=10$ GeV",
            )
        # low scale bands
        band_lo_low = self.search(
            lambda e: np.isclose(e["hadronicS"], 10.0 ** 2) and e["orderFlag"] == "LO"
        )
        lefts = np.array([])
        rights = np.array([])
        all_data = None
        for idx in band_lo_low:
            fn = str(dir) % (10, idx, "xt")
            if pathlib.Path(fn).exists():
                data = np.loadtxt(fn)
                lefts = data[:, 0]
                rights = data[:, 1]
                if all_data is None:
                    all_data = data[:, 2]
                else:
                    all_data = np.c_[all_data, data[:, 2]]
        ax.plot(
            0.5 * (rights + lefts),
            np.max(all_data, axis=1)
            / (rights - lefts)
            * config.norm_sigma
            * config.nbTimesGeV2,
            # xerror=(data[:, 1] - data[:, 0]),
            # align="edge",
        )
        ax.plot(
            0.5 * (rights + lefts),
            np.min(all_data, axis=1)
            / (rights - lefts)
            * config.norm_sigma
            * config.nbTimesGeV2,
            # xerror=(data[:, 1] - data[:, 0]),
            # align="edge",
        )
        # high scale centrals
        central_lo_high = self.search(
            lambda e: np.isclose(e["hadronicS"], 200.0 ** 2)
            and np.isclose(e["cMuF2"], 1)
            and np.isclose(e["cMuR2"], 1)
            and e["orderFlag"] == "LO",
        )
        fn = str(dir) % (10, central_lo_high[0], "xt")
        if pathlib.Path(fn).exists():
            data = np.loadtxt(fn)
            ax.bar(
                data[:, 0],
                data[:, 2]
                / (data[:, 1] - data[:, 0])
                * config.norm_sigma
                * config.nbTimesGeV2
                / 8.0,
                (data[:, 1] - data[:, 0]),
                align="edge",
                label=r"LO, $\sqrt{S}=200$ GeV",
            )
        central_nlo_high = self.search(
            lambda e: np.isclose(e["hadronicS"], 200.0 ** 2)
            and np.isclose(e["cMuF2"], 1)
            and np.isclose(e["cMuR2"], 1)
            and e["orderFlag"] == "NLO"
        )
        fn = str(dir) % (10, central_nlo_high[0], "xt")
        if pathlib.Path(fn).exists():
            data = np.loadtxt(fn)
            ax.bar(
                data[:, 0],
                data[:, 2]
                / (data[:, 1] - data[:, 0])
                * config.norm_sigma
                * config.nbTimesGeV2
                / 8.0,
                (data[:, 1] - data[:, 0]),
                align="edge",
                label=r"NLO, $\sqrt{S}=200$ GeV",
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
        path : str
            histogram path
        n_bins : int
            number of bins
    """

    def __init__(self, path, n_bins):
        super().__init__(path, config.setup_fig10)
        self.n_bins = n_bins

    def append(self, cfg, tag):
        cur_flags = config.default_flags.copy()
        cur_flags.update(config.flags[cfg["orderFlag"]])
        self.r.append(
            {
                "objArgs": objArgs,
                "projection": cfg["proj_"],
                "Q2": config.Q2,
                "hadronicS": 10.0 ** 2,
                "lambdaQCD": config.lambdaQCD,
                "pdf": (cfg["pdf"], 0),
                "IntegrationConfig": config.intCfg,
                "run": "F",
                "activateHistograms": [
                    (
                        FullyDiffHistT.HAQRapidity,
                        self.n_bins,
                        str(self.path) % (11, tag, "y"),
                        # 1e-3,
                        # 1,
                    )
                ],
                "mu2": (2.0, 0.0, 0.0, 0.0),
                "flags": cur_flags,
                "tag": tag,
                "msg": f"S_h = {cfg['hadronicS']}, muf2 = {cfg['cMuF2']}, mur2 = {cfg['cMuR2']}",
                "IntegrationOutput": True,
            }
        )

    def plot(self, dir, fp, show=False):
        """
        Draw figure.

        Parameters
        ----------
            fp : str
                template file name
            show : bool
                show plot on screen?
        """
        print("[INFO] Plotting ...")
        # fig = plt.figure()
        # fig.suptitle(f"Figure 10")
        # ax = fig.add_subplot(111)
        # ax.set_xlabel("$x_T$")
        # ax.set_ylabel(r"$d\Delta\sigma_{\gamma p}^c/dx_T$ [nb]")
        # ax.set_xscale("log")
        # fig.legend()
        # fig.savefig(fp)
        # if show:
        #     plt.show()
        # plt.close(fig)


datadir = dir / "fig%d-%d-%s.dat"

# Figure 10
n_xt_bins = 20
fig10 = Fig10Runner(datadir, n_xt_bins)
# fig10.run(1)
fig10.plot(datadir, dir / "fig10.pdf")
