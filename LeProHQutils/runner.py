# -*- coding: utf-8 -*-

from multiprocessing import Pool

from LeProHQpp import (
    DynamicScaleFactors,
    FullyDiffLeptoProduction,
    InclusiveLeptoProduction,
)


class Runner:
    """Wrapper to run several jobs in parallel"""

    def __init__(self, callee):
        """Constructor"""
        self.args = []
        self.callee = callee
        self.output = []

    def clear(self):
        """
        Clear state
        """
        self.args.clear()
        self.output.clear()

    def append(self, e):
        """
        Append an element
        @param e element
        """
        self.args.append(e)

    def run(self, n_jobs=None):
        """
        Compute all elements
        @param n_jobs number of parallel threads
        """
        print("Computing %d points ..." % len(self.args))
        with Pool(n_jobs) as p:
            self.output = p.map(self.callee, self.args)
        return self.output


def global_setter(o, p):
    """Set global vars"""
    if "run" not in p:
        raise KeyError("no called function set")
    # global getter and setter
    if "projection" in p:
        o.setProjection(p["projection"])
    if "nlf" in p:
        o.setNumberOfLightFlavours(p["nlf"])
    if "Q2" in p:
        o.setQ2(p["Q2"])
    if "m2" in p:
        o.setM2(p["m2"])
    # partonic setter
    if "partonicEta" in p:
        o.setPartonicEta(p["partonicEta"])
    if "partonicS" in p:
        o.setPartonicS(p["partonicS"])
    # hadronic setter
    if "pdf" in p:
        o.setPdf(*p["pdf"])
    if "muR2" in p:
        o.setMuR2(DynamicScaleFactors(*p["muR2"]))
    if "muF2" in p:
        o.setMuF2(DynamicScaleFactors(*p["muF2"]))
    if "mu2" in p:
        o.setMu2(DynamicScaleFactors(*p["mu2"]))
    if "lambdaQCD" in p:
        o.setLambdaQCD(p["lambdaQCD"])
    if "alphaSByLHAPDF" in p:
        o.setAlphaSByLHAPDF(*p["alphaSByLHAPDF"])
    if "xBjorken" in p:
        o.setXBjorken(p["xBjorken"])
    if "hadronicS" in p:
        o.setHadronicS(p["hadronicS"])
    if "flags" in p:
        for k in p["flags"]:
            setattr(o.flags(), k, p["flags"][k])
    # leptonic setter
    if "alphaEM" in p:
        o.setAlphaEM(p["alphaEM"])
    if "polarizeBeams" in p:
        o.setPolarizeBeams(p["polarizeBeams"])
    if "leptonicS" in p:
        o.setLeptonicS(p["leptonicS"])
    if "Q2min" in p:
        o.setQ2min(p["Q2min"])
    if "Q2minByHVQDIS" in p:
        o.setQ2minByHVQDIS(p["Q2minByHVQDIS"])
    # int config has to be set AFTER flags!
    run = p["run"]
    if "IntegrationConfig" in p:
        method = run if isinstance(run, str) else run[0]
        for k in p["IntegrationConfig"]:
            setattr(o.getIntegrationConfig(method), k, p["IntegrationConfig"][k])


def global_run(o, p):
    """Run global functions"""
    run = p["run"]
    # find run
    if isinstance(run, str):
        p["res"] = getattr(o, run)()
    elif isinstance(run, tuple) or isinstance(run, list):
        p["res"] = getattr(o, run[0])(*run[1:])


def global_post_run(o, p):
    """Post-process run stuff"""
    p["IntegrationOutput"] = o.getIntegrationOutput()
    if "msg" in p:
        print(p["msg"])


# thread worker for Inclusive
def run_inclusive(p):
    """Inclusive working method"""
    # setup
    o = InclusiveLeptoProduction(*p["objArgs"])

    # setters
    global_setter(o, p)
    if "Delta" in p:
        o.setDelta(p["Delta"])

    # runners
    global_run(o, p)
    # post process
    global_post_run(o, p)
    return p


InclusiveRunner = Runner(run_inclusive)
"""Wrapper for InclusiveLeptoProduction to run several jobs in parallel"""


def run_fully_diff(p):
    """Fully differential working method"""
    # setup
    o = FullyDiffLeptoProduction(*p["objArgs"])

    # setters
    global_setter(o, p)
    if "xTilde" in p:
        o.setXTilde(p["xTilde"])
    if "omega" in p:
        o.setOmega(p["omega"])
    if "deltax" in p:
        o.setDeltax(p["deltax"])
    if "deltay" in p:
        o.setDeltay(p["deltay"])

    if "activateHistograms" in p:
        for k in p["activateHistograms"]:
            o.activateHistogram(*k)
    # run
    global_run(o, p)
    # post process
    global_post_run(o, p)
    return p


FullyDiffRunner = Runner(run_fully_diff)
"""Wrapper for FullyDiffLeptoProduction to run several jobs in parallel"""
