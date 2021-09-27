# -*- coding: utf-8 -*-

import re

from joblib import Parallel, delayed

from LeProHQpp import (DynamicScaleFactors, FullyDiffLeptoProduction,
                       InclusiveLeptoProduction)


# provide runner
class AbstractRunner:
    """abstract wrapper to run several jobs in parallel"""

    def __init__(self, callee):
        """constructor"""
        self.l = []
        self.callee = callee

    def append(self, e):
        """
        append an element
        @param e element
        """
        self.l.append(e)

    def run(self, n_jobs=None):
        """
        computes all elements
        @param n_jobs number of parallel threads
        """
        print("Computing %d points ..." % len(self.l))
        return Parallel(n_jobs=n_jobs, backend="threading")(
            delayed(self.callee)(e) for e in self.l
        )


def global_setter(o, p):
    """set global vars"""
    if "run" not in p:
        raise KeyError("no called function set")
    run = p["run"]
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
    if "IntegrationConfig" in p:
        method = run if type("") == type(run) else run[0]
        for k in p["IntegrationConfig"]:
            setattr(o.getIntegrationConfig(method), k, p["IntegrationConfig"][k])


def global_run(o, p):
    """run global functions"""
    run = p["run"]
    # find run
    if type("") == type(run):
        if "F" == run:
            p["res"] = o.F()
        elif "sigma" == run:
            p["res"] = o.sigma()
        else:
            if None == re.match("^[cd][gq](Bar[RF]?)?[01]_[VA][VA]$", run):
                raise ValueError("unknown function: " + run)
            p["res"] = getattr(o, run)()


def global_post_run(o, p):
    """post run stuff"""
    if "IntegrationOutput" in p:
        p["IntegrationOutput"] = o.getIntegrationOutput()
    if "msg" in p:
        print(p["msg"])


# thread worker for Inclusive
def run_inclusive(p):
    """inclusive working method"""
    # setup
    o = InclusiveLeptoProduction(*p["objArgs"])

    # setters
    global_setter(o, p)
    if "Delta" in p:
        o.setDelta(p["Delta"])

    # runners
    global_run(o, p)
    run = p["run"]
    if tuple == type(run) or list == type(run):
        if "dF_dHAQTransverseMomentum" == run[0]:
            p["res"] = o.dF_dHAQTransverseMomentum(run[1])
        elif "dF_dHAQTransverseMomentumScaling" == run[0]:
            p["res"] = o.dF_dHAQTransverseMomentumScaling(run[1])
        elif "dF_dHAQRapidity" == run[0]:
            p["res"] = o.dF_dHAQRapidity(run[1])
        elif "dF_dHAQFeynmanX" == run[0]:
            p["res"] = o.dF_dHAQFeynmanX(run[1])
        else:
            # Util.pWarn(p)
            raise ValueError("unknown function: " + run)
    # post process
    global_post_run(o, p)
    return p


class InclusiveRunner(AbstractRunner):
    """wrapper for InclusiveLeptoProduction to run several jobs in parallel"""

    def __init__(self):
        """constructor"""
        super().__init__(run_inclusive)


def run_fully_diff(p):
    """fully differential working method"""
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


class FullyDiffRunner(AbstractRunner):
    """wrapper for FullyDiffLeptoProduction to run several jobs in parallel"""

    def __init__(self):
        """constructor"""
        super().__init__(run_fully_diff)
