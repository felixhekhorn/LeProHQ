import numpy as np

from LeProHQpp import Projection as proj
import runner

# global parameters
nlf = 3
m2 = 1.5**2
Delta = 1e-6
Q2 = 10
lambdaQCD = 0.194
pdf = ("MorfinTungB",0)
intCfg = {"verbosity": 1}
default_flags = {"usePhoton": True,"usePhotonZ": False,"useZ": False}
flags = { # split by different contributions
    "LO": {"useLeadingOrder": True , "useNextToLeadingOrder": False}, # LO result
    #"NLO": {"useLeadingOrder": True , "useNextToLeadingOrder": True}, # NLO result
}


projs = [proj.F2, proj.FL]

# rapidity configuration
rap_bins = np.linspace(-1., 1., 11)
rap_setups = (
    (.1,2.),(.01,3.5),(.001,4.5),(.0001, 5.5)
)

#	mu2 = (4.,-1.,0.,4.,)

objArgs = (nlf, m2, Delta)

r = runner.InclusiveRunner()
def append_rapidity(r, proj_, xBj, y0, flag_label, tag):
    for bin_, y in enumerate(y0 * rap_bins):
        cur_flags = default_flags.copy()
        cur_flags.update(flags[flag_label])
        r.append({
            "objArgs": objArgs, "projection": proj_, "Q2": Q2, "xBjorken": xBj, "lambdaQCD": lambdaQCD, "pdf": pdf, "IntegrationConfig": intCfg,
            "run": ("dF_dHAQRapidity", y), "flags": cur_flags, "bin": bin_, "tag": tag, "IntegrationOutput": True
            })


def append_all_rapidities(r):
    for j, p in enumerate(projs):
        for k, (xBj, y0) in enumerate(rap_setups):
            for l, flag_label in enumerate(flags.keys()):
                append_rapidity(r, p, xBj, y0, flag_label, (j,k,l))

append_all_rapidities(r)

l = r.run(-2)

def reorder(l):
    m = np.full((len(projs),len(rap_setups),len(flags),len(rap_bins)), np.nan)
    for e in l:
        m[(*e["tag"],e["bin"])] = e["res"]
    return m

print(reorder(l))
