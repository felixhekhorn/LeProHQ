import numpy as np
import matplotlib.pyplot as plt

import yaml

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
    "NLO": {"useLeadingOrder": True , "useNextToLeadingOrder": True}, # NLO result
}


projs = [proj.F2, proj.FL]

# rapidity configuration
n_rap_bins = 5

#	mu2 = (4.,-1.,0.,4.,)

objArgs = (nlf, m2, Delta)

class RapidityRunner:
    setups = (
        (.1,2.),(.01,3.5),(.001,4.5),(.0001, 5.5)
    )

    def __init__(self, n_bins):
        self.r = runner.InclusiveRunner()
        self.m = None
        self.bins = np.linspace(-1., 1., n_bins)

    def append(self,proj_, xBj, y0, flag_label, tag):
        for bin_, y in enumerate(y0 * self.bins):
            cur_flags = default_flags.copy()
            cur_flags.update(flags[flag_label])
            self.r.append({
                "objArgs": objArgs, "projection": proj_, "Q2": Q2, "xBjorken": xBj, "lambdaQCD": lambdaQCD, "pdf": pdf, "IntegrationConfig": intCfg,
                "run": ("dF_dHAQRapidity", y), "flags": cur_flags, "bin": bin_, "tag": tag, "IntegrationOutput": True
                })

    def append_all(self):
        for j, p in enumerate(projs):
            for k, (xBj, y0) in enumerate(self.setups):
                for l, flag_label in enumerate(flags.keys()):
                    self.append(p, xBj, y0, flag_label, (j,k,l))

    def run(self, n_jobs = None):
        self.append_all()
        l = self.r.run(n_jobs)
        self.m = self.reorder(l)

    def reorder(self, l):
        m_raw = np.full((len(projs),len(self.setups),len(flags),len(self.bins)), np.nan)
        m = {
            "res": m_raw.copy(),
            "error": m_raw.copy()
        }
        for e in l:
            m["res"][(*e["tag"],e["bin"])] = e["res"]
            m["error"][(*e["tag"],e["bin"])] = e["IntegrationOutput"].error
        return m

    def dump(self, fn):
        out = {
            "sorting": ["projections", "setups", "flags", "bins"],
            "projections": [str(p) for p in projs],
            "setups": [list(s) for s in self.setups],
            "flags": [f for f in flags],
            "bins": self.bins.tolist(),
            "res": self.m["res"].tobytes(),
            "error": self.m["error"].tobytes()
        }
        with open(fn,"w") as o:
            yaml.safe_dump(out,o)

    def load(self, fn):
        with open(fn, "r") as o:
            mm = yaml.safe_load(o)
        self.m = {}
        shape = (len(mm["projections"]),len(mm["setups"]),len(mm["flags"]),len(mm["bins"]))
        self.m["res"] = np.frombuffer(mm["res"]).reshape(shape)
        self.m["error"] = np.frombuffer(mm["error"]).reshape(shape)

    def plot(self, fp):
        fn = fp % 11
        fig = plt.figure()
        fig.suptitle(f"x = {self.setups[0][0]:f}")
        ax = fig.add_subplot(111)
        ax.set_xlabel("y (rapidity)")
        ax.set_ylabel("$dF_2(x,Q^2,m_c^2,y)/dy$")
        ax.semilogy(-self.bins * self.setups[0][1], self.m["res"][0,0,0], label="LO")
        ax.semilogy(-self.bins * self.setups[0][1], self.m["res"][0,0,1], label="NLO")
        ax.set_ylim([1e-6,1e-2])
        ax.set_xlim([-2.,2.])
        fig.legend()
        fig.savefig(fn)
        plt.show()
        plt.close(fig)

rr = RapidityRunner(n_rap_bins)
#rr.run(-2)
#rr.dump("rap.yaml")
rr.load("rap.yaml")
rr.plot("fig%d.pdf")
