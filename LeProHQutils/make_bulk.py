"""Generate bulk data files."""

from runner import InclusiveRunner

# global parameters
Delta = 1e-6

nlf = 3
m2 = 1.0
intCfg = {"verbosity": 1}

q2 = 1.0
partonic_s = 8.0*1.0
e = {
    "objArgs": (nlf, m2, Delta),
    "Q2": q2,
    "partonicS": partonic_s,
    "run": "cq1_VV"
}
InclusiveRunner.append(e)
print(InclusiveRunner.run())

xi_min = 0.01
xi_max = 2.5e3
eta_min = 0.09
eta_max = 1e6
