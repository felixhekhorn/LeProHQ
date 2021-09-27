import numpy as np

from LeProHQpp import Projection as proj

nlf = 3
m2 = 1.5 * 1.5
Q2 = 1e-4  # fake photo-production
lambdaQCD = 0.2

default_flags = {"usePhoton": True, "usePhotonZ": False, "useZ": False}
flags = {  # split by different contributions
    "LO": {"useLeadingOrder": True, "useNextToLeadingOrder": False},  # LO result
    # "NLO": {"useLeadingOrder": True , "useNextToLeadingOrder": True}, # NLO result
}
intCfg = {"verbosity": 1}

setup_fig10 = []
for orderFlag, pdf in [("LO", "GRSV96STDLO"), ("NLO", "GRSV96STDNLO")]:
    for cMuF2 in [0.25, 1.0, 4.0]:
        for cMuR2 in [0.25, 1.0, 4.0]:
            setup_fig10.append(
                {
                    "hadronicS": 10.0 ** 2,
                    "orderFlag": orderFlag,
                    "pdf": pdf,
                    "cMuF2": cMuF2,
                    "cMuR2": cMuR2,
                }
            )
setup_fig10.append(
    {"hadronicS": 200.0 ** 2, "orderFlag": "LO", "pdf": pdf, "cMuF2": 1.0, "cMuR2": 1.0}
)
setup_fig10.append(
    {
        "hadronicS": 200.0 ** 2,
        "orderFlag": "NLO",
        "pdf": pdf,
        "cMuF2": 1.0,
        "cMuR2": 1.0,
    }
)


aem = 1.0 / 137.0
norm_sigma = (4.0 * np.pi * np.pi * aem) / Q2
nbTimesGeV2 = 3894.0 / 10000.0 * 1e6

setup_fig11_fig12 = [
    {"proj_": proj.x2g1, "pdf": "GRSV96STDLO", "orderFlag": "LO"},
    {"proj_": proj.x2g1, "pdf": "GRSV96STDNLO", "orderFlag": "NLO"},
    {"proj_": proj.F2, "pdf": "GRV94NLO", "orderFlag": "NLO"},
]
