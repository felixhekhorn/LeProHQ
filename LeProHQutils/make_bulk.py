"""Generate bulk data files."""

import numpy as np
from runner import InclusiveRunner

# global parameters
Delta = 1e-6
nlf = 3
m2 = 1.0
intCfg = {"verbosity": 1}
fnc = "cq1"
proj = "F2"
cc = "VV"
path = f"{fnc}-{proj}_{cc}-bulk.dat"

# grid parameters
lnxi_min = np.log(0.01)
lnxi_max = np.log(2.5e3)
n_xi = 2  # 41
lnxis = np.linspace(lnxi_min, lnxi_max, n_xi)
lneta_min = np.log(0.09)
lneta_max = np.log(1e6)
n_eta = 3  # 51
lnetas = np.linspace(lneta_min, lneta_max, n_eta)

# collect grid
for lnxi in lnxis:
    for lneta in lnetas:
        q2 = np.exp(lnxi) * m2
        partonic_s = 4.0 * m2 + 4.0 * m2 * np.exp(lneta)
        e = {
            "objArgs": (nlf, m2, Delta),
            "Q2": q2,
            "partonicS": partonic_s,
            "projection": proj,
            "run": f"{fnc}_{cc}",
            "IntegrationConfig": intCfg,
        }
        InclusiveRunner.append(e)
# Run!
grid = InclusiveRunner.run()
# resort
res_grid = np.array(list(map(lambda e: e["res"], grid))).reshape(n_xi, n_eta)
# insert borders
bulk_grid = np.insert(np.insert(res_grid, 0, lnxis, axis=1), 0, [0.0, *lnetas], axis=0)
# save
np.savetxt(path, bulk_grid)
