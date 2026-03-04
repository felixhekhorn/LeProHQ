"""Generate bulk data files."""

import numpy as np
from runner import InclusiveRunner

# global parameters
nlf = 3
m2 = 1.0
Delta = 1e-6
intCfg = {"verbosity": 1}
fnc = "cq1"
proj = "FL"
cc = "VV"
path = f"{fnc}-{proj}_{cc}-bulk.dat"

# grid parameters
lneta_min = np.log(9e-1)
lneta_max = np.log(1e5)
n_eta = 51
lnetas = np.linspace(lneta_min, lneta_max, n_eta)
lnxi_min = np.log(1e-2)
lnxi_max = np.log(1e5)
n_xi = 41
lnxis = np.linspace(lnxi_min, lnxi_max, n_xi)

# collect grid
for lneta in lnetas:
    for lnxi in lnxis:
        q2 = np.exp(lnxi) * m2
        e = {
            "objArgs": (nlf, m2, Delta),
            "Q2": q2,
            "partonicEta": np.exp(lneta),
            "projection": proj,
            "run": f"{fnc}_{cc}",
            "IntegrationConfig": intCfg,
        }
        InclusiveRunner.append(e)
# Run!
grid = InclusiveRunner.run(5)
# resort
res_grid = np.array(list(map(lambda e: e["res"], grid))).reshape(n_eta, n_xi)
# insert borders
bulk_grid = np.insert(np.insert(res_grid, 0, lnetas, axis=1), 0, [0.0, *lnxis], axis=0)
# save
np.savetxt(path, bulk_grid)
