import sys
import pathlib

import numpy as np

here = pathlib.Path(__file__)
sys.path.append(here.parents[2] / "LeProHQpy")

import LeProHQ
import hqcoef

#fs = (("cg0","F2"),"c2log")
#fs = (("dq1","F2"),"d2nloq")
#fs = (("cg1","F2"),"c2nlog")
fs = (("cgBar1","F2"),"c2nlobarg")
#fs = (("cq1","F2"),"c2nloq")
#fs = (("cqBarF1","F2"),"c2nlobarq")
xis = np.geomspace(1e-3,1e3,30)
etas = np.geomspace(1e-3,1e3,30)
threshold = 1e-2

print(f"xi\t\teta\t\tLeProHQ\t\thqcoeff\t\trel. err")
c = 0
for xi in xis:
    for eta in etas:
        me = LeProHQ.__getattribute__(fs[0][0])(fs[0][1],"VV",xi, eta)
        other = hqcoef.__getattribute__(fs[1])(eta,xi)
        rel_err = np.abs(me/other - 1.)
        if rel_err < threshold:
            continue
        print(f"{xi:.2e}\t{eta:.2e}\t{me:.2e}\t{other:.2e}\t{rel_err}")
        c += 1
tot = len(xis) * len(etas)
print(f"{c}/{tot} points of {fs[0][0]}@{fs[0][1]} exceed the threshold of {threshold}")
