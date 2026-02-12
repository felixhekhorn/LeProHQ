"""Compare to PLB347-143."""

import numpy as np
import matplotlib.pyplot as plt

import LeProHQ
import hqcoef

# fs = (("cg0","F2"),"c2log")
# fs = (("dq1","F2"),"d2nloq")
# fs = (("cg1","F2"),"c2nlog")
# fs = (("cgBar1","F2"),"c2nlobarg")
fs = (("cq1", "F2"), "c2nloq")
# fs = (("cqBarF1","F2"),"c2nlobarq")
xis = np.geomspace(1e-3, 1e3, 30)
etas = np.geomspace(1e-3, 1e3, 30)
threshold = 1e-2

# collect data
me = []
other = []
for xi in xis:
    for eta in etas:
        me.append(LeProHQ.__getattribute__(fs[0][0])(fs[0][1], "VV", xi, eta))
        other.append(hqcoef.__getattribute__(fs[1])(eta, xi))
me = np.array(me).reshape(len(xis), len(etas))
other = np.array(other).reshape(len(xis), len(etas))

# plot
fig, ax = plt.subplots(1, 1)
cax = ax.imshow(
    np.log(np.abs(me / other - 1.0)),
    extent=(
        np.log(etas).min(),
        np.log(etas).max(),
        np.log(xis).min(),
        np.log(xis).max(),
    ),
    origin="lower",
)
ax.set_xlabel(r"$\log(\eta)$")
ax.set_ylabel(r"$\log(\xi)$")
cbar = fig.colorbar(cax)
cbar.set_label("rel. Error")
fig.savefig(f"{fs[0][0]}-{fs[0][1]}.pdf")
plt.close(fig)
# print(me)
#         rel_err = np.abs(me / other - 1.0)
#         if rel_err < threshold:
#             continue
#         print(f"{xi:.2e}\t{eta:.2e}\t{me:.2e}\t{other:.2e}\t{rel_err}")
#         c += 1
# tot = len(xis) * len(etas)
# print(f"{c}/{tot} points of {fs[0][0]}@{fs[0][1]} exceed the threshold of {threshold}")
