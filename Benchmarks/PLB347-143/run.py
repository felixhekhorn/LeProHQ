"""Compare to PLB347-143."""

import numpy as np
import matplotlib.pyplot as plt

import LeProHQ
import hqcoef

# fs = (("cg0","F2"),"c2log")
# fs = (("dq1","F2"),"d2nloq")
# fs = (("cg1","F2"),"c2nlog")
# fs = (("cgBar1","F2"),"c2nlobarg")
# fs = (("cq1", "F2"), "c2nloq")
fs = (("cq1", "FL"), "clnloq")
# fs = (("cqBarF1","F2"),"c2nlobarq")
xis = np.geomspace(1e-3, 1e3, 30)
etas = np.geomspace(1e-3, 1e3, 30)
# threshold = 1e-2

# collect data
me = []
other = []
for xi in xis:
    for eta in etas:
        me.append(getattr(LeProHQ, fs[0][0])(fs[0][1], "VV", xi, eta))
        other.append(getattr(hqcoef, fs[1])(eta, xi))
me = np.array(me).reshape(len(xis), len(etas))
other = np.array(other).reshape(len(xis), len(etas))

# plot
fig, axs = plt.subplots(2, 2)
for data, lab, axs in zip(
    [me / other - 1.0, me - other], ["log10(rel. Error)", "log10(abs. Error)"], axs
):
    cax = axs[0].imshow(
        np.log10(np.abs(data)),
        extent=(
            np.log10(etas).min(),
            np.log10(etas).max(),
            np.log10(xis).min(),
            np.log10(xis).max(),
        ),
        origin="lower",
    )
    axs[0].set_xlabel(r"$\log_{10}(\eta)$")
    axs[0].set_ylabel(r"$\log_{10}(\xi)$")
    cbar = fig.colorbar(cax)
    cbar.set_label(lab)
    axs[1].hist(np.log10(np.abs(data)).flatten(), 30)
    axs[1].set_xlabel(lab)
    axs[1].set_ylabel("Pixel")


fig.tight_layout()
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
