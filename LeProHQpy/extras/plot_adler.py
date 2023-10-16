import numpy as np
from scipy.integrate import quad
from LeProHQ import dq1, Adler
from matplotlib import pyplot as plt


ETA_MAX = 1e9

def dq(z, proj, cc, xi):
    if xi * (1. - z) / z <= 4 :
        return 0.0
    eta = xi / 4.0 * (1.0 / z - 1.0) - 1.0
    eta = min(eta, ETA_MAX)
    r = (
        xi / np.pi
        / z
        * (4.0 * np.pi) ** 2
        * dq1(proj,cc, xi, eta)
    )
    return -r

xis = np.geomspace(1e-4,1e8,30)
proj,cc = "F2", "VV"

from_quad = []
from_mma = []
for xi in xis:
    from_quad.append(quad(dq,0.,1.,args=(proj, cc, xi)))
    from_mma.append(Adler(proj, cc, xi))

from_quad = np.array(from_quad)
from_mma = np.array(from_mma)

# print(xis)
# print(from_mma)
print(from_quad)

fig = plt.figure()
fig.suptitle(f"{proj=},{cc=},{ETA_MAX=}")
ax0 = fig.add_subplot(3, 1, (1, 2))
ax0.semilogx(xis, from_quad[:,0],label="quad")
ax0.semilogx(xis, from_mma,label="MMa")
# plt.setp(ax0.get_xticklabels(), visible=False)
ax0.legend()
ax1 = fig.add_subplot(3, 1, 3, sharex=ax0)
ax1.semilogx(xis, from_quad[:,0] / from_mma, label="quad/MMa")
ax1.legend()
fig.savefig(f"adler.pdf")
plt.close(fig)
