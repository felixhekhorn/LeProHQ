"""Mellin moments for 2602.04455.

Get the necessary python packages via `pip install numpy scipy pandas LeProHQ`.
"""

from scipy.integrate import quad
import LeProHQ
import numpy as np
import numpy.typing as npt
import pandas as pd

# see caption of table 2
m2 = 8.0  # GeV^2
Q2 = 49.0  # GeV^2
xi = Q2 / m2  # Eq. 2
ns = np.arange(2.0, 22.0 + 1.0, 2.0)


def ker(z: float, n: float, fn, proj: str, xi: float) -> float:
    """Mellin kernel."""
    eta = xi / 4.0 * (1.0 / z - 1.0) - 1.0  # Eq. 2
    # The paper focuses only on the vector-vector coupling ~ VV
    return z ** (n - 2.0) * fn(proj, "VV", xi, eta)


def tabulate(fnc, norm) -> npt.NDArray[np.float64]:
    """Generate a benchmark table"""
    data = [ns]
    # loop projections
    for proj in ["FL", "F2"]:
        # loop n
        data.append(
            norm
            * np.array(
                [
                    quad(ker, 0.0, 1.0 / (1.0 + 4.0 / xi), args=(n, fnc, proj, xi))[0]
                    for n in ns
                ]
            )
        )
    return np.array(data).T


def cg1(proj: str, cc: str, xi: float, eta: float) -> float:
    """Set mu2 = Q2"""
    return LeProHQ.cg1(proj, cc, xi, eta) + LeProHQ.cgBar1(proj, cc, xi, eta) * np.log(
        xi
    )


def cq1(proj: str, cc: str, xi: float, eta: float) -> float:
    """Set mu2 = Q2"""
    return LeProHQ.cq1(proj, cc, xi, eta) + LeProHQ.cqBarF1(proj, cc, xi, eta) * np.log(
        xi
    )


# there are some global normalizations flowing around
norm_lo = xi / np.pi
norm_nlo = xi * 4.0**2 * np.pi

# Collect tables
tabs = {
    "2": tabulate(cq1, norm_nlo),
    "3": tabulate(LeProHQ.dq1, norm_nlo),
    "4": tabulate(LeProHQ.cg0, norm_lo),
    "5": tabulate(cg1, norm_nlo),
}
# Save them
for lab, tab in tabs.items():
    np.savetxt(f"tab{lab}.txt", tab)

# Compare LeProHQ vs KMS vs RSvN
for lab, tab in tabs.items():
    ref_tab = np.loadtxt(f"tab{lab}-ref.txt")
    for j, fk in enumerate(["FL", "F2"]):
        df = pd.DataFrame()
        df["N"] = ns
        df["LeProHQ"] = tab[:, 1 + j]
        df["KMS"] = ref_tab[:, 1 + 2 * j]
        df["RSvN"] = ref_tab[:, 2 + 2 * j]
        df["LeProHQ/KMS"] = df["LeProHQ"] / df["KMS"]
        df["LeProHQ/RSvN"] = df["LeProHQ"] / df["RSvN"]
        df["KMS/RSvN"] = df["KMS"] / df["RSvN"]
        print(f"Checking {fk} in Table {lab}")
        print(df)
        print()
