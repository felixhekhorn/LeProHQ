import numpy as np
import matplotlib.pyplot as plt
from LeProHQ import cg1

m = 1.5
Q2 = 10
cc = "VV"
xi = Q2 / (m**2)

projections = ["x2g1"]
titles = {
    "F2": "cg1 vs Eta for proj=F2",
    "FL": "cg1 vs Eta for proj=FL",
    "x2g1": "cg1 vs Eta for proj=x2g1",
}
eta_values = np.linspace(0.085, 0.105, 100)  # np.logspace(-2, -0.5, 100)
for proj in projections:
    cg1_values = [cg1(proj, cc, xi, eta) for eta in eta_values]
    cg1_values = np.nan_to_num(cg1_values)
    plt.figure(figsize=(8, 6))
    plt.plot(eta_values, cg1_values, "o", label=f"cg1 vs eta for proj={proj}")
    plt.xscale("log")
    plt.xlabel("Eta (log scale)")
    plt.ylabel("cg1")
    plt.title(f"{titles[proj]} (log scale) for cc={cc}, m={m}, Q^2={Q2}")
    plt.grid(True, which="both", ls="--")
    plt.legend()
plt.savefig("test.pdf")
