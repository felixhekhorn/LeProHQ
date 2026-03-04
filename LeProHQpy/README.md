# LeProHQpy

(Un-)Polarized Lepto-Production of Heavy Quarks.

This is the stand-alone Python wrapper for the fully-inclusive coefficient functions.

To see this implementation of the coefficient functions in action, i.e. actual structure functions, please use [yadism](https://github.com/NNPDF/yadism).

Here is a snippet of the actual usage:

```python
import LeProHQ

proj = "F2"  # F2 structure function
cc = "VV"  # vectorial-vectorial coupling
m2 = 2.25  # [GeV^2] heavy quark mass
Q2 = 2.25  # [GeV^2] virtuality
xi = Q2 / m2  # []
eta = 1.0  # [] (partonic) distance from threshold (in units of threshold)

LeProHQ.cq1(proj, cc, xi, eta)
```
