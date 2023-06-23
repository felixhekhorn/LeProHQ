# LeProHQ++


## Build

We use the [meson build system](https://mesonbuild.com/). Make sure to have all external dependencies listed below available on your system.

To compile the code run

```sh
meson setup builddir
meson compile -C builddir
```

### External Dependencies
- [Python](https://www.python.org/) >= 3.9
- [gsl](https://www.gnu.org/software/gsl/) >= 2.7
- [Relativistic Kinematics (RK)](https://rk.hepforge.org/) >= 1.7
- [Dvegas](https://dvegas.hepforge.org/) >= 2.0.3 with patched [`dvegas.h`](https://github.com/felixhekhorn/LeProHQ/blob/main/Patches/dvegas.h.patch) and [`dvegas.cpp`](https://github.com/felixhekhorn/LeProHQ/blob/main/Patches/dvegas.cpp.patch)
- [boost](https://www.boost.org/) >= 1.74 - used modules are `filesystem` and the headers for `format.hpp`
- [lhapdf](https://lhapdf.hepforge.org/) >= 6.4

## Run

The build process generates a Python library `LeProHQpp` at `builddir/src/LeProHQpp.so`.
(Note, this is the Python binding of the full fledged program and not only the fully inclusive coefficient functions,
which instead are provided by [LeProHQpy](https://github.com/felixhekhorn/LeProHQ/tree/main/LeProHQpy)
through the Python library `LeProHQ`.)

Here, a short snippet which demonstrates how to use this library:
```py
import LeProHQpp
from LeProHQpp import Projection as proj

# input configuration
nlf = 3
m2 = 1.51**2
Delta = 1e-6

# compute 2*x*g_1(x=0.1, Q2=10) at LO
obj = LeProHQpp.InclusiveLeptoProduction(nlf, m2, Delta)
obj.setProjection(proj.x2g1)
obj.setPdf("DSSV_REP_LHAPDF6", 0)
obj.setAlphaSByLHAPDF("DSSV_REP_LHAPDF6", 0)
obj.setQ2(10.0)
obj.setXBjorken(0.1)
obj.flags().useNextToLeadingOrder = False
print(obj.F())
```

[LeProHQutils](https://github.com/felixhekhorn/LeProHQ/tree/main/LeProHQutils) contains a convenience wrapper to this core C++ implementations
providing a mapping between dictionaries and the respective calls. This is handy when running very often with different configurations - see
the [Benchmarks](https://github.com/felixhekhorn/LeProHQ/tree/main/Benchmarks) for examples.
