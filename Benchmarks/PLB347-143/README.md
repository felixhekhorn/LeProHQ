# PLB347-143

The authors provide a fast (tabulated) implementation of the fully inclusive coefficient functions [here](https://inspirehep.net/literature/380446) that is used e.g. in [OPENQCDRAD](https://www-zeuthen.desy.de/~alekhin/OPENQCDRAD/), [APFEL](https://apfel.mi.infn.it/) or [xFitter](https://www.xfitter.org/xFitter/).

## Running

- compile the Python interface to `hqcoeff`:
  - get the `hqcoeff` code from inside xfitter: https://gitlab.cern.ch/fitters/xfitter/-/blob/master/src/ABM/hqcoef.f
  - apply the provided patch: `hqcoeff.patch`
  - run the  `f2py.sh` script to compile (the provided `hqcoef.pyf` should extract all necessary bindings from Fortran)
- adjust `run.py` to your needs and run it. It prints a direct comparison, showing the differences


## References

```bibtex
@article{Riemersma:1994hv,
    author = "Riemersma, S. and Smith, J. and van Neerven, W. L.",
    title = "{Rates for inclusive deep inelastic electroproduction of charm quarks at HERA}",
    eprint = "hep-ph/9411431",
    archivePrefix = "arXiv",
    reportNumber = "SMU-HEP-94-25, ITP-SB-94-59, INLO-PUB-16-94",
    doi = "10.1016/0370-2693(95)00036-K",
    journal = "Phys. Lett. B",
    volume = "347",
    pages = "143--151",
    year = "1995"
}
```
