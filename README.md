# LeProHQ
(Un-)Polarized Lepto-Production of Heavy Quarks.

[![thesis](https://img.shields.io/badge/DOI-%20%20%20%2010.15496%2Fpublikation--34811-blue)](https://inspirehep.net/literature/1757437)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5220831.svg)](https://doi.org/10.5281/zenodo.5220831)

## Citation Policy

Please cite
```bibtex
@phdthesis{Hekhorn:2019nlf,
    author = "Hekhorn, Felix",
    title = "{Next-to-Leading Order QCD Corrections to Heavy-Flavour Production in Neutral Current DIS}",
    eprint = "1910.01536",
    archivePrefix = "arXiv",
    primaryClass = "hep-ph",
    doi = "10.15496/publikation-34811",
    school = "Tubingen U., Math. Inst.",
    year = "2019"
}
```

## Structure of this Repository

- [LeProHQ++](https://github.com/felixhekhorn/LeProHQ/tree/main/LeProHQ%2B%2B) contains the core C++ implementations that compute the actual structure functions.

- [LeProHQutils](https://github.com/felixhekhorn/LeProHQ/tree/main/LeProHQutils) contains set of helper functions.

- [LeProHQpy](https://github.com/felixhekhorn/LeProHQ/tree/main/LeProHQpy) contains the stand-alone Python wrapper for the fully inclusive coefficient functions as is required by [yadism](https://n3pdf.github.io/yadism/).

- [Patches](https://github.com/felixhekhorn/LeProHQ/tree/main/Patches) contains a list of patches that are either required to run the C++ implementation or to run some legacy PDF codes.

- [Benchmarks](https://github.com/felixhekhorn/LeProHQ/tree/main/Benchmarks) contains a list of benchmark results to verify the implementations.
