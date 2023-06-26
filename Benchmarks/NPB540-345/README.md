# NPB540-345

The photo-production limit of heavy quark production was computed [here](https://inspirehep.net/literature/473492) in which they investigate the inclusive distributions: i.e. the rapidity distribution and the Feynman momentum distribution of the heavy anti-quark.

## Running

- You need to have the `GRSV` and `GRV` pdfs compiled (as used in the original paper)
- to use the inclusive implementation:
  - edit (the bottom of) `run-Inclusive.py`
  - run the file
  - eventual intermediate results will be stored as a single YAML file in the `Inclusive/` subdirectory
- to use the fully differential implementation
  - edit (the bottom of) `run-FullyDiff.py`
  - run the file
  - the histograms will be stored as `.dat` files in the `FullyDiff/` subdirectory
- each implementation will create a figure as `.pdf` that should be similar to the original in the respective subdirectory

## References

```bibtex
@article{Bojak:1998zm,
    author = "Bojak, I. and Stratmann, M.",
    title = "{Photoproduction of heavy quarks in next-to-leading order QCD with longitudinally polarized initial states}",
    eprint = "hep-ph/9807405",
    archivePrefix = "arXiv",
    reportNumber = "DO-TH-98-12, DTP-98-36",
    doi = "10.1016/S0550-3213(98)00751-2",
    journal = "Nucl. Phys. B",
    volume = "540",
    pages = "345--381",
    year = "1999",
    note = "[Erratum: Nucl.Phys.B 569, 694--694 (2000)]"
}
```
