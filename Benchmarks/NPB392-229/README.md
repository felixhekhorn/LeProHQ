# NPB392-229

In parallel with the original NLO heavy quark contributions [published here](https://inspirehep.net/literature/335018), the authors also [published a first study](https://inspirehep.net/literature/339361) in which they investigated the inclusive distributions: i.e. the rapidity distribution and the transverse momentum distribution of the heavy anti-quark.

## Running

- You need to have the `MorfinTungB` pdf compiled (as used in the original paper)
- to use the inclusive implementation:
  - edit (the bottom of) `run-Inclusive.py`
  - run the file
  - eventual intermediate results will be stored as a single YAML file in the `Inclusive/` subdirectory
- to use the fully differential implementation
  - edit (the bottom of) `run-FullyDiff.py`
  - run the file
  - the histograms will be stored as `dat` files in the `FullyDiff/` subdirectory
- each implementation will create a figure as `pdf` that should be similar to the original in the respective subdirectory

## References

```
@article{Laenen:1992zk,
    author = "Laenen, Eric and Riemersma, S. and Smith, J. and van Neerven, W. L.",
    title = "{Complete O (alpha-s) corrections to heavy flavor structure functions in electroproduction}",
    reportNumber = "ITP-SB-92-09",
    doi = "10.1016/0550-3213(93)90201-Y",
    journal = "Nucl. Phys. B",
    volume = "392",
    pages = "162--228",
    year = "1993"
}
@article{Laenen:1992xs,
    author = "Laenen, Eric and Riemersma, S. and Smith, J. and van Neerven, W. L.",
    title = "{O(alpha-s) corrections to heavy flavor inclusive distributions in electroproduction}",
    reportNumber = "FERMILAB-PUB-92-271-T, ITP-SB-92-41",
    doi = "10.1016/0550-3213(93)90202-Z",
    journal = "Nucl. Phys. B",
    volume = "392",
    pages = "229--250",
    year = "1993"
}
```