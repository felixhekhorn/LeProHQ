# NPB540-345

The photo-production limit of heavy quark production was computed [here](https://inspirehep.net/literature/) and the authors also [published a subsequent study](https://inspirehep.net/literature/) in which they investigated the inclusive distributions: i.e. the rapidity distribution and the Feynman momentum distribution of the heavy anti-quark.

## Running

- You need to have the `GRSV` and `GRV` pdfs compiled (as used in the original paper)
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
TODO
```