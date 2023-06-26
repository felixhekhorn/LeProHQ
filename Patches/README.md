# Patches

We provide several patch files to improve and fix third party codes.

If the patch fails take a look to the actual changes and do them yourself (since this refers to third party code).

## DSSV

The original code is private and can be retrieved from the [authors](https://inspirehep.net/literature/1291054).
`DSSV_gluon_update.patch` changes the access to the underlying grid files to accept direct file path.
This is needed for the [`PdfWrapper`](https://github.com/felixhekhorn/LeProHQ/blob/main/LeProHQ%2B%2B/src/Pdf/PdfWrapper.cpp).

Note, that later a [MC replica version was published](https://inspirehep.net/literature/1722265).

## GRSV

The original code is private and can be retrieved from the [authors](https://inspirehep.net/literature/398625).
`grsv.patch` changes the access to the underlying grid files to accept direct file path.
This is needed for the [`PdfWrapper`](https://github.com/felixhekhorn/LeProHQ/blob/main/LeProHQ%2B%2B/src/Pdf/PdfWrapper.cpp).

## Dvegas

The original code can be obtained from [here](https://dvegas.hepforge.org/).
`dvegas.h.patch` and `dvegas.cpp.patch` fix some exception declarations and turns the empty integration into an exception (rather than a failure).
