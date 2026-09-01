# Staging release notes

Version 0.1.0 is a pre-publication staging release prepared from the MATLAB
R2025a verification pass of 2026-09-01.

It adds deterministic drivers and reference CSV files for every numerical
table, makes figure output work in a standalone checkout, fixes previous-
iterate bookkeeping in the multiplier-free Reynolds solver, aligns smoothed
cavity diagnostics with the central-path threshold, and documents the
sampling conventions for discontinuous pressures.

The final certification pass adds executable checks for quantitative claims
in the manuscript prose and a complete safeguarded 48-case Newton matrix. It
also records the corrected Newton residual-order description used in the
submission manuscript.

Before public release:

- replace provisional citation metadata with the final paper DOI or preprint
- create the GitHub repository and set its URL in `CITATION.cff`

The source is licensed under the BSD 3-Clause License.
The quick reproduction test passed from the final local repository location
after the license and metadata update.
