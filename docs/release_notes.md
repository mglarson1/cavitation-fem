# Review-stage release notes

## 0.2.0 -- September 8, 2026

This snapshot accompanies Tables 1--14 of the revised manuscript. It replaces
the old Table 2 averaged strip errors with raw two-dimensional Reynolds
errors, adds checked unsmoothed QP references and smoothing comparisons,
and adds the original and pressure-compatible end-condition studies.
The new Stokes solver options preserve the original default formulation.

The new comparisons remove the original large cavity-length ratio after
end-condition calibration, while retaining a smaller front-position
difference on the refined reference domain. Grid and domain sensitivity
remain explicit limitations. Neither smoothing convergence nor pressure
agreement certifies a free-boundary error bound.

The current code requires Optimization Toolbox for its QP studies. Reference
summary CSVs cover all manuscript tables and supporting diagnostics. Dense
profiles and MAT fields are regenerated, rather than included in Git.
No paper DOI or public archival DOI has been assigned by this update.

## 0.1.0 -- September 1, 2026

The earlier snapshot supplied Tables 1--10, all seven figure drivers,
prose-claim checks and the original/safeguarded 48-case Newton matrices.
It is preserved in Git history. The code remains under BSD-3-Clause.
