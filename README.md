# Structure-preserving finite element methods for cavitation

MATLAB implementation and reproducibility package for the manuscript
*Structure-Preserving Augmented Lagrangian Finite Element Methods for
Cavitation in Reynolds and Stokes Flows* by Peter Hansbo and Mats G. Larson.

The code implements nodal and multiplier-free Reynolds discretizations,
stabilized Crouzeix--Raviart methods in two and three dimensions,
Taylor--Hood comparison calculations, and an undamped smoothed Newton method.
The regression drivers reproduce all ten manuscript tables and all seven
vector figures without random inputs.

## Requirements

- MATLAB R2020a or newer
- MATLAB R2025a for the reference run committed here

No additional MATLAB toolboxes are required by the numerical solvers.

## Quick start

From the repository root, run

```matlab
addpath('matlab')
reproduce('quick')
```

This checks Tables 1--6 and the quoted steep-pit mesh sensitivity. On the
machine used for the reference run it takes about one minute.

The available modes are

```matlab
reproduce('quick')    % Tables 1--6
reproduce('tables')   % Tables 1--10, including all 48 Newton cases
reproduce('figures')  % all seven manuscript figures
reproduce('full')     % all tables and figures
reproduce('certify')  % tables plus prose-claim and safeguard audits
```

The archival table run is intentionally expensive. The 48-case Newton matrix
dominates its runtime and took about 30 minutes on the reference machine.
CSV files are written to `matlab/results`; figures are written to `figures`.
Set the environment variable `CAVITATION_FIGURE_DIR` to choose another figure
directory.

## Repository layout

```text
matlab/             solvers, mesh utilities, and reproduction drivers
results/reference/  committed MATLAB R2025a reference CSV files
figures/            generated vector figures (not committed by default)
docs/               protocols, scope, and release notes
```

Each table driver contains numerical assertions against the committed values.
The extraction rules that matter for discontinuous pressure fields and cavity
lengths are documented in `docs/numerical_protocols.md`.
The complete MATLAB R2025a run is summarized in
`docs/verification_summary.md`.

## Reproducibility status

- Tables 1--10 pass their regression checks in MATLAB R2025a
- all 48 smoothed-Newton cases converge in the documented undamped mode
- all quantitative claims stated only in the manuscript prose pass their
  dedicated regression checks
- the safeguarded 48-case matrix reproduces the undamped peak pressures to
  within `4.72e-9` and the cavity lengths to roundoff
- Tables 6 and 10 use 8001 mid-plane samples on the stated interior window
- the seven figure drivers generate vector PDF files
- no random numbers or external data files are used

The optional residual-decrease safeguard in `solvedisccrn.m` backtracks in
four cases with `delta/r = 0.5` and `s = 1e-4`, but produces the same reported
peak pressures and cavity lengths. The manuscript and reference matrix use
undamped Newton steps. See `docs/numerical_protocols.md` for details.

## Citation

Use `CITATION.cff` for the software citation. Bibliographic details for the
paper can be added there after publication.

## License

The software is released under the BSD 3-Clause License. See `LICENSE`.
