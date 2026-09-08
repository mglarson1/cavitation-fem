# Structure-preserving finite element methods for cavitation

MATLAB implementation for *Structure-Preserving Augmented Lagrangian Finite
Element Methods for Cavitation in Reynolds and Stokes Flows*, by Peter Hansbo
and Mats G. Larson. This 0.2.0 review-stage snapshot accompanies the
38-page manuscript revised on September 8, 2026, with Tables 1--14 and
seven figures.

The implementation includes nodal and multiplier-free Reynolds methods,
jump-stabilized deviatoric Crouzeix--Raviart Stokes methods in two and three
dimensions, Taylor--Hood and full-gradient comparisons, smoothed Newton,
and checked unsmoothed quadratic-programming references.

## Requirements and commands

The reference environment is MATLAB R2025a. The original finite element
solvers use base MATLAB. Optimization Toolbox (`quadprog`) is required for
the QP verification and pressure-end studies, including `tables`, `full`,
`verification`, `pressure-ends` and `certify`. Earlier MATLAB versions are
not certified for this snapshot. No random numbers or external datasets
are used.

From the repository root:

```matlab
addpath('matlab')
reproduce('quick')          % Tables 1--6 and steep-pit mesh sensitivity
reproduce('tables')         % Tables 1--14, including 48 Newton cases
reproduce('figures')        % All seven figures
reproduce('domain')         % Original end-location study, Table 12
reproduce('pressure-ends')  % Couette calibration and fronts, Tables 13--14
reproduce('verification')   % Raw 2D Reynolds, QP, smoothing and end checks
reproduce('full')           % All tables and figures
reproduce('certify')        % Tables plus prose and safeguard checks
```

Generated data go to `matlab/results/`, which is ignored by Git. Reference
CSV files are in `results/reference/`. Runtime and roundoff-level residuals
can differ between runs; compare physical diagnostics and stated tolerances
rather than demanding identical file bytes. Figure output goes to `figures/`;
set `CAVITATION_FIGURE_DIR` to override that location. Full certification is
expensive: it includes the original and safeguarded 48-case Newton matrices,
three-dimensional solves and the new local/global refinement studies.

## Interpreting the results

Tables 6--9 retain the original boundary conditions. Table 8 includes two
explicitly unconverged full-gradient active-set rows; their last iterates
are not accepted nonlinear solutions. Table 11 uses a checked unsmoothed
QP reference to separate pressure error from cavity error. Small central-
path products are not resolved in relative terms even when cavity samples
agree with the reference.

Table 12 demonstrates end-location sensitivity with fixed local meshes and
observation window. Tables 13--14 prescribe transverse velocity and normal
traction at Stokes ends. These conditions reproduce exact flat Couette flow;
under regular traces they impose the same ambient pressure as Reynolds.
The resulting cavity comparison remains sensitive to domain ends and mesh
alignment. An unchanged front on nested grids is not a free-boundary error
estimate. This numerical boundary variant is distinct from the boundary
spaces analyzed in the manuscript's original Stokes theorem.

See `docs/numerical_protocols.md`, `docs/verification_summary.md`, and
`docs/release_notes.md` for definitions, evidence and remaining limitations.
The original annotated manuscript and exploratory scripts are not part of
this curated code snapshot. MAT fields and dense profiles are regenerated
by the drivers; only summary reference CSVs are committed.

## Citation and license

Use `CITATION.cff` for software attribution. The review-stage repository is
<https://github.com/mglarson1/cavitation-fem>; publication metadata can be
added after acceptance. The code is distributed under the BSD 3-Clause
License in `LICENSE`.
