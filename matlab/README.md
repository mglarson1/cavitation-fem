# MATLAB reproduction guide

The files in this directory reproduce the finite element calculations and
vector figures of *Structure-Preserving Augmented Lagrangian Finite Element
Methods for Cavitation in Reynolds and Stokes Flows* (31 pages, Tables 1--12,
Figures 1--7, revision of September 10, 2026). The reference results were
computed with MATLAB R2025a and R2025b. No random numbers or external data
files are used. Optimization Toolbox (`quadprog`) is required for Tables
8--10 and 12.

From this directory, run

```matlab
reproduce('quick')          % Tables 1--4
reproduce('tables')         % Tables 1--12 and the prose checks of Sections 5.6 and 5.8
reproduce('figures')        % Figures 1--7
reproduce('korn')           % Korn constants of Remark A.6 (about 6 GB memory)
reproduce('full')           % tables, figures and Korn constants
reproduce('domain')         % Table 10
reproduce('pressure-ends')  % Table 9 and the flat Couette calibration
reproduce('verification')   % 2D Reynolds errors, QP references and end studies
reproduce('certify')        % tables, Korn constants and all prose checks
```

Table output is written as CSV files under `matlab/results`. In the
manuscript tree figures are written to `tex/figures`; in the standalone
release they are written to `figures`. Set `CAVITATION_FIGURE_DIR` to
override the destination. The 48-run smoothed Newton matrix dominates the
runtime of `tables`, about 30 minutes on the reference computer. The
36x12x12 three-dimensional solve and the Korn constants need about 6 GB of
memory.

## Tables and quantitative claims

| Manuscript output | Driver | Main CSV in `results/` |
|---|---|---|
| Table 1 | `reproduce_reynolds_tables.m` | `table01_reynolds_refinement.csv` |
| Table 2 | `reproduce_reynolds_tables.m`, `verify_reynolds_2d_convergence.m` | `table02_reynolds_convergence.csv` |
| Section 5.3, largest eigenvalues | `reproduce_reynolds_tables.m` | `table03_reynolds_stability_threshold.csv` |
| Table 3 | `reproduce_reynolds_tables.m` | `table04_reynolds_method_comparison.csv` |
| Table 4 | `reproduce_stokes_core_tables.m` | `table05_stokes_channel_refinement.csv` |
| Table 5 | `reproduce_stokes_3d_table.m` | `table09_stokes_3d_verification.csv` |
| Table 6 | `reproduce_constitutive_table.m` | `table07_constitutive_law_comparison.csv` |
| Section 5.6, mechanical pressure | `verify_mechanical_pressure.m` | `prose_mechanical_pressure.csv` |
| Table 7 | `reproduce_discretization_table.m` | `table08_discretization_comparison.csv` |
| Table 8 | `reproduce_model_comparison.m` | `model_comparison_pressure_normal_flow.csv` |
| Table 9, flat Couette calibration | `reproduce_pressure_ends.m` | `pressure_ends_2026-09-08/table14_front_refinement.csv`, `flat_channel.csv` |
| Table 10 | `reproduce_domain_length.m`, then `reproduce_end_sensitivity.m` | `domain_length_2026-09-08/table12_domain_length.csv`, `end_sensitivity.csv` |
| Table 11 | `reproduce_solver_table.m` | `solver_comparison.csv` |
| Section 5.8, 48 smoothed Newton runs | `reproduce_smoothed_newton_table.m` | `table10_smoothed_newton_all_48_cases.csv` |
| Table 12 and QP references | `verify_stokes_qp_reference.m`, `verify_smoothing_against_qp.m` | `stokes_qp_reference.csv`, `smoothing_against_qp.csv` |
| Remark A.6 | `verify_korn_constants.m` | `korn_constants.csv` |
| Other quantitative prose | `verify_manuscript_claims.m` | `prose_*.csv` |

The CSV names keep the numbering of the September 8 release, whose tables
were ordered differently; the table above gives the current correspondence.
`reproduce_stokes_core_tables.m` also writes the traction-free four-pit
comparison (`table06_reynolds_stokes_comparison.csv`) and the steep-pit
mesh sensitivity, which enter Tables 7 and 10 and the text of Section 5.7.
Each driver checks its values against the manuscript and stops with an
error on a regression.

## Figures

| Figure | File | Driver |
|---|---|---|
| 1 | `pit.pdf` | `makepit.m` |
| 2 | `pitlambda.pdf` | `makepitlambda.m` |
| 3 | `pitline.pdf` | `makepitline.m` |
| 4 | `channel.pdf` | `makechannel.m` |
| 5 | `stokes3dfine.pdf` | `makestokes3d.m` |
| 6 | `pitcompareends.pdf` | `makepitcompare.m` |
| 7 | `pitfieldends.pdf` | `makepitfield.m` |

All figure drivers call `configure_submission_figures.m`, which forces white
axes and black labels independently of the MATLAB desktop theme. In the
three-dimensional figure each boundary face is drawn once, since coplanar
faces drawn twice interleave in vector output.

## Solvers

The proposed-method solvers are `solvereynolds.m` (nodal mixed Reynolds),
`solvereynolds2.m` (multiplier-free stabilised Reynolds), `solvedisccrs.m`
(two-dimensional jump-stabilised Crouzeix--Raviart Stokes, with active-set,
smoothed Newton and quadratic programming options) and `solvedisccrs3.m`
(three dimensions, direct saddle-point solves). `solvedisccrs3_ip.m` solves
the same three-dimensional frozen systems by the iterated penalty method with
a nested dissection Cholesky factor and is used for the 36x12x12 mesh;
`box3dcase.m` sets up that geometry. `solvedisccr.m`, `solvedisccrn.m` and
`solvedisc.m` supply the full-gradient Crouzeix--Raviart, smoothed Newton and
Taylor--Hood comparisons.

The active-set iterations stop when the set repeats. `solvedisccr.m` accepts
an optional seventh argument `tol` for the sign rule of Section 4.1, which
stops as soon as a frozen solution has `p >= -tol` and `div u >= -tol` on every
element. It is used with `tol = 1e-12` for full-gradient Crouzeix--Raviart on
the two finer steep-pit meshes, where the set changes only by elements with
pressure and divergence at rounding level and never repeats.

## Pressure-normal-flow ends

`pressure_end_pit.m` and `run_pressure_end_case.m` prescribe the transverse
velocity and the normal traction at the Stokes ends (the boundary part
Gamma^E of the manuscript). The transverse Crouzeix--Raviart face mean is
zero, the normal face mean is free, and the boundary penalty acts only on the
prescribed component, through the `dirichlet_components` option of
`solvedisccrs.m`. These problems are solved with `solver='qp'`, a convex
quadratic program followed by working-set refinement and checks of
stationarity, signs and complementarity (`solve_cavitation_qp.m`).

Sampling for Tables 8--10 uses spacing 0.000375 on 4 < x < 20 and the upper
mid-plane trace of the Stokes pressure. A Reynolds interval is cavitated when
both endpoint nodes are active, a Stokes sample when its element is inactive.
The reformation front is the rightmost cavitated sample relative to x = 12.
The traction-free studies of Tables 6, 7 and 11 use 8001 samples on [0,24]
restricted to 4 < x < 20, and Taylor--Hood samples are cavitated when
`p <= 1e-8`.
