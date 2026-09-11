# Cavitation finite element reproduction code

This repository accompanies *Structure-Preserving Augmented Lagrangian Finite
Element Methods for Cavitation in Reynolds and Stokes Flows*, by Peter Hansbo
and Mats G. Larson. Version 0.3.0 follows the September 10, 2026 manuscript:
31 pages, 12 tables and seven figures.

The MATLAB sources implement the nodal Reynolds method, jump-stabilized
Crouzeix--Raviart Stokes methods, comparison discretizations and numerical
verification drivers. The September 11 audit preserves Peter's numerical
code; it adds verification records and updates the repository documentation.

## Running the code

Use MATLAB R2025a or R2025b. Optimization Toolbox (`quadprog`) is required
for the QP reference and pressure-normal-flow studies. The finest 3D solve
and Korn calculation each need approximately 6 GB of memory. No external
data files are needed to run the computations.

From the repository root:

```matlab
addpath('matlab')
reproduce('quick')       % Tables 1--4
reproduce('tables')      % Tables 1--12 and associated comparisons
reproduce('figures')     % Seven figures; see the export note below
reproduce('korn')        % Numerical Korn constants
reproduce('certify')     % Full tables, Korn constants and prose checks
```

See [the reproduction guide](matlab/README.md) for the complete mapping of
current tables to drivers and CSV files. The historical CSV filenames keep
their original table numbers; they are not the current manuscript numbering.
Generated CSVs and fields go to `matlab/results/`; figures go to `figures/`
in this standalone repository. `CAVITATION_FIGURE_DIR` overrides figure output.
The 48-case Newton matrices dominate the full certification runtime.

## Reference results and audit

`results/reference/` contains the manuscript's reference summaries, including
the updated sign-rule, fine 3D, mechanical-pressure, end-condition and Korn
results. Runtime and roundoff-sized fields need not match byte for byte on
another computer. See [verification summary](docs/verification_summary.md)
and [September 11 audit](docs/numerical_audit_2026-09-11.md) for independent
solver checks and known discrepancies.

Four qualifications accompany this snapshot. Table 10 prints 2.1472 where
the supplied value 2.14714760459147 rounds to 2.1471. Also, batch figure
export in MATLAB R2025a omitted legends in two audit outputs; the manuscript's
existing figures have the legends. Inspect regenerated figure legends before
using the PDFs. These issues were not corrected in the numerical code
or manuscript during the audit. The optional sign-rule exit also leaves
`info.pnorm` at the previous iteration's value; it is not the final pressure
change in that case and is not used by the reported table diagnostics.
Finally, unrecognized `reproduce` modes silently do nothing; use the documented
mode names. All documented modes exercised in the audit completed.

The numerical Korn values are finite-mesh lower bounds, not a numerical
proof of uniform coercivity. Cavity predictions retain the documented mesh,
end-condition, observation-window and classification qualifications.

## Citation and license

Citation metadata is in `CITATION.cff`. The code is distributed under the
BSD 3-Clause license in `LICENSE`. The existing repository remains private
during manuscript review; this update does not register a DOI or change
repository visibility.
