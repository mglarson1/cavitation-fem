# Manuscript reference summaries

These 52 CSV files accompany the September 10, 2026 manuscript (12 tables).
See `../../matlab/README.md` for the current table-to-file mapping: legacy
CSV filenames retain earlier table numbers.

The September 11 repository refresh updates Table 7's full-gradient sign-rule
rows and Table 5's fine 3D row, and adds `solver_comparison.csv`,
`model_comparison_pressure_normal_flow.csv`, `end_sensitivity.csv`,
`prose_mechanical_pressure.csv` and `korn_constants.csv`.

The supplied summaries are preserved rather than replaced with rerun timing
or roundoff values. The Table 10 value 2.14714760459147 is retained faithfully;
the manuscript prints 2.1472, one last-digit rounding discrepancy documented
in the audit. No reference value was changed to silence a regression.

Dense profiles, MAT fields and generated figures are not stored here.
Reproduction drivers regenerate them under the documented output directories.
