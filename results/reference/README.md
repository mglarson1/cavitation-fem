# Reference summary data

These CSVs accompany the September 8, 2026 manuscript (Tables 1--14).
Root CSVs contain the original studies and QP/smoothing verification;
`domain_length_2026-09-08/` and `pressure_ends_2026-09-08/` contain the
new end-location and front studies. The old Table 2 data are superseded by
raw 2D errors; Git history preserves the earlier snapshot.

The `unstabilized` identifier in Table 8 means full-gradient CR without a
jump penalty. Its two `converged=0` rows are last-iterate diagnostics.
All current manuscript qualifications apply to these data. Timing and
roundoff fields are recorded observations, not exact regression targets.
Run the drivers to generate complete MAT fields and dense profiles under
`matlab/results/`; see the root README and numerical protocols.
