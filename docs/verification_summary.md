# Verification summary -- September 11, 2026

The September 10 revision completed `reproduce('certify')` on MATLAB R2025a,
including all current tables, both 48-case Newton matrices, the fine 3D
solve, Korn constants and prose checks. Independent direct/QP comparisons,
residual checks and a Korn energy assembly check support the new solvers.

Of 369 audited numerical table cells, 368 agree at their displayed precision.
The optional sign-rule return also has a stale `info.pnorm` diagnostic; it
does not affect the table solutions. Table 10 has one last-digit rounding discrepancy: 2.14714760459147 is printed
as 2.1472 rather than 2.1471. Fifty updated prose comparisons pass.

Unrecognized reproduction modes silently return without doing work; use the
documented mode names.

Seven figure drivers complete, with a recorded batch-export legend issue
in Figures 5 and 6 on MATLAB R2025a. Fresh standalone calls to the two
drivers export their legends correctly. Existing manuscript figures remain
unchanged. See the full audit for scope, measured errors and limitations.

The source set contains 59 MATLAB files, including eight new helpers/drivers;
52 CSVs preserve the manuscript reference summaries. Source and result
checksums are retained with the local audit; the code update preserves
Peter's numerical implementation.

Full detail: [numerical audit](numerical_audit_2026-09-11.md).
