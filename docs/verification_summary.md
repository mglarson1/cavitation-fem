# Verification summary

The reference calculations were run with MATLAB R2025a on 2026-09-01. Every
mesh and parameter configuration used in Tables 1--10 was executed, including
the complete 48-case smoothed-Newton matrix. The standalone checkout also
passed `reproduce('quick')` from its own directory. A subsequent certification
pass verified the quantitative claims stated only in the prose and reran the
complete 48-case matrix with the optional Newton safeguard.

## Checked outputs

- Tables 1--4: Reynolds refinement, manufactured-solution convergence,
  stability threshold, and comparison of the two Reynolds formulations
- Tables 5--6: two-dimensional Stokes refinement, Reynolds--Stokes comparison,
  and steep-pit mesh sensitivity
- Table 7: all eight Taylor--Hood constitutive-law runs
- Table 8: all nine discretization and mesh combinations
- Table 9: all three three-dimensional meshes
- Table 10: all 48 mesh, smoothing, and geometry combinations, plus four
  focused undamped-Newton checks
- Figures 1--7: generated as vector PDF files
- Prose-only checks: Reynolds KKT and complementarity residuals, gamma
  invariance, projected Gauss--Seidel agreement, stabilization sensitivity,
  active-set cycling, the fine gentle-pit comparison, three-dimensional gamma
  invariance, and the Newton residual history

The regression drivers check convergence, iteration counts, peak pressures,
cavity measures, sign conventions, complementarity or central-path defects,
and the mesh-specific values quoted in the manuscript.

## Corrections made during the audit

The executable audit identified two manuscript-facing numerical corrections.
The manufactured-solution errors in Table 2 were replaced by values obtained
from the documented high-order integration rule. In Table 8, the
unstabilized computation on the 2048-element mesh converges in nine active-set
iterations; failure after 100 iterations begins on the two finer meshes.
The same iteration count was corrected in the active-set column of Table 10.

Several previously implicit extraction conventions are now fixed in
`numerical_protocols.md`, notably the upper pressure trace in Table 5, the
one-cell pressure tube in Table 9, and the central-path cavity threshold in
Table 10.

The prose audit also made the Reynolds active-set location precise by stating
the mesh level and using the x-projection of all active nodes. It replaced an
unsupported claim of a final Newton residual order of 1.97: the certified
superlinear phase has observed orders 1.37--1.59, followed by one 2.71
estimate and tolerance-dominated values.

## Newton-mode qualification

The reference matrix uses undamped Newton steps, as stated in the manuscript.
With the optional residual-decrease safeguard, 44 of 48 cases accept full
steps throughout. The four gentle-pit cases with `s = 1e-4` use seven to nine
backtracks. The complete safeguarded matrix showed windowed peak-pressure
differences no larger than `4.72e-9` and cavity-length differences no larger
than `2.22e-16`. The reference mode remains undamped because that is the
algorithm stated in the manuscript.
