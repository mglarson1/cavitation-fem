# Numerical protocols -- September 10 manuscript

The authoritative table/driver/CSV mapping is in `matlab/README.md`. The
manuscript now has 12 tables; the previous 14-table arrangement is historical.

- Tables 1--3 retain the nodal mixed Reynolds, raw two-dimensional manufactured
  solution, and multiplier-free stabilized comparison. The stability
  eigenvalues formerly in Table 3 now appear in Section 5.3.
- Table 4 retains the 2D jump-stabilized deviatoric Stokes channel. Its sampled
  pressure uses the upper trace at longitudinal cell centers with x < 2.5.
- Table 5 adds the 36x12x12 tetrahedral grid. `solvedisccrs3_ip` assembles the
  same forms as the direct solver and uses iterated penalty solves with r=1e4
  and nested-dissection Cholesky factors. The pressure peak is over a one-cell
  tube, not the offset line used by Figure 5.
- Table 6 compares constitutive laws with Taylor--Hood. The separate mechanical
  pressure check uses stabilized CR, so its cavity counts and lengths are
  not the Table 6 values. Mechanical pressure is p-(2/3+lamfac)*div(u), mu=1.
- Table 7 applies a sign tolerance of 1e-12 to full-gradient CR only on the
  two finer meshes. The coarse full-gradient row retains repeated-set
  stopping. Table 11 applies the sign rule on all three meshes, explaining
  eight versus nine iterations for the same coarse problem.
- Tables 8--10 use pressure-normal-flow ends: prescribed transverse velocity
  face means and free normal means, with penalty on the prescribed component.
  The comparison uses a fixed 4 < x < 20 window, step 0.000375 and the upper
  Stokes mid-plane trace. Reynolds intervals are cavitated when both endpoint
  nodes are active; Stokes samples use the inactive pressure set. Fronts are
  the rightmost cavitated sample relative to x=12; these are sampled fronts,
  not exact interfaces. Table 10 retains traction-free comparisons separately.
- Table 11's smoothed Newton runs use s=1e-8 and residual tolerance 1e-10.
  Table 12 instead uses residual tolerance 1e-14; this explains the differing
  iteration counts. Its pressure error uses an element-area weighted norm
  against the independently assembled QP reference. Small central-path
  products are not resolved in relative terms even when sampled cavities agree.
- The Korn experiment assembles discontinuous P1 trace-free strain and full
  jump energies on the stated box and walls/inflow boundary subset. Its
  finite-grid eigenvalues are lower bounds on the relevant constants.

Figures 3, 5, 6 and 7 were updated. Figure 7 averages CR vertex evaluations
for visualization and draws an interpolated zero-velocity contour; this
postprocessing is not used for the pressure/cavity table diagnostics.
Reversed flow was independently confirmed in the raw CR field by the audit.
Figure 5 samples at y=z=0.501, avoiding the mesh-aligned trace ambiguity.

Legacy default-stopping diagnostics and the original 48-case matrix are
retained. Repeated-set nonconvergence in those checks is distinct from the
new sign-rule solutions reported in Tables 7 and 11.
