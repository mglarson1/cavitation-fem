# Numerical audit of Peter's September 10 revision

Audit date: September 11, 2026. Manuscript: the accompanying TeX manuscript and its
31-page PDF, with numerical examples v13 and Korn appendix v5. This audit
preserved the manuscript, numerical source files, figures and saved results.
It concerns the implemented numerical experiments; it is not a new proof
review of all manuscript theorems.

## Conclusion and findings

No defect affecting the reported solutions was found in the tested revised
experiments. A stale returned diagnostic is documented below. The full
MATLAB R2025a certification completed, and independent assembly,
solver, residual and raw-field checks support the new results. This is
verification of the stated experiments, not a guarantee for arbitrary
geometries, tolerances or solver options.

Four minor issues remain, recorded without editing Peter's manuscript or code:

1. **Table 10, inlet-only extension, traction-free Stokes peak:** the table
   prints 2.1472; the saved and fresh value is approximately 2.14714760459147,
   which rounds to **2.1471** at four decimals. The driver's 6e-5 regression
   tolerance accepts the printed number, but strict displayed-precision
   auditing does not. The difference is about 5.24e-5 and does not alter the
   end-sensitivity conclusion. All other 368 audited numerical table cells
   agree at the printed precision.
2. **Batch PDF export:** `reproduce('figures')` completed, but the MATLAB
   R2025a batch exports of Figures 5 and 6 omitted the visible legends.
   Fresh standalone calls to `makestokes3d` and `makepitcompare` produced the
   correct legends. The existing manuscript figures have the legends and
   were preserved. Inspect batch-generated PDFs before replacing them;
   use standalone driver runs for these two figures. No numerical values
   changed in these checks. The precise rendering cause was not established.
3. **Stale sign-rule diagnostic:** in `solvedisccr.m`, the optional sign-rule
   break (line 142) precedes the update of `pnorm` (line 145), so returned
   `info.pnorm` describes the previous iteration. On 2048 elements it is
   0.00270502762613 instead of the actual final change 2.17785994712e-5;
   on 8192 elements it is 7.00696056301e-5 instead of 9.20295886877e-11.
   Pressure, velocity, iteration count, sign-based stopping and manuscript
   diagnostics are unaffected. The field should not be read as the final
   pressure change when the sign rule fires. The production source was
   preserved; `sign_rule_diagnostic.csv` records the independent check.

4. **Silent invalid reproduction mode:** `reproduce.m` no longer checks
   that its mode is recognized. A call to `reproduce('audit_invalid_mode')`
   returns successfully without performing work. A typo can therefore look
   like a successful run. This was confirmed in MATLAB R2025a; all documented
   modes used in this audit executed normally. No production fix was made.

## Changed code and experiment mapping

The update comprises eight changed and eight new MATLAB sources. The
standalone snapshot contains 59 MATLAB sources and 52 reference CSVs.

- `solvedisccr` gains the optional sign rule. Default repeated-set behavior
  is preserved; Tables 7 and 11 deliberately use different stopping rules
  on their coarsest full-gradient row, giving nine and eight iterations.
- `solvedisccrs3_ip` changes the frozen linear-system solver, with the same
  element and jump matrices as the direct 3D implementation. `box3dcase`
  supplies the shared geometry and face-mean boundary data.
- `pressure_end_pit`, model and end-sensitivity drivers consistently use
  transverse essential velocity and normal natural traction at the ends.
  The transverse penalty mask agrees with the prescribed components.
- The mechanical-pressure driver uses the stated three-dimensional stress
  trace with mu=1. Its stabilized CR numbers are distinct from Table 6's
  Taylor--Hood constitutive-law comparison.
- The Korn driver evaluates the discontinuous P1 quotient, conforming
  quotient and H10 control on the stated box and boundary subset.
- Figures 3, 5, 6 and 7 are updated. CSV filenames intentionally retain
  their historical numbering; `matlab/README.md` maps the 12 current tables.

## Executed verification

1. A fresh `reproduce('certify')` completed in an isolated copy using
   MATLAB R2025a and Optimization Toolbox. This includes all 12 tables,
   both 48-case Newton matrices, the finest 3D case, Korn constants,
   pressure-end/domain audits and the existing prose checks. No saved
   working result was overwritten. Legacy repeated-set cycling checks are
   retained as diagnostics and do not invalidate the new sign-rule rows.
2. All 369 numerical table cells were compared against the saved datasets
   and fresh outputs: 368 pass strict rounding; the single Table 10
   discrepancy above is identical in both audits. Fifty updated prose
   comparisons pass against both datasets.
3. The new 3D solver was compared with direct solves on 6x2x2, 9x3x3 and
   18x6x6 meshes, with penalty parameters 1e2 and 1e4. Assembled matrices
   agree exactly. The maximum velocity difference is below 1.1e-11 and
   pressure difference below 1.9e-10; free-velocity stationarity residuals
   are below 6.9e-12.
4. The finest 36x12x12 iterative solve was checked separately. Stationarity
   is 8.90e-12, minimum divergence -2.47e-14, minimum pressure zero and
   complementarity 1.44e-13. Figure 5's sampled departure from 6-2x is
   0.22797309, consistent with its printed 0.23.
5. Sign-rule solutions on 2048 and 8192 elements were compared against the
   independently assembled QP reference. The matrices agree exactly;
   maximum pressure differences are below 7.1e-13, stationarity below
   2.9e-15, and every sampled cavity classification agrees. Peaks reproduce
   1.98695238 and 2.09338476 to the eight digits quoted in the manuscript.
6. Korn matrices on the 9x3x3 mesh were checked by independent direct
   integration of a deterministic piecewise affine field. Relative energy
   errors are 6.39e-16 for the strain/jump form and zero for the gradient
   form. The eigenvalue 6.3742156807 has generalized residual 3.28e-12.
   The full four-mesh Korn driver also passes. Growing finite-grid values
   are lower bounds and do not themselves establish an h-uniform bound.
7. Reversed flow was verified before graphical averaging: minimum raw CR
   vertex velocity -0.014119, minimum element mean -0.013529 and minimum
   face velocity -0.013654. Figure 7's averaged contour is visualization,
   not a free-boundary or error-estimation diagnostic.
8. All seven figure drivers completed in the audit folder. Updated numerical
   pages and figure outputs were inspected. Figure 3 is pixel-identical
   to its manuscript source at the checked rendering scale. Legends in
   standalone Figures 5 and 6 are correct; the batch issue is recorded above.
9. Source dependency analysis and inspection of dynamically evaluated calls
   support the curated standalone source set. A separate quick regression
   and dependency check also passed from that curated copy. MATLAB and Optimization
   Toolbox are required. The active PDF has 31 pages; its existing build
   log contains no LaTeX errors, unresolved references or layout warnings.

## Preservation and repository update

SHA-256 checks confirm that original manuscript files, MATLAB files and
saved results were unchanged by this audit. No section version or wrapper
was edited; no submission package was replaced. The previous submission
package remains the September 8 snapshot and does not represent this revision.

The GitHub update copies Peter's numerical changes and supplied reference
summaries, updates the table mapping and release documentation to version
0.3.0, and includes this audit and compact evidence CSVs. Existing reference
values are not adjusted to force agreement. Repository synchronization is
recorded separately after the commit is pushed; no visibility change or
journal submission is part of this audit.

Detailed logs, probe code, source checksums and visual checks are retained
in the local `review/audit_2026-09-11/` folder. The independent probes expose
assembled forms in audit-only copies; they do not modify production sources.
