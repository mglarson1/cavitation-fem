# Numerical protocols

This file records conventions that are easy to lose when a pressure field is
discontinuous or a cavity is represented by an active set.

## Reynolds verification

Table 1 uses the nodal active-set solution and reports the fraction of active
nodes. Table 2 uses a one-element-wide strip for the manufactured
one-dimensional obstacle problem. The two pressure traces are averaged to
cancel the equal and opposite diagonal-splitting errors, and the resulting
piecewise linear error is integrated with high-order quadrature. Table 3
detects quadratic Dirichlet nodes from boundary-edge topology rather than
from endpoint flags, which would incorrectly classify two interior diagonal
midpoints.

## Two-dimensional channel

The CR/P0 pressure in Table 5 has two traces on the mesh-aligned centerline.
The reported pressure is the upper trace sampled at longitudinal cell centers
for `x < 2.5`. The cavity fraction is the exact fraction of inactive elements
by area.

## Common-pit comparison

Tables 6--8 use 8001 uniformly spaced physical mid-plane samples on
`0 <= x <= 24` and retain only `4 < x < 20`. Reynolds intervals are
cavitated when both endpoint nodes are active. CR/P0 Stokes samples are
cavitated when their containing elements are inactive. Taylor--Hood samples
use the threshold `p <= 1e-8` because complementarity is weak in that space.

The steep-pit mesh ratios compare every Stokes mesh with the fixed Table 6
Reynolds cavity length `6.219`, exactly as stated in the manuscript.

## Three-dimensional pressure

The CR/P0 pressure in Table 9 has no unique value on the mesh-aligned line
`y = z = 0.5`. The peak is the maximum over the one-cell tube surrounding
that line, restricted to `x < 2.5`. The cavity fraction is computed by
tetrahedral volume; these structured meshes have uniform volumes.

## Smoothed Newton method

Table 10 uses the central-path cavity criterion
`p < sqrt(s/gamma)` and 8001 mid-plane samples on `4 < x < 20`. The complete
matrix contains four meshes, six values of `s`, and both pits with
`delta/r = 0.5` and `1`, for 48 cases.

The reference calculations use undamped Newton steps from the prescribed
boundary velocity and zero pressure. All 48 converge. The optional
residual-decrease safeguard accepts full steps in 44 cases. In the four
`delta/r = 0.5`, `s = 1e-4` cases it uses seven to nine backtracks, while the
windowed peak differs from the undamped result by at most `4.8e-9` and the
reported cavity length is unchanged. The central-path identity is checked to
the nonlinear solver tolerance and its absolute defect is stored in the CSV
output.
