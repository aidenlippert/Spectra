# Automatic determinant references on two H4 Hamiltonians

The H4 rectangle now has an exact Schur interval of width **1.0076066967e-10 Hartree** using 26 selected determinants. Using the existing variational witness to guide selection reduces this to 17 determinants at the same width. The square does not reach the target under the 32-determinant budget: the best tested selector gives width **0.1304555622 Hartree**.

This is a finite reference-selection experiment. Construction and replay explicitly assemble the full 70-state, four-electron/eight-mode matrix. The result does not remove configuration enumeration or establish scalable molecular reference discovery.

## Exact certificate

Partition the rational Hamiltonian into retained determinants P and complement Q. Let

\[
\gamma=\min_{i\in Q}\left(H_{ii}-\sum_{j\in Q,j\ne i}|H_{ij}|\right).
\]

Gershgorin's theorem gives H_QQ>=gamma I. For b<gamma, write K=H_PQ H_QP. The sufficient lower certificate is

\[
H_{PP}-bI-\frac{K}{\gamma-b}\succ0.
\]

Indeed (H_QQ-bI)^(-1)<=I/(gamma-b), so the checked matrix is a conservative lower bound on the exact Schur complement. Rational LDL verifies positivity. The strict condition gamma>b is checked before division; a positive retained matrix with a negative denominator cannot certify anything through this argument.

The upper endpoint is recomputed from the integer amplitude witness in ascending bitstring order. Replay reconstructs H, P, Q, gamma, K, the Schur pivots, and the Rayleigh quotient. It does not trust stored energy or gamma fields. Invalid particle sectors, duplicate/nonphysical retained bitstrings, zero or malformed witnesses, and failed positivity are rejected.

## Selection comparison

Initially remove the Q row with the smallest exact Gershgorin lower bound, breaking ties by ascending bitstring, until gamma exceeds the valid upper by 0.01 Hartree. This needs five retained determinants for the rectangle and eighteen for the square.

Three rules then extend P, with a common cap of 32:

- **Row:** continue removing the worst Q Gershgorin row.
- **Coupling:** use a numerical lowest vector of the current scalar-Schur effective matrix and retain the Q determinant with largest squared coupling to it. Numerics only select candidates; exact positivity accepts bounds.
- **Witness:** prioritize coupling to the P restriction of the already available integer variational witness. This consumes more information than the coupling selector and must not be described as reference discovery from H alone.

| Fixture / selector | Best retained size | Final greedy size | Certified width (Ha) | Stop |
|---|---:|---:|---:|---|
| Rectangle / row | 30 | 32 | 0.03314526408676067 | Dimension cap |
| Rectangle / coupling | 26 | 26 | 1.0076066967035776e-10 | Accuracy target |
| Rectangle / witness | 17 | 17 | 1.0076066967035776e-10 | Accuracy target |
| Square / row | 31 | 32 | 0.2357547089902512 | Dimension cap |
| Square / coupling | 32 | 32 | 0.1304555622232512 | Dimension cap |
| Square / witness | 32 | 32 | 0.28635766589125117 | Dimension cap |

The best previous lower is preserved: replacing part of a scalar Schur bound by an exact retained block does not make this particular sequence of scalar approximations monotonically tighter. Several trial widths increase, as the histories record.

Observed complete constructions take about 0.49–1.11 seconds on this run; dependency-free replays take 0.04–0.09 seconds. These small-instance observations are not asymptotic performance evidence. All six exports have independent `python -S` replay receipts bound to the original rational fixtures. The source witnesses come from prior finite-basis molecular calculations; the old SOS lower certificates are not used to construct these new lower bounds.

The square exposes the next reference problem: a scalar complement bound can give poor lower energies even after the selected space is large. The next diagnostic separates leakage in the physical ground component from artificially low Schur estimates in other components, and tests component-specific complement bounds at the same retained size.

## Invalidated early probe

The first worker probe omitted gamma>b, hard-coded five retained indices, and trusted a stored upper. Its claimed 4.57e-24 width is invalid and withdrawn. `h4_receipt_invalidated.json` preserves that failure with an explicit invalid status. The subsequent five-state `h4_receipt.json` is a superseded diagnostic, not an accepted standalone certificate. Only the six fixture/selector subdirectories above contain accepted exports.

The current tests include the precise negative-denominator counterexample, an independent ordered-occupation CAR matrix oracle, all six accepted replays, false lower bounds, changed Hamiltonians, and malformed upper witnesses. Molecular intervals refer to the exported rational finite-basis electronic Hamiltonians; integral, basis, geometry, and model errors are separate.

The exact square diagnostic finds no missing upper-witness support in either final coupling- or witness-selected P. Normalized ||QH psi||^2 is 5.12148e-18 for the coupling selection and exactly zero for witness selection. This rules out omitted amplitudes of that witness as the explanation for the wide lower interval. The next attack is a component-specific complement bound; it must still certify every Q-only component. These figures are in `square_witness_leakage.json`.

The full marginal suite passes 153 tests in 49.337 seconds after these changes.

The subsequent [component-response experiment](marginal_molecular_component_response.md) closes the unperturbed square interval to 1.00251e-10 Ha using exact moments. It also tests connected hopping and identifies the remaining interblock approximation loss. The scalar results above remain the controls for that comparison.
