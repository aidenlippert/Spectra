# H6: complete reduced construction from the Hamiltonian

The H6 spin-reduced certificate no longer depends on a previously constructed full-sector proof. A new run imports only the original rational Hamiltonian, mode count, and particle count. It constructs its own reference, upper witness, complementary blocks, factors, and response directions. The resulting global original-H interval has width **1.931992366627265e-10 Hartree**, using P32 and 28 response directions.

All stages stay within 400 spin-zero source configurations. There is no inherited P, upper, FCI vector, component partition, gap, or factor certificate. The constructor still explicitly enumerates the entire spin-zero sector and diagonalizes complement blocks as large as 200 for numerical proposals. It has not removed exponential growth.

## Construction and acceptance

1. Read only `modes`, `particles`, and `hamiltonian`. Optional input fields containing states, witnesses, response recipes, and factors are ignored. The numerical budgets are separate algorithm parameters.
2. Construct the nearby spin-symmetric Hs using the CAR Casimir projection. Verify its Hermiticity, particle conservation, and exact spin commutators. Recompute the coefficient-norm difference from the original H.
3. Start from the determinant with the lowest N orbital indices occupied. For even N this is spin zero. Repeatedly solve the small retained Ritz problem and add the determinant with largest external residual magnitude, up to P32. A vanished residual can stop selection early; it does not establish that the selected state is the ground state.
4. Starting from that selection's own integer witness, perform 20 sparse physical Krylov upper-refinement steps on Hs. Exact Rayleigh evaluation accepts the upper. No independent FCI witness is supplied.
5. Enumerate only spin-zero determinants. Discover the Q action graph and its connected components, then verify complete coverage and absence of intercomponent coupling.
6. Use the smallest numerical block eigenvalue minus 10^-4 Ha to propose gamma, rounded downward to a rational decimal. Require gamma above the independently evaluated upper. Construct new fixed-point LDL factors and verify their residual margins exactly. Numerical eigenvalues are never accepted as lower endpoints.
7. Generate optimized inverse-response directions using sparse Q actions. Reconstruct exact Gram matrices and accept every interval only through exact complement and Schur checks.
8. Transfer the Hs lower endpoint to the original H using the recomputed operator-norm allowance, and evaluate the final upper on the original H. Export a complete certificate only if its original-H width meets the requested target.

The original-H error allowance is 5.4e-11 Ha. The internal response target leaves room for two such allowances: one for transferring the lower and one for the possible difference between the two witness Rayleigh values. The final checker recomputes the actual original-H upper, rather than assuming that allowance is saturated.

The resulting proof uses the existing exact spin-reduction verifier. Its logical justification remains global SU(2) reduction followed by a P/Q split inside spin zero; P itself need not be spin adapted. The nearby-H perturbation allowance protects the original Hamiltonian even if it is not exactly spin symmetric.

## Recorded H6 result

| Quantity | Result |
|---|---:|
| Full fixed-number sector dimension | 924 |
| Enumerated spin-zero sector dimension | 400 |
| Selected P dimension | 32 |
| Complement block dimensions | 200, 168 |
| Own upper-refinement directions / support | 20 / 200 |
| Certified complement threshold | -6.26489910104 Ha |
| Response directions | 28 |
| Hs response interval width | 1.3919923647405865e-10 Ha |
| Original-H interval width | 1.931992366627265e-10 Ha |
| Largest numerical eigensolve | 200 |
| Final certificate bytes | 1,643,705 |

The recorded construction timer is about 223.08 seconds, starting after input validation and spin projection. It includes the own upper refinement, complement discovery and factorization, response proposals, repeated exact response replays, and the final original-H replay. This is not a benchmark against another chemistry package or an asymptotic speedup claim.

Construction uses 400 distinct source configurations under Hs and 200 under the original H; the latter are a subset of the same 400 configurations. Thus there are 600 distinct operator/source pairs, not merely 400 operator evaluations. Repeated acceptance checks evaluate some pairs again. Final replay reports these categories separately. On Hamiltonians that also break Sz conservation, the final original-H witness actions can reference additional configurations; the accounting explicitly permits and reports this.

## Tests that protect the provenance claim

A small exactly spin-symmetric hopping system is constructed twice, once with bare input and once with deliberately false P, upper, gap, block, and response metadata. Both runs produce identical certificates. During construction, tests intercept file reads to permit only the chosen input and newly generated output files. They reject the inherited-certificate reducer, check every action source is spin zero, constrain determinant enumeration to the spatial-orbital combination problem, and limit numerical eigensolves to the small test's three-dimensional blocks.

A disconnected test starts from the wrong exact eigenstate. Its selected reference cannot support a complement threshold above the upper, so construction records a failure and exports no final certificate. An oversized spin-zero sector is rejected before enumeration. These are genuine budget/refusal gates, not evidence of a general reference-selection theorem.

The direct H6 certificate is independently replayed with `python -S`, without numerical packages. Its original Hamiltonian matches the fixture exactly, and its final P matches the newly generated selection record.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m experiments.marginal_spin_constructor \
  --source results/marginal_h6/fixture.json \
  --output results/marginal_h6/direct_spin

python -S -m experiments.marginal_spin_reduction \
  --verify results/marginal_h6/direct_spin/certificate.json
```

The constructor preserves existing output directories; reproduction requires a fresh output path.

## A concrete transfer obstacle and a route around it

A separate probe adds spin-independent hopping of strength 1/1000 between spatial orbitals 0 and 1. It preserves SU(2), but joins the two complementary components into one 368-dimensional block. The current direct block-factor interface rejects that block because its maximum dimension is 256. This is an implementation/representation budget, not evidence that the physical gap disappears.

The perturbation has four CAR coefficients of magnitude 1/1000, so its exact coefficient-norm bound is eta=0.004 Ha. Independently replaying the unperturbed complement proof gives the rigorous transferred inequality

\[
Q(H_s+\Delta H)Q\succeq(\gamma-\eta)Q.
\]

The shifted threshold remains **0.06415952517180076 Ha above** the perturbed Hamiltonian's evaluated upper witness. This supplies a usable prospective complement bound without factoring the joined 368-state block. This initial probe records exact arithmetic and an independently checked reference gap. Complete perturbed intervals and regenerated response certificates are now supplied by the subsequent connected transfer implementation.

The subsequent [connected transfer implementation](marginal_h6_connected_transfer.md) now retains the certified reference blocks and regenerates the physical upper and response at strengths1/1000 and1/100, reaching original-H widths2.05690e-10 and1.54576e-10 Ha without factoring the joined368-state block. A separate exact four-square norm bound halves the conservative hopping norm; its use on a stronger H6 transfer remains to be tested. Beyond that finite transfer lie the unresolved growth of the reference, response, explicit spin sector, factor sizes, and integer precision. General marginal representability, the asymmetric quartic optimum, and a uniform chemistry compiler remain open.
