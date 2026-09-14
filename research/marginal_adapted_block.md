# Adapting the physical block to the connected-chain objective

The adapted block gives a strictly better certified upper energy for the
half-filled open Hubbard chain with U=4 and t=1. At one million sites the exact
certificate implies the outward-rounded interval

\[
-0.611636\le E_0/N\le -0.55820017099.
\]

The exact interval width is approximately 0.0534358290063282 per site. The
preceding fixed-block filter gave an upper density of -0.5581139722444636;
the improvement is approximately 0.0000861987492082 per site. This closes a
specific variational optimization step, while leaving a substantial gap to
the lower certificate.

## The effective operator

For an eight-site block, write T_L,T_R for its positive internal endpoint
hopping operators and D_L,D_R for endpoint doublons. Let

\[
Q_i=D_i-\tfrac12 n_i+\tfrac12 I,\qquad
w=\frac{\eta^2}{1+\eta^2},
\]

\[
K_\eta=H_8+w\left[\tfrac12(T_L+T_R)-8(Q_L+Q_R)\right].
\]

This is an actual number-conserving operator, with particle-hole symmetry
modulo fixed total number. The verifier checks its CAR expression and signed
symmetries independently. In the preserved half-filled symmetry sector,
each endpoint spin occupation equals one half, so

\[
e_8+g_\eta=\langle K_\eta\rangle+4w-
\frac{2\eta}{1+\eta^2}.
\]

This is the energy per block in the infinite-chain limit. A finite chain of
q blocks retains the exact count of q-1 interfaces:

\[
E_q=q e_8+(q-1)g_\eta.
\]

Optimizing K_eta therefore accounts for both the block interior and the
energy of joining its boundaries. It need not minimize a finite chain's
energy; the replay uses the finite formula without making that claim.

## The accepted state and its restricted limit

The new physical state is

\[
\phi=\frac{p(K_\eta+4I)\phi_0}{\|p(K_\eta+4I)\phi_0\|},
\qquad \eta=0.230522,
\]

where phi_0 is the independently validated original eight-site state.
The nine rational coefficients of p are frozen in
`results/marginal_graded_hubbard8/adapted_block/certificate.json`.

Exact Horner recurrence constructs a nonzero state. Its original-H energy,
effective-K energy, and two endpoint pair RDMs are recomputed. An independent
four-site CAR contraction verifies the filter contribution from those RDMs.
The resulting values are approximately

\[
e_8=-4.235081379150620,\qquad
g_\eta=-0.230521832973418.
\]

The isolated block energy rises slightly, while joining the blocks recovers
more energy. Their sum improves.

At this specific eta, form the physical vectors
v_j=(K_eta+4I)^j phi_0 for 0<=j<=8. Exact integer arithmetic proves

\[
[\langle v_i|K_\eta|v_j\rangle]
-(-4.22965897)[\langle v_i|v_j\rangle]\succeq0.
\]

This 9-by-9 PSD certificate bounds every nonzero Rayleigh quotient in this
span. The achieved thermodynamic energy density is within approximately
2.17316116414e-9 of that restricted lower limit. It does not bound arbitrary
block states, larger polynomial spans, different eta, or the full chain's
ground energy. The full-chain lower bound is separately supplied by the
fresh all-Fock six-site window certificate.

Numerical discovery used a 1,260-dimensional signed-orbit matrix. Exact
replay also reaches 1,260 orbit amplitudes. This is a substantial finite
eight-site calculation, not evidence that a general many-body search has
collapsed. The global chain state and its exponentially large normalization
are never expanded.

## Why assembly remains physical

The polynomial preserves definite total spin populations and the checked
particle-hole symmetry. The contact occupations are checked exactly after
adaptation. The existing linear-filter merging proof therefore applies.

For a remote even, spin-number-preserving observable O, the active cut
hopping satisfies

\[
\langle hO+Oh\rangle=0,\qquad
\langle hOh\rangle=\langle O\rangle.
\]

The second identity uses the partner contact's half occupations and the
orthogonality of directed hopping sectors. For F=I-eta*h, division by the
exact norm 1+eta^2 preserves the remote expectation. This is a specific
nonunitary-filter identity; disjoint support alone would not preserve a
reduced state. At each merge the entire previously joined cluster supplies
the fixed total spin populations needed for the next step.

## A proved obstruction for the next filter family

The natural extension F=I-a*h+b*h^2 introduces a new fourth-moment term.
For remote O, [O,h]=0, so the terms linear in b still reduce through the
existing second-moment identity. The b^2 term does not generally do so.

If the partner contact has half occupation in each spin and doublon
probability d, exact CAR contraction gives

\[
\operatorname{Tr}_B(\rho_B h^4)
=(4-6d)I+(12d-3)C_A,
\qquad C_A=1-n_{A\uparrow}-n_{A\downarrow}+2D_A.
\]

The operator is scalar for arbitrary A precisely when d=1/4. Otherwise a
quadratic filter can change remote observables correlated with contact charge.
For the explicit two-site block states, written with ascending occupation
bit labels,

\[
\phi_A\propto |3\rangle+|12\rangle+|9\rangle-|6\rangle,
\quad \phi_B\propto |9\rangle-|6\rangle,
\]

both contact spins initially have occupation one half. Applying I+h^2 at
the middle interface changes the remote A-site doublon probability from
1/4 to 2/11. The unnormalized product and filtered norms are exactly 8 and
44. `discovery/quadratic_filter_obstruction.py` replays both the partial
operator identity and this counterexample.

Thus the next broader filter must propagate a boundary state or functional
that retains these correlations. Reusing only the scalar merge shift would
produce incorrect energies. A finite boundary transfer contraction is the
next concrete target; its implementation, certified rounding, accuracy and
cost remain to be established.

An initial exact diagnostic already validates a 16-by-16 norm transfer,
e^T G(BG)^(q-1)e, against direct determinant CAR application. Here G is the
block's endpoint double-layer tensor, B reshuffles F-dagger F, and e traces
the two unacted outer sites. Twelve comparisons cover one through four
two-site blocks and linear, quadratic, and mixed filters. Ordinary local
occupation tensors suffice: adjacent hopping's preceding Jordan-Wigner
strings cancel as an operator identity. This norm-only diagnostic does not
yet implement energy insertions or certify large-chain numerical rounding.
Its replay is `discovery/boundary_norm_transfer_probe.py`.

## Reproduction and validation

```sh
OPENBLAS_NUM_THREADS=1 python -S results/marginal_graded_hubbard8/discovery/adapted_block.py
python -S results/marginal_graded_hubbard8/discovery/quadratic_filter_obstruction.py
```

The adapted certificate replay took 66.960 seconds, including a fresh
six-site lower proof, the original filtered baseline, exact polynomial-span
PSD verification, and independent replays at 8, 16, 24, 64 and one million
sites. The previous rotation-family and transferred-target certificates
were also freshly replayed after the shared endpoint reader was refactored.

The focused 24 tests passed in 31.613 seconds. They include a nonconstant
polynomial checked using direct determinant CAR action, the exact matching
million-site baseline, lower-bound refusal gates, and prior two- and
three-block physical assembly checks. All 449 regression tests passed in
320.018 seconds. The 13 adapted-replay source/input hashes and 11 freshly
replayed rotation/transfer hashes matched afterward. Full validation results
are recorded in `results/marginal_final_validation.json`.

General N-representability, accuracy-versus-cost guarantees, molecular
Coulomb transfer, thermal and dynamic prediction, and synthesis remain open.
