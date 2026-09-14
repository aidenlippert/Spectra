# A one-body commutator extension of the Hamiltonian-derived dictionary

The exact H4 dual obstruction applies to the depth-zero creator-channel
dictionary. It does not apply after adding new generators. The dual checker
explicitly refuses proposals claiming a nonzero `one_body_steps` depth.

Root implemented an optional extension in `commutator_dictionary.py`:
start from each cubic creator-channel piece of `[H,a_i]`, then add the cubic
part of `[h1,p]`, where `h1` is the degree-two part of the original H. Keep
the quadratic baseline, initial pieces and adjoints. Repeated applications
form a bounded operator Krylov dictionary. Default depth zero is unchanged.
No physical state, old certificate, dual witness or upper bound is used to
construct the new generators. Every square coefficient through degree six
is matched. A focused algebra test verifies a genuinely new direction, not
just a repeated old dictionary.

The **actual depth-one H4 run** is
`results/certificate_scaling/commutator_krylov/h4_true_one_step/`:

| Quantity | Result |
|---|---:|
| Polynomial Gram entries | 8,304 |
| Symbolic word-pair products | 270,364 |
| Coefficient-map nonzeros | 274,580 |
| Construction / solve / export | 6.653 / 3.033 / .649 s |
| Total discovery/export | 10.335 s |
| Exact interval width | 4.89771659e-5 Ha |
| Exact exported factor rows / nonzeros | 311 / 10,922 |
| Certificate bytes | 387,865 |
| Factor denominator bit length | 109 |

Independent `python -S` two-sided replay passed .0016 Ha in .422 seconds.
Floating Clarabel status was `optimal_inaccurate`; it was not the acceptance
criterion. The earlier depth-zero creator run failed at .003564 Ha.

This is an accuracy gain from new H-only generators. It is **not yet a
performance win**: although the Gram matrix has fewer entries than the full
cubic dictionary's 17,776, the polynomial coefficient map is much denser
than that dictionary's 10,256 nonzeros and H4 total wall time is larger than
the previous full-cubic run. Large rational denominators also increase
certificate bytes. All these costs must count in the scaling test.

A depth-one H6 run is active on warm host A under an isolated source snapshot
`results/lambda_runs/scs_refinement/krylov_source.tar.gz`. Its 120-second solver
budget is separate from map construction/export; an external 900-second cap
bounds the experiment. No H6 result is claimed yet.

The agent's earlier directory `h4_one_step_v2` used neither an implemented
Krylov extension nor the new CLI flag. Its receipt is annotated as a repeated
depth-zero run and is excluded from this comparison.

Open question: can low depth yield accurate certificates on larger active
spaces while products/maps remain economical? A positive finite sweep alone
would still not prove an efficiently checkable sufficient condition or a
general scaling theorem.

## A support-preserving diagonal driver

Root also tested `--diagonal-driver`, using only the diagonal of the original
one-body H in the same orbital basis. For `D=sum_i e_i n_i`, the commutator
with a canonical word is exactly the same word multiplied by
`sum_(c,i in w) (2c-1)e_i`. This avoids support growth and has an independent
exact CAR comparison test. It is a basis-dependent proposal, not a discarded
part of the physical Hamiltonian: the certificate still matches the full H.

H4 (`h4_diagonal_one_step`) passed independent `python -S` replay at interval
width **6.51245170e-5 Ha**. Gram entries remained 8,304, while word-pair work
dropped from 270,364 to 120,892 and map nonzeros from 274,580 to 167,104.
Construction/solve/export were 3.231/2.790/.585 seconds, total 6.605 seconds.
The certificate has 311 factors, 10,930 nonzeros, 388,289 bytes and a 109-bit
factor denominator. This is one timing sample, not a robust speedup claim.
The next bounded comparison is the same diagonal driver on H6; it has not
been run yet. Twenty focused tests now pass across the changed modules.

## Completed extension and contraction experiments (supersedes live notes above)

All H6 experiments described here are terminal and their intervals have separate
standard-library replays. Values are Ha for the same frozen rational H.

| Run | Construction s | Solve wrapper s | Total s | Exact interval width |
|---|---:|---:|---:|---:|
| H4 diagonal, numerical contraction, Clarabel | .681 | 3.395 | 4.466 | .0000651119 |
| H4 diagonal, contraction, conditioned SCS30 | .663 | 30.188 | 31.084 | .000106215 |
| H6 full one-body driver, exact map, Clarabel120 | 466.428 | 171.678 | 656.378 | 6.63125 |
| H6 diagonal, exact map, Clarabel120 | 153.612 | 141.007 | 308.600 | .353744 |
| H6 diagonal, contraction, SCS180 | 14.516 | 235.500 | 261.977 | 18.1236 |
| H6 diagonal, contraction, conditioned SCS180 | 11.437 | 230.709 | 253.275 | 4.24923 |

These H6 failures are finite solver/export outcomes, not exact cone obstructions.
The existing depth-zero dual obstruction does not apply to these enriched cones.
No additional unmotivated long solve has been started for this dictionary.

`polynomial_gram_contraction.py` computes the exact monomial CAR map once and
then contracts numerical polynomial coefficients into it in batches. The final
SOS exporter remains rational and every physical coefficient is replayed. No
small coefficient is thresholded out. The backend reduces H6 symbolic word-pair
products from 3190092 to 200340, retaining 43740 Gram entries and 4839648 map
nonzeros. Construction dropped from 153.6 to 14.5 seconds in one comparison;
this is a construction improvement, not an accuracy or complete-work speedup.
Both outer products, their sum, sparse-dense products and peak primary temporary
array bytes are counted in receipts. Primary-array bytes are not total RSS.

A whole-map H4 comparison checked 3365964 entries against exact rational Gram
construction; the largest numerical discrepancy was 1.776e-15. H4 preserved
its independently verified accuracy with this backend. A separate tiny test
checks every map coefficient, including off-diagonal and degree-six terms.

Optional row conditioning divides each equality row by max(1, its map norm).
It leaves the unscaled residual objective unchanged and transforms saved dual
values back to physical coefficient units. Its integration control passes H4;
it did not make H6 accurate. Singular-value diagnostics of the original H4/H6
polynomial coefficient bases find full numerical row rank at tolerance 1e-10,
with estimated condition numbers 124–416 (H4) and 145–350 (H6). This does not
establish exact rank or a causal explanation of the H6 convergence failure.

Artifacts: `results/certificate_scaling/commutator_krylov/`, including each
run's `pre_solve.json`, `receipt.json`, `certificate.json`, `interval.json`,
`map_equivalence_h4.json` and `basis_conditioning.json`. Every discovery/map/
export cost remains charged. The focused campaign suite now has 32 passing
tests, including the separate dressed-quadratic branch; no full regression
claim is made.
