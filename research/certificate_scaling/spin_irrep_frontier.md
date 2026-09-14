# Exact SU(2) discovery backend and molecular scaling experiment

The spin-invariant discovery backend is implemented. Independent exact replay
passes chemical accuracy on H4 and H6; H8 SCS300 fails and H10 remains in progress. The previous
unrestricted chosen mixed-cubic H4/H6/H8/H10 ladder already passes .0016 Ha.
The question here is whether symmetry reduces total discovery cost at comparable
accuracy. Smaller PSD matrices alone do not answer that question.

## Exact basis construction and scope

`spin_basis.py` groups word spans by their spatial creation and annihilation
multisets. Degree at most three gives local spin patterns of dimension at most
eight. Exact rational raising-map kernels supply highest weights; divided-power
lowering supplies every descendant. Local rank checks prove completeness of
the original word span. No occupation-sector basis is enumerated.

`spin_parity.py` computes the full commuting GF(2) parity space directly on
spatial orbital bits, then lifts each bit to both spin modes. Simply filtering
noncommuting generators would miss valid combinations. All equivalent-copy
cross terms are retained in the invariant PSD multiplicity matrices.

For p_r=ad(S-)^r(p_0)/r!, invariant contractions use weights proportional to
1/binom(2j,r). Ranks j=0,1/2,1,3/2 use integer weights [1], [1,1], [2,1,2],
[3,1,1,3]. The exporter repeats rows by these weights and returns the existing
rational factor certificate format. Exact replay re-expands the ORIGINAL H.

`spin_invariant_discovery.py` uses the full chosen mixed-cubic word dictionary,
with its selected local pure triples; it is not every possible pure cubic word.
A dictionary that cuts a spin pair and is not spin-closed refuses rather than
silently expanding. The selected H4/H6/H8/H10 dictionaries are closed.

The frozen H fixtures have tiny SU(2) defects from rational rounding. Earlier
exact Casimir projection checked nearby invariant Hs with coefficient-l1
perturbations 2e-12, 5.4e-11, 2.36e-10, 4.74e-10 Ha respectively. The current
backend solves against the original H directly and charges all defects in the
full residual. No projected H is substituted in the final benchmark.

Symmetry averaging preserves the invariant zero-residual SOS cone. It does not
prove that the coordinatewise coefficient-l1 residual objective has the same
optimum as the unrestricted formulation; that norm is not SU(2)-invariant.

## Measured dimensions and completed solves

| Fixture | Previous Gram entries | Invariant Gram entries | Largest invariant PSD block |
|---|---:|---:|---:|
| H4 | 17,776 | 5,088 | 32 |
| H6 | 217,268 | 56,996 | 108 |
| H8 | 1,206,464 | 320,192 | 256 |
| H10 | 4,674,900 | 1,224,400 | 500 |

These are entry counts, not stored certificate bytes or asymptotic accuracy
bounds. Basis timing and corrected full ladder-action counts are in
`results/certificate_scaling/spin_irrep/decomposition_ladder_corrected_counts.json`.
The original count file undercounted validation actions; its dimensions were
correct. Full construction/solve wall times always included those actions.

| Run | Source seconds | Independent interval width (Ha) | Outcome |
|---|---:|---:|---|
| H4 SCS, 30 s solver cap | 6.298 | 7.898699e-7 | Pass, coefficient replay |
| H4 Clarabel, 60 s cap | 3.206 | .00169400542 | Fail, coefficient replay |
| Same H4 Clarabel + spectral residual proof | +.791 proof | .000375913598 | Pass, spectral replay |
| H6 SCS, 180 s cap | 209.352 | .000068110864 | Pass, coefficient replay |
| H6 Clarabel, 120 s cap | 185.020 | 7.189853963 | Fail, coefficient replay |

H4/H6 SCS factor rows are 703/1,581 and factor nonzeros 17,008/135,804;
certificate sizes are 176,459/1,034,974 bytes. Map nonzeros are 27,456/391,588.
The H4 previous monomial map had 10,256 nonzeros: spin compression increases
map fill. Independent replay takes .611/6.906 seconds respectively, beyond
source times. Solver statuses and caps are diagnostics; only exact intervals
count as passes. H6 source downloads hash-match all six remote files.

H8 SCS300 on instance A and H10 SCS600 on B use frozen archive
`spin_discovery_v2_source.tar.gz`, SHA256
`7b45026d6fe407fa00e9cef636fa317e1bc9a8658d249654d26d4ec5057fcf10`.
External whole-process caps are 900/1500 seconds. Version two adds early
sector/budget refusal and solver iteration/options metadata; numerical algebra
is unchanged from the completed H6 runs. Output root on both hosts is
`/home/ubuntu/spectra-spin-discovery-v2/results/spin_runs/`.

## Validation and remaining obligations

52 focused tests pass across the new spin modules, interacting parent control,
checkpoint/discovery, polynomial contraction and residual/certificate helpers.
This is not a complete repository regression. New tests include independent
spin commutators, all descendant ranks, cross-copy maps, original-H positive
controls with symmetry defects, and invalid input/export refusal.

The exact molecular upper witnesses remain exponential-cost FCI validation
artifacts; their discovery is not removed by this backend. A passing finite
ladder and the fixed-degree polynomial dictionary do not prove that required
accuracy persists at growing system size. Symmetry is a reduction mechanism;
a chemistry-relevant structural accuracy/discovery theorem remains open.

## Completed H8 boundary and subsequent map optimization

H8 SCS300 completed 2,500 iterations in 555.279 seconds total: construction
147.973, solve304.662, export102.644. Exact coefficient replay gives width
6.276704644 Ha. A separate spectral residual proof improves this to
0.699480442455 Ha, still FAIL. Proof generation costs90.081s and independent
spectral interval replay costs93.802s, in addition to the source run. The
negative Gram eigenvalue mass before clipping was.0876091. All costs and
failed outcomes remain in the campaign receipt.

H10 construction took681.690s with10,040,532 map nonzeros. The numerical
solver cap stays600s. Its outer ceiling was extended from1500 to2100s to
allow exact export after the measured construction overrun. A bounded guard
supervises the child while suspending only the earlier timeout timer and
resumes that timer when the child completes or reaches the new ceiling.
The original frozen v2 numerical code is unchanged during this run.

Profiling H6 showed repeated polynomial canonicalization and exact Fraction
arithmetic dominate column construction. An initial fused accumulator was
9% slower and was rejected. The adopted `spin_gram_columns.py` precomputes
channel adjoints and uses the exact permutation sign of already-canonical
monomials for adjoint maps. It retains exact Hermiticity checks. All28,894 H6
columns, indices and work counts match an independent CAR reference exactly.
With both caches cleared before each builder, construction took16.024s for
the reference and9.452s for the revised builder. This single local comparison
is a column-construction measurement, not a total solver speed claim.
Two dedicated tests also compare every H4 column and rational cross-copy
polynomials. The backend now calls this helper for future runs. Frozen v3
sourceSHA8058408841be83b73ad9383bbc9fe4f155de8be33723152e40aa0bff58f063ac
is saved locally and has not replaced the running cloud source.
