# Hidden orbital structure: recovery and a molecular obstruction

The active chemistry-scaling goal remains open. This campaign implements a
Hamiltonian-only hidden-basis detector, validates exact certificates on a
known solvable interacting family, and derives a quantitative test that
excludes that particular approximation route on the molecular ladder.

## Recovery from the Hamiltonian

`hidden_density_basis.py` constructs the quartic commutator map
`A -> [Q(A),V]` on real symmetric one-body matrices. Numerical SVD and a
generic linear combination of kernel generators propose a simultaneous
eigenbasis. Acceptance requires exact rational orthogonality, exact recovery
of the density interaction and path, and an exact ratio-family match.
An exact finite-field rank minor also provides a rational nullity upper bound.
No generating rotation, state, or SOS factors are passed to discovery.

The controls use the previously verified biased fermionic ratio chain hidden
by two dense rational Householder rotations. The new transport certificate
contains the observed Hamiltonian, recovered rotation, and inner CAR
certificate. Its standard-library replayer checks orthogonality and the
Hamiltonian transformation through the exterior-square representation, then
replays the inner lower certificate and fixed-number structured upper proof.
The orbital unitary preserves particle number and spectrum. No many-body
matrix or state vector is enumerated.

| Modes | Inner factor rows | Certificate bytes | Discovery, seconds | Transport replay, seconds |
|---:|---:|---:|---:|---:|
| 4 | 3 | 4,562 | 0.316 | 0.002 |
| 6 | 5 | 19,244 | 0.071 | 0.006 |
| 8 | 7 | 60,955 | 0.285 | 0.019 |
| 10 | 9 | 146,850 | 1.888 | 0.119 |
| 12 | 11 | 321,277 | 4.854 | 0.272 |
| 14 | 13 | 587,345 | 11.032 | 0.576 |
| 16 | 15 | 1,053,465 | 22.973 | 1.114 |

All seven intervals are exactly `[0,0]`. Sizes 4–8 ran locally; 10–16 ran
on the warm A10 hosts, using CPU NumPy/SciPy. These are single runs, with
first-use import overhead included in the M4 discovery time. They are not
cross-machine speed comparisons or a GPU benchmark. Fixture construction is
recorded separately and excluded from discovery.

The observed dense Hamiltonian already has order M^4 coefficients; the
rotation has M^2 entries and the inner factors have linear support. The map
has polynomially many entries, but its dense linear algebra is costly.
Numerical rational reconstruction is bounded and may refuse; no completeness
or uniform polynomial-time discovery guarantee is proved for arbitrary
coefficient heights or small joint eigenvalue gaps. This is recognition of a
known solvable family, not a new solution of general interacting matter.

## The molecular condition fails, with an exact witness

For the frozen canonical linear H4/H6/H8/H10 fixtures, exact modular rank
minors prove real-symmetric commutant nullity at most 3. The required dimensions
would be 8, 12, 16, and 20 respectively. Thus there is no exact real orbital
density basis for those quartic interactions. This is stronger than one
failed search for a rotation.

The quantitative obstruction also covers complex orbital unitaries and
arbitrary one-body-multiplier number ideals. Let P remove the quartic
one-body lift from the interaction. On the full Hermitian one-body domain,
construct C(A)=[Q(A),P(V)]. If its M smallest squared singular values are
lambda_1,...,lambda_M, then every decomposition

    V = U V_density U† + L(X) + R

satisfies

    ||R||_F² >= (lambda_1 + ... + lambda_M)/4.

The Frobenius norm here is the quartic coefficient matrix on the two-particle
wedge. It is invariant under orbital unitaries and no larger than the
coefficient L1 norm. This is not the norm of the fixed-N many-body operator.

The proof uses the M orthonormal occupation projectors, the fact that each
quartic term changes at most four occupations, and the minimum-trace
principle. Full details, including the diagonal metric, imaginary sector and
number-ideal projection, are in `approximate_commutant_obstruction.md`.

`commutant_obstruction.py` replaces floating eigenvalue claims with exact
integer Gram matrices and fraction-free inertia. At threshold t=1/10000,
if k eigenvalues lie below t, it proves the squared distance is at least
`max(0,M-k)*t/4`. Every pivot division and sign is exact; zero leading pivots
are refused. The norm is rounded downward to a rational number.

| Chain | Modes | Negative inertia, real + imaginary | Certified distance floor, Ha |
|---|---:|---:|---:|
| H4 | 8 | 3 + 1 | 0.010000000000 |
| H6 | 12 | 3 + 1 | 0.014142135623 |
| H8 | 16 | 3 + 1 | 0.017320508075 |
| H10 | 20 | 6 + 2 | 0.017320508075 |

These are the stronger **number-ideal-quotiented, complex-unitary** results.
Without quotienting, H10's floor is 0.02 Ha; that number must not be reused
for the quotient case. All four floors exceed 0.0016 Ha. Therefore replacing
the interaction by a rotated density-density interaction and paying a small
coefficient residual cannot meet that residual budget, even after optimizing
over one-body number-ideal multipliers.

This does not rule out accurate energy predictions, higher-body number-ideal
multipliers, alternate residual norms, or general noncommuting SOS factors.
It is a structural exclusion for frozen rational Hamiltonians and the stated
representation. A lower bound on coefficient residual is not a lower bound
on ground-energy error.

## Verification and provenance

`hidden_campaign_replay.py`, run with `python -S`, recomputes the seven
transport proofs, four exact rank witnesses, and four quotient obstructions.
It compares the full inertia records and pivot hashes, checks fixture hashes,
and writes `results/certificate_scaling/hidden_density_basis/validation.json`.
All 15 proof replays passed in 107.64 seconds. All 16 focused tests passed.
Focused tests cover independent CAR commutators, genuine mixed orbital
rotations, transport rejection, exact congruences with known inertia,
projection idempotence and covariance, and known approximation bounds.

The initial eight remote discovery/detection jobs all exited successfully.
All 22 downloaded files matched their remote SHA-256 manifests. The executed
initial snapshot is `results/lambda_runs/hidden_density_basis/source.tar.gz`.
Later proof code is frozen separately in `final_source.tar.gz`. One preliminary
H10 real-obstruction run failed while converting huge pivots to decimal for
hashing; hexadecimal hashing fixed that portability issue. It was rerun, and
the successful complex and quotient results are the accepted evidence.

Both existing A10 hosts remain warm at the previously verified combined
rate of $2.58/hour. No instances were launched or terminated in this campaign.

## Next attack

The molecular certificates must retain material noncommuting interaction
structure. The next practical experiment is to improve discovery and exact
export for the noncommuting cubic H6 certificate: its previous numerical
objective was within the target, but exact residual correction destroyed
that accuracy. Test whether better cone conditioning
and solve precision produce a genuinely accurate exact certificate, then
study its factor structure and test transfer to H8. A numerical objective
alone is not an accepted lower bound. This is a prerequisite experiment for
finding compressible noncommuting factors, not itself a scaling theorem.
