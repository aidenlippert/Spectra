# Symbolic extraction: a number-ideal identity replaces the sector check

## Result

Four successful hopping certificates have been converted into explicit CAR
polynomial identities of the form

    H-bI = sum_alpha B_alpha† B_alpha + (Nhat-N) X + R.

The new lower-bound verifier uses exact rational polynomial arithmetic and does
not construct a Fock basis, a fixed-number sector matrix, or eigenvalues.
Its soundness argument applies to arbitrary finite mode count M and declared
particle number N. The nontrivial extracted examples remain six-mode examples.

All four original identities and one symmetry-compressed variant replay. The
matched certificate's uncorrected polynomial residual has coefficient l1 norm
approximately 11.93; its explicit degree-four number multiplier reduces that to
0.000064370835. This identifies the number-sector relation that the previous
matrix check had used implicitly.

## The algebraic bridge

Use normal-ordered, number-conserving ladder monomials. If X is Hermitian and
number conserving, it commutes with Nhat, and (Nhat-N)X vanishes on the target
sector. It may have large coefficients and need not be positive.

Every B_alpha is supplied by a rational Gram factor, so B_alpha†B_alpha is
positive by construction. The checker expands the identity exactly using
canonical anticommutation relations. For the resulting Hermitian residual

    R = sum_w r_w W_w,

each W_w is a product of norm-at-most-one ladder operators. Consequently

    ||R|| <= eta := sum_w |r_w|,
    H >= (b-eta) I on the N-particle sector.

This is a direct operator inequality. No assumption that the proposed multiplier
is optimal, that the numerical solve converged, or that the proof system is
complete enters its validity. Failure to fit the ideal produces a larger R and
a weaker lower bound rather than an unsupported certificate.

The relevant representation is therefore **positive fermionic squares modulo
the fixed-number ideal**, with an explicitly bounded remainder. This supplies
the symbolic version of the previously implicit finite-sector equality.

## How X was extracted

The imported factors are the unchanged rational factors from the prior hopping
experiment. Their exact Gram products are accumulated and CAR-normal-ordered.
For M=6, the proposer searches X in the 142-dimensional real Hermitian space of
body rank at most two: 1 scalar, 21 Hermitian one-body generators, and 120
Hermitian two-body generators.

Multiplication by Nhat-3 maps those generators into body rank at most three.
The numerical proposer solves a coefficient-space least-squares system with up
to 662 canonical monomial rows and rounds coefficients to denominator 10^10.
The accepted result comes only from exact rational re-expansion, never from the
least-squares residual or reported rank. No universal completeness theorem for
this multiplier ansatz is assumed.

The largest polynomial space has size

    sum_(k=0..3) binomial(M,k)^2 = O(M^6)

at this fixed body degree. That removes exponential sector enumeration from
this verifier, but does not make sixth-power storage inexpensive or prove that
fixed degree will suffice for other Hamiltonians. Certificate coefficient bit
length, factor count, sparsity, and discovery cost remain separate issues.

## Verified bounds

The upper endpoints below are the exact rational variational upper bounds
already verified in the previous experiment. The new procedure replaces the
lower-bound check only; it does not replace that upper-state computation.

| Case | t | Symbolic lower | Previously verified upper | Interval width |
|---|---:|---:|---:|---:|
| Matched, full mixed factors | 0.2 | 0.550936271779 | 0.551000400320345 | 0.0000641285 |
| Asymmetric, full mixed factors | 0.2 | 0.540804119066 | 0.540877181133642 | 0.0000730621 |
| Matched, eight selected pairs | 1 | -1.645960112550 | -1.645751311064574 | 0.0002088015 |
| Asymmetric, eight selected pairs | 1 | -1.692917759912 | -1.692748540472140 | 0.0001692194 |

These intervals are wider than the earlier sector-row-norm intervals because
the dimension-independent monomial coefficient bound is more conservative.
They still retain useful accuracy without sector enumeration in the lower
verifier. Units are those of the model's unit interaction, not hartrees.

## Symmetry compression of the multiplier

For the matched model, average X over simultaneous permutations of the three
flavors and exchange of the two triples (12 orbital permutations). Discard
averaged coefficients smaller than 10^-6 in magnitude. This is a proposal:
the checker recomputes the full residual after the operation, charging every
discarded or changed term.

The result has 64 monomials organized into ten signed permutation orbits.
The certificate stores ten representative words with rational coefficients and
the orbital permutations. The verifier reconstructs the signed orbit sums,
checks Hermiticity and number conservation, and replays the identity. The ten
orbits have sizes 1,6,6,6,12,6,3,12,6,6; adjoint-related orbits share a
coefficient in this example.

This compressed multiplier gives

    lower = 0.550912917080333...,
    eta   = 0.000087725533666...,
    width = 0.000087483240012... .

Only the multiplier has been compressed this way. The full matched certificate
still has 242 positive-factor rows and 7,124 nonzero factor coefficients.
Those factors are not replaced by the ten multiplier orbits. The asymmetric
case and the selected-pair cases have different coefficient counts, recorded
in the machine-readable receipts.

## Verification and independent checks

- The exact verifier replays five certificates with standard-library Python.
  Running Python with `-S` disables site packages; neither NumPy, SciPy nor
  CVXPY is needed for verification or orbit expansion.
- Seven new tests cover the exact number identity, the exact earlier triple
  identity, incorrect particle number, scalar/factor corruption, invalid
  coefficients, non-Hermitian or nonconserving multipliers, saved hopping
  replay, and orbit corruption.
- A simple number identity is checked at M=1000,N=500, as well as empty and
  filled edge sectors. This checks the generic representation and absence of
  enumeration; it is not a 1000-mode interacting-material calculation.
- Independent occupation-state action checks reproduce the four imported
  Hamiltonian matrices. Those are tests, not part of the symbolic checker.
- The existing six CAR tests and seven hopping tests remain part of the
  focused verification. Together with the new tests they total 20.
- An independent agent reviewed the exact verifier and replayed both the
  uncompressed and compressed matched certificate without finding a material
  soundness error. This is code review, not formal proof-assistant verification.

The verifier is restricted to real rational coefficients and homogeneous-charge
operator dictionaries of ladder degree at most three. H and X have degree at
most four. It rejects malformed sectors, invalid operator indices, non-rational
encodings and ambiguous multiplier representations. A well-formed poor proof is
reported with its honest residual and lower bound.

## What remains unsolved

The symbolic bridge is implemented. The certificate search still used
six-mode sector SDPs, and the upper state still came from numerical exact
diagonalization followed by an exact rational Rayleigh quotient. Neither
discovery nor the complete primal/dual workflow has been made scalable.

The next concrete target is to search directly in coefficient space for

    H-bI = w†Qw + (Nhat-N)X + R,

using small independently weighted operator dictionaries. That would move
discovery onto the same algebra already used by this verifier. Growth of the
necessary dictionary and multiplier degree, rather than proof-checking by
sector enumeration, would then be the principal representation question.

## Reproduction

From the Spectra root, extraction uses the existing NumPy environment:

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m experiments.marginal_symbolic

Verification and compression use standard-library Python:

    python3 -S -m experiments.marginal_symbolic --verify results/marginal_symbolic/full_mixed_matched.json
    python3 -S -m experiments.marginal_symbolic --compress results/marginal_symbolic/full_mixed_matched.json
    python3 -S -m unittest tests.test_marginal_hunt_car tests.test_marginal_symbolic -v

The complete focused suite additionally uses the prior numerical environment:

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m unittest tests.test_marginal_hunt_car tests.test_marginal_hopping tests.test_marginal_symbolic -v

[Implementation](../experiments/marginal_symbolic.py),
[new tests](../tests/test_marginal_symbolic.py),
[extraction receipts](../results/marginal_symbolic/summary.json),
[symbolic replay](../results/marginal_symbolic/replay.json),
[compressed certificate](../results/marginal_symbolic/matched_symmetry_compressed.json),
[orbit representatives](../results/marginal_symbolic/matched_symmetry_orbits.json),
[independent review](marginal_symbolic_review.md).
