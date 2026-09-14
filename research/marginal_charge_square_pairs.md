# Charge-square pair certificates and a stronger-coupling transfer test

Four compact consistency corrections improve the exact energy lower bounds
on the three existing targets. The same representation also certifies a new
million-site target at W=1, ten times the previous positive range-two coupling.
All cases use U=4, t=1, V=1/2, an open half-filled chain, and N=1,000,000.

| W | Accepted lower/site, approximately | Physical upper/site, approximately | Gain over preceding lower |
|---|---:|---:|---:|
| 0 | -0.6432643847662889 | -0.6106763470511881 | 9.027242567648e-5 |
| +1/10 | -0.6439121977013573 | -0.6114511605830021 | 9.902535166736e-5 |
| -1/10 | -0.6428845043456355 | -0.6099015335193740 | 8.238359661880e-5 |
| +1 | -0.6617398444509426 | -0.6184244823693281 | New target |

Each lower certificate was replayed exactly on the complete six-site Fock
space. Each physical upper was recomputed, including an exact 24-site
contraction checked against its 160-bit enclosure before million-site
transfer. The three previous uppers are unchanged. No GPU ran.

## The four compact corrections

Define q_i=n_up,i+n_down,i-1 and p_i=q_i^2. The operator p_i is the indicator
of an empty or doubly occupied site. On five sites the four reflection-odd
pair operators are

    Y_ij = p_i p_j - p_(4-j) p_(4-i),
    (i,j) in {(0,1),(0,2),(0,3),(1,2)}.

The certificate adds a rational linear combination of Y_left-Y_right.
Each translated sum cancels exactly on the periodic chain, preserving the
physical Hamiltonian. The existing opening deduction is then applied to
the physical Hamiltonian. At W=0 no range-two or quartic interaction is
added to the target.

The v8 certificate uses a new bounded `charge_square_pair_telescope` field.
The four operators require eight charge-square products in total. Individually
the five-site diagonals have 384, 256, 384 and 256 nonzero entries; every local
telescope has norm 2. A separate exact rank replay over all 4096 determinants
proves that these four directions are independent of the six quadratic
charge telescopes, giving combined rank ten.

The separate sparse diagonal remains the same nine-shape span with 64
nonzero entries. All 94 local PSD blocks, their total dimension 4096, and
maximum local dimension 200 are unchanged. The four new coefficients were
optimized together with the previously free profiles, quadratic correction,
projector penalties and sparse coefficients. The four searches used
115, 115, 114 and 119 matrix evaluations, respectively, below the existing
250-evaluation cap. Fresh CAR reconstruction and the numerical affine model
agreed within 9.06e-14. Numerical results were proposals only; exact PSD
acceptance used freshly reconstructed integer matrices and optional
congruence witnesses.

## Exact ceilings for the expanded family

A new v4 family verifier requires all four charge-square pair moments to
vanish exactly, along with the previous quadratic and profile constraints.
The exact mixture-size cap increases from 20 to 24 for the nine-shape span;
older versions retain their earlier caps. All four actual dual witnesses
use 24 physical integer vectors with positive rational weights.

| W | Periodic gap from accepted lower to expanded-family ceiling/site |
|---|---:|
| 0 | 1.1348724240303176e-7 |
| +1/10 | 8.162929920975288e-8 |
| -1/10 | 9.313586930791479e-8 |
| +1 | 8.77759853651351e-8 |

These ceilings cover all four charge-square pair coefficients, all six
quadratic charge coefficients, all reflected mean-correct nearest and
range-two density profiles, and the existing sparse span. The projector
sources, ratio and overlap ceilings are fixed. Exact trace, positivity,
profile moments and both fidelity inequalities are checked. A previous
quadratic-only witness is explicitly rejected under the new pair constraints.

These are ceilings on attainable LOWER certificates for a fixed family.
They are not physical energy upper bounds, global representability proofs,
or limits on other sources, correction spans or window sizes.

## Transfer at W=1

The W=1 calculation retains the same projector sources, sparse shapes,
compact operators and physical trial-state recipe. Only scalar lower
parameters were retuned; the physical energy of the trial state was freshly
evaluated for the new Hamiltonian. The resulting physical interval has
width approximately 0.0433153621/site, compared with 0.0325880377/site at
W=0. Thus transfer succeeds at certification but does not preserve accuracy.

For context, the replay's historical nearest-model interval enlarged by the
operator-norm perturbation bound gives approximately
[-1.6435234842, 0.3893216529]. Direct certification closes about 97.87% of that
loose perturbative width. This comparison is against a norm bound, not an
exact solution or a competitive independent algorithm.

This is a stronger-coupling one-dimensional transfer example. It does not
establish transfer to arbitrary molecular interactions, arbitrary ranges,
higher dimensions or requested-accuracy scalability.

## Remaining consistency obstruction

The exact charge-polynomial probe again checks the complete 52-element
particle-hole-even, reflection-odd basis of functions of five site charges.
The six quadratic moments and four charge-square pair moments vanish;
the other 42 moments remain nonzero in every accepted dual mixture.
For the three original targets, the strongest normalized separator is

    Y = p2 p3 p4 - p0 p1 p2.

Its local telescope has norm 2 and its five-site diagonal has 192 nonzero
entries. Its expectations are approximately 0.00110700926, 0.00112470385 and
0.00108309356 for W=0,+1/10,-1/10. It tests consistency of triples of empty
or doubly occupied sites. No certificate using this triple correction has
been optimized or accepted yet.

The full reflection-odd indicator-only function space on five binary p_i
has 12 dimensions: two single-site directions, four pairs, four triples and
two quadruples. The single-site directions are already covered by the
quadratic profile constraints; the four pair directions are now covered.
Testing the remaining six indicator directions together is a concrete next
bounded experiment. Mixed signed-charge and noncommuting constraints would
still remain beyond that indicator-only span.

These are necessary translation-consistency tests on positive local
mixtures. A nonzero moment disproves the relevant consistency condition;
it does not by itself yield another energy improvement or a new physical
positivity theorem.

## Artifacts and verification

The root is `results/marginal_graded_hubbard8/charge_square_pairs/`.
`combined_summary.json` compares all four targets. Each target directory
contains the v8 energy certificate and optional congruence witnesses,
`range_two_replay.json`, the accepting v4 family certificate and replay,
and `charge_polynomial_overlap.json`. All 85 recorded source/input hashes
per target were checked. `charge_square_pair_rank.json` records the exact
rank calculation, with 28 source hashes checked.
`indicator_basis_summary.json` selects the six remaining indicator directions
from the accepted charge-polynomial probe receipts. Earlier executed sources are saved in
`pre_square_pair_sources/manifest.json`.

Focused tests passed: 22 energy/correction tests in 6.34 seconds and nine
family tests in 11.43 seconds. They cover independent full-Fock expansion,
exact diagonal PSD acceptance and rejection, particle-hole and reflection
symmetry, periodic cancellation, malformed inputs, rejection of the old
witness under the new constraint, and retention of the old mixture caps.
The full regression passed: 570 tests and 96 subtests in 416.07 seconds.
The result is recorded in `results/marginal_final_validation.json`.

The general representability and requested-accuracy scalability questions
remain open. This checkpoint supplies stronger finite-family certificates,
exact limits on further scalar tuning within that family, a farther transfer
example and a new consistency obstruction.
