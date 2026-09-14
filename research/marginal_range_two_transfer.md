# Exact transfer beyond nearest-neighbor interactions

The construction now has accepted million-site intervals for

    H = U sum_i n_i,up n_i,down
        - t sum_i,spin (c†_i,spin c_i+1,spin + h.c.)
        + V sum_i q_i q_i+1 + W sum_i q_i q_i+2,
    q_i = n_i,up + n_i,down - 1,

on the open, half-filled one-dimensional chain. Here U=4, t=1 and V=1/2.
The additional W interaction extends the Hamiltonian beyond nearest neighbors.

| W | Certified lower/site, approximately | Physical upper/site, approximately | Interval width |
|---|---:|---:|---:|
| +1/10 | -0.6442040721424114 | -0.6114511605830021 | 0.03275291155940926 |
| -1/10 | -0.6431229490016361 | -0.6099015335193740 | 0.03322141548226207 |

The original projector sources, nine diagonal shapes, block trial state and
boundary filter recipe were transferred. Scalar lower profiles, penalties
and shape coefficients were retuned. No new projector source or trial-state
optimization was used. Both energy intervals passed a fresh standard-library
replay in approximately22.2 seconds, including the physical upper calculation.

## Lower-bound accounting

The v6 certificate retains all94 six-site Fock/reflection blocks and the
existing overlapping-projector inequalities. The new term is diagonal:

    sum_{j=0}^3 w_j q_j q_j+2,
    w_j = w_3-j,  sum_j w_j = 5 W.

The numerical searches used the fixed uniform profile w_j=5W/4. The lower
verifier also accepts explicitly supplied signed, reflected profiles with
the same exact sum. Summing all translated six-site windows and dividing
by five reproduces the physical W coupling. The centered local onsite
shift vanishes on the global half-filled sector, as in the earlier proof.

Opening the periodic chain removes one nearest-neighbor bond and two
next-nearest density bonds. The norm deduction per chain is

    2t + |V| + 2|W|.

The added density factors have norm one because each q_i has eigenvalues
-1,0,1. Exact tests enumerate every eight-site determinant to check the
window counting and the two removed density bonds. Separate all-Fock tests
compare the local lower against an independently calculated diagonal minimum,
including signed local coefficients and refusals above the true minimum.

The64-entry telescoping format and maximum200-dimensional local PSD blocks
are unchanged. Range-two fields are refused in the old energy version.

## Fresh physical upper

Each next-nearest density term touches at most one boundary filter when
blocks have at least four sites. A dressed term spans at most two adjacent
blocks, so the existing48-dimensional norm/energy transfer recurrence can
be reused. The new compiler evaluates six-site cut patches containing the
pairs(0,2),(1,3),(2,4),(3,5), with the filter on sites2 and3. Four-site partial
traces of the physical block supply these patches; supplied RDMs are not
accepted as physical input.

The upper implementation was tested against independently expanded physical
filtered states on4,8 and12 sites, comparing both exact norm and exact energy.
These tests include positive and negative W and exact outward enclosures.
For the actual H8 block recipe, the replay checks a24-site exact contraction
inside its enclosure, then evaluates the million-site chain at160 bits.

For comparison, enlarging the old accepted interval by the elementary norm
bound |W|(N-2)/N on each side gives width about0.232849. The new certified
intervals close85.93% and85.73% of that conservative width. This comparison
measures improvement over a norm-perturbation estimate, not a scaling law.

## Limits of the two fixed-profile relaxations

New18-vector local PSD mixtures certify the following gaps between each
accepted periodic lower and its fixed-family ceiling:

| W | Remaining fixed-profile family gap/site |
|---|---:|
| +1/10 | 9.446437040882777e-8 |
| -1/10 | 5.122156012051345e-8 |

The exact verifier reuses all six nearest-profile cancellations, nine
diagonal cancellations, normalization and fidelity inequalities. It then
adds the exact expectation of the specified range-two local profile.
These caps fix that profile; they do not cap optimization over every
mean-correct range-two profile. The old nearest-only family verifier and
comparison drivers refuse range-two targets.

The distinction is consequential. The two accepted dual mixtures have
nonzero expectations of the profile direction[1,-1,-1,1]: approximately
-0.00705543 and-0.00476930. An exact4096-state identity identifies this as

    Y(first five sites) - Y(last five sites),
    Y = q0 q2 - q2 q4.

Thus the direction sums to zero under translation and is a necessary
consistency constraint. Its five-site diagonal has320 nonzero entries, but
its formula has only two charge products. The existing v6 profile field can
represent this change without expanding the64-entry sparse telescope.
Optimizing that free profile and certifying its larger-family limit remain
next steps; no resulting energy improvement is claimed here.

## Validation and artifacts

The534-test marginal regression passed in368.15 seconds. A subsequent
19-test family supplement passed in14.01 seconds; it covers the new
fixed-profile wrapper and the final explicit nearest-only refusal guard.
The accepting energy, family-limit and separator receipts have verified
source/input hashes. Modified earlier sources are preserved with hashes.

Production modules:

* `experiments/marginal_range_two_density.py`: exact local profile and charge diagonal.
* `experiments/marginal_projector_extendibility.py`: v6 lower acceptance.
* `experiments/marginal_range_two_transfer.py`: physical upper compilation.
* `experiments/marginal_range_two_family_limit.py`: fixed-profile family cap.

Artifacts are under `results/marginal_graded_hubbard8/range_two_density/`.
The root `combined_summary.json` compares both targets. Each W directory
contains `range_two_replay.json`, `range_two_family_limit_replay.json`,
`free_profile_separator.json`, the certificates and optional PSD witnesses.

This establishes two range-two density-transfer examples within one
dimension. Arbitrary long-range or molecular integrals, higher dimensions,
general representability and end-to-end requested-accuracy scalability
remain unproved. Fixed-width proof replay and a million-site instance do
not establish those broader claims. No GPU was launched for this work.
