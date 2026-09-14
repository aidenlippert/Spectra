# Exact limit of the nine-correction relaxation

This records the earlier fixed support. Newer accepted support replacements
and transfer results are in `research/marginal_shape_selection.md`.

The accepted nine-correction periodic lower certificate is
`-0.6436570569233962` per site. An independently replayed local PSD mixture
now caps every lower certificate in that fixed family at approximately
`-0.6436569911639332`. The exact remaining gap is
`6.575946298619209e-8` per periodic site. Thus further tuning within these
sources, ratio, ceilings and nine shapes can produce at most that gain.

This is a bound on the relaxation's best LOWER certificate. It is not a
physical ground-energy upper bound. The million-site open physical interval
remains approximately `[-0.6436595569233962, -0.6106763470511881]`.

## Exact dual argument

The witness is a convex mixture of18 physical six-site integer vectors.
Each has at most400 determinant amplitudes. Positive rational weights and
unit trace make the local density matrix PSD without a numerical eigenvalue
test. A standard-library replay constructs the physical CAR actions afresh
and checks all of the following exactly:

* All six reflected profile derivatives have zero expectation.
* Every one of the nine supplied diagonal differences Y_left-Y_right has
  zero expectation.
* Half-projector expectation is at most its fixed ceiling theta_half.
* Expectation of Phalf + (1/2) Pcharged is at most its fixed joint ceiling.

For nonnegative alpha and beta, taking the expectation of the local
positivity inequality then gives

    ell <= E[K0] + alpha E[Phalf] + beta E[Phalf + r Pcharged],
    (ell - alpha theta_half - beta theta_joint)/5 <= E[K0]/5.

All real coefficients of the nine shapes and all reflected mean-correct
profiles are covered by the zero moments. The source projectors, ratio and
ceilings remain fixed. The driver checks those values against the accepted
energy certificate and reconstructs its diagonal correction as an exact
linear combination of the nine shapes. The union contains60 entries.

The new bounded verifier allows at most9+k mixture sources for k<=9 shapes:
one trace row, six profile rows, k diagonal rows and two fidelity rows.
The proposal needed18 positive source weights. The older v1 verifier's
16-source cap is unchanged. No PSD matrix-size cap was raised.

## Construction and replay

Numerical discovery generated244 candidate states from eight active
reflection/spin sectors, solved a small linear program and reconstructed
its18-column basis with rational arithmetic. That proposal alone was
explicitly nonaccepting. Its discovery time was12.830 seconds.

`experiments/marginal_diagonal_family_limit.py` independently accepted the
physical mixture. The driver
`results/marginal_graded_hubbard8/discovery/diagonal_family_limit_replay.py`
replayed it and matched the energy family in1.891 seconds, using only the
standard library. No numerical solver is imported on this accepting path.

Artifacts are in
`results/marginal_graded_hubbard8/joint_projector/signed_density/symmetric_diagonals_9/extra_31/`:

* `diagonal_family_limit_proposal.json`: untrusted discovery output.
* `diagonal_family_limit_certificate.json`: exact input to the new verifier.
* `diagonal_family_limit_replay.json`: accepted limit, rational gap and hashes.
* `remaining_diagonal_overlap_moments.json`: exact remaining separators.

Focused old/new family tests passed, including refusal of the formerly
accepted dual when the violated diagonal condition is added. Its failure
is a mathematical distinction between the two relaxations, not a numerical
optimizer status. Validation details are in `results/marginal_final_validation.json`.

## What still fails physically

The new mixture still violates93 of the120 available symmetry-invariant
diagonal consistency directions. All nine selected moments cancel exactly.
The largest remaining direction is basis index11, with exact nonzero moment
approximately0.0039845523827132835. Its eight-entry shape gives a necessary
translation-consistency separator; it is not yet a stronger energy bound.

The current energy certificate uses60 of the64 allowed entries. Adding every
missing shape would exceed the accepted format. Possible next experiments
are exchanging existing shapes or combining them into a more effective
bounded-support correction; neither has been tested here. Other projector
sources or ratios could also escape this fixed-family limit.

One-dimensional negative-V parameter transfer has an earlier exact success.
Generic molecular or long-range transfer, higher dimensions, global
representability and end-to-end requested-accuracy scalability remain
unproved. The finite successes do not solve general quantum chemistry.
