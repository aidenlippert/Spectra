# Stronger compact certificates from a free range-two profile

Allowing one additional scalar profile parameter produced stronger exact
certificates on the original nearest-neighbor model and both range-two
density targets. The million-site open intervals are:

| W | Lower/site, approximately | Recomputed physical upper/site, approximately | Gain over the preceding lower |
|---|---:|---:|---:|
| 0 | -0.6434299980096883 | -0.6106763470511881 | 9.548622846432e-5 |
| +1/10 | -0.6440671196088968 | -0.6114511605830021 | 1.3695253351456e-4 |
| -1/10 | -0.6430629666121880 | -0.6099015335193740 | 5.998238944816e-5 |

All three targets use U=4, t=1, V=1/2, the same projector sources and nine
sparse diagonal shapes. Each physical upper was evaluated afresh, including
the24-site exact contraction enclosure; it agrees with its previous value.
The improvement is entirely in the lower certificates.

## The compact consistency correction

The six-site range-two profile now has the form

    [5W/4 + delta, 5W/4 - delta, 5W/4 - delta, 5W/4 + delta].

Its sum remains5W, so changing delta preserves the translated Hamiltonian.
The added local term is

    delta (q0 q2 - q1 q3 - q2 q4 + q3 q5)
      = delta (Y_left - Y_right),
    Y = q0 q2 - q2 q4.

This identity was previously checked on every one of4096 six-site
determinants. The translated terms cancel exactly. At W=0 the entire
range-two profile is such a correction: the physical Hamiltonian gains no
range-two interaction. The winning W=0 value is delta=-39877/250000.

The five-site diagonal Y has320 nonzero entries, but its formula contains
only two charge products. The existing v6 profile field represents the
correction with four rational values. The64-entry sparse diagonal field,
all94 local Fock blocks and the maximum200-dimensional local PSD size remain
unchanged. This demonstrates that compact operator formulas can add useful
constraints that do not fit a sparse-entry representation.

The W=+1/10 profile is approximately
[-0.067653,0.317653,0.317653,-0.067653]. For W=-1/10 it is approximately
[-0.250332,0.000332,0.000332,-0.250332]. Signed local values are admissible
because their translated sum is checked exactly and local PSD is verified
on the full Fock space.

## Numerical construction and exact replay

Three bounded searches used74,75 and74 matrix evaluations, respectively,
under the same250-evaluation limit. The additional parameter was optimized
together with six nearest-profile parameters, two projector penalties and
nine sparse-shape coefficients. Fresh physical matrix reconstruction and
the affine model agreed within3.02e-14.

All proposals were rounded to exact rational certificates and replayed with
optional integer congruence witnesses. Complete lower-and-upper replays took
approximately29–30 seconds each while other bounded jobs were running.
These observations are not a controlled performance comparison. No GPU ran.

## Expanded-family ceilings

The old fixed-profile local PSD mixtures were not valid dual witnesses for
this larger family: they had nonzero expectations of the new profile
direction. A new accepting v2 family format requires that extra expectation
to vanish exactly. Its scalar constraint count increases from18 to19,
allowing at most19 physical mixture sources for nine shapes. The older
fixed-profile cap and all PSD dimension caps are unchanged.

Fresh19-vector rational mixtures certify these remaining periodic gaps:

| W | Gap to the expanded-family ceiling/site |
|---|---:|
| 0 | 9.794610292439077e-8 |
| +1/10 | 7.815154859454668e-8 |
| -1/10 | 1.1177187240864765e-7 |

All six nearest-profile derivatives, all nine sparse-shape moments and the
additional range-two profile moment are exactly zero. Trace, positivity and
both fidelity inequalities are exact. The caps cover every reflected,
mean-correct range-two profile for the fixed physical coupling, sources,
ratio, ceilings and sparse-shape span. They are limits on attainable LOWER
certificates, not physical energy upper bounds or universal relaxation limits.

## Remaining compact quadratic separator

An exact probe enumerated all six independent reflection-odd quadratic
charge polynomials on five sites, including charge squares. Five cancel
exactly, as implied by the now-enforced profile constraints. The remaining
one is

    Y = q0 q3 - q1 q4,
    Y_left - Y_right = q0 q3 - 2 q1 q4 + q2 q5.

Its expectations in the three accepted dual mixtures are approximately
-0.00478440,-0.00418119 and-0.00533049, respectively. Thus these local PSD
mixtures still fail necessary translation consistency. Y has416 nonzero
diagonal entries but only two charge products. Optimizing a certificate
with this additional correction remains untested; these moments are
separators, not claims of another energy improvement.

## Artifacts and scope

Outputs are under `results/marginal_graded_hubbard8/free_range_two_profile/`.
The root `combined_summary.json` compares all three cases. Each case directory
contains the energy certificate, optional PSD witnesses, `range_two_replay.json`,
`range_two_family_limit_certificate.json`, `range_two_family_limit_replay.json`
and `quadratic_charge_overlap.json`. Source/input hashes were checked after
acceptance. Earlier executed source versions are preserved in
`pre_free_profile_sources/`.

The focused tests include a formerly accepted fixed-profile dual that must
now fail, invariance under large signed profile changes, and the exact
one-source increase permitted by the added constraint. The full542-test regression and96 subtests passed in392.73 seconds;19
focused tests passed in13.85 seconds. Details are recorded in
`results/marginal_final_validation.json`.

These are stronger certificates and two range-two transfer examples within
one dimension. Arbitrary molecular or long-range interactions, higher
dimensions, general representability and requested-accuracy scalability
remain unproved. A fixed-family ceiling does not bound other source choices,
other correction spans or the physical ground energy from above.
