# Exact transfer under coherent-amplitude contraction

A checked H6 weighted-row certificate can now transfer across an entire hopping interval without solving a new positivity problem. The new verifier also handles targets whose spectator-dependent amplitudes change sign, a case the direct polynomial compiler refuses. Its conclusion is a complement lower bound; accurate ground-energy intervals still require separate witnesses and response calculations.

## The sufficient theorem

Let a reference certificate establish K_ref(s)≥0 with positive metric v on the full valence complement Q. Suppose a target Hamiltonian has the same diagonal, the same sector, and

\[
|H^{target}_{st}|\le |H^{ref}_{st}|\quad(s\ne t).
\]

Then its weighted Q-row numerator is at least the reference numerator for every s∈Q. The reference complement threshold remains valid. This argument specifically requires a weighted-row reference; the implementation rejects arbitrary PSD or unknown reference families.

No target sign assumption is needed. For a coherently grouped reference amplitude A with checked enclosure [lo,hi] and a constant shift Δ, magnitude contraction is equivalent to

\[
(A+\Delta)^2-A^2=\Delta(2A+\Delta)\le0.
\]

The maximum of this expression over the enclosure occurs at hi when Δ>0 and at lo when Δ<0. The verifier checks that rational endpoint expression exactly. It reconstructs all changed grouped amplitudes, refuses nonconstant shifts, checks the diagonal, verifies the reference proof, and binds the result to the actual target Hamiltonian and requested threshold.

Every convex interpolation H_θ=(1−θ)H_ref+θH_target, 0≤θ≤1, inherits magnitude contraction by convexity of absolute value. Thus the accepted endpoint certificate proves the whole segment. This is an analytic interval theorem; no parameter grid or sampling supplies the guarantee.

## The H6 interval

The tested perturbation adds δ times the spin-independent hopping between the first two localized orbitals. The four affected directed spin amplitudes all have the exact reference enclosure

\[
[-4274038918973/25000000000000,
 -14275927873493/100000000000000].
\]

Consequently the conservative contraction test permits

\[
0\le\delta\le14275927873493/50000000000000
=0.28551855746986.
\]

The reference threshold is −783/125=−6.264 Ha. This interval endpoint is sufficient; no maximality or optimality is claimed. The preserved metric and reference numerator were constructed previously. Reusing their proof does not remove that original discovery cost.

Independent replay artifacts cover δ=1/50, δ=3/50, δ=3/20 and the rational interval endpoint. At δ=3/20, the affected amplitude enclosures cross zero and the existing direct fixed-sign compiler rejects the target. The contraction verifier accepts it by comparing absolute magnitudes instead. The endpoint certificate also covers every intermediate δ analytically.

Only the δ=1/50 result is integrated with an existing verified ground-energy witness and response basis. Its energy interval remains the previously certified interval. The larger perturbations are Q-gap results, without newly constructed ground-energy intervals.

## Integration and validation

The energy-gap dispatcher admits a dedicated `joint_polynomial_monotone_transfer_v1` recipe containing only its kind and explicit reference certificate. It injects the current outer Hamiltonian and threshold, refuses hidden overrides, and requires the complete valence retained space. Reference sector mismatch, unverified reference kinds, increasing hopping magnitudes, changed diagonal terms and unsupported spectator-dependent shifts are refused.

Focused tests cover contractions, zero hopping, sign reversal, a target with genuinely mixed-sign spectator amplitudes, actual-H and threshold corruptions, a new edge, corrupt reference positivity and incomplete valence coverage. The mixed-sign test independently confirms that direct compilation refuses the target before exact transfer succeeds.

Source: `experiments/marginal_monotone_transfer.py`. Test: `tests/test_marginal_monotone_transfer.py`. Artifacts: `results/marginal_h6/polynomial_metric/bounded_quotient/monotone_*/` and `monotone_family.json`.

The remaining target is transfer under changes that increase relevant couplings or alter diagonal correlations, while discovering suitable metrics and response spaces with controlled cost. [The bounded-degree quotient construction](marginal_bounded_quotient.md) addresses a different part of that pipeline.

Final validation: **360 full marginal tests pass**. The current progress record accounts for all nine new independent standard-library replays.
