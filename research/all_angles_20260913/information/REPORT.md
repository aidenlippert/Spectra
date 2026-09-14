# Information / reduced-variable routes (45–52)

## Bounded attempt

`attempt_information_routes.py` exactly enumerates the two-particle sector of a four-site spinless fermion ring with nearest-neighbour hopping `t=1` and repulsion `U=2`. The dense diagonalization is an oracle only: it validates the diagnostic and is not used to construct a certificate. The run gives `E0=-1.2360679775` and a ground 1-RDM printed in `results/.../run.txt`.

The constructive test uses only diagonal occupation marginals. Two mixed ensembles, formed from translated adjacent-pair and alternating-pair configurations, both have occupation vector `(1/2,1/2,1/2,1/2)`, while their interaction energies differ. Thus a diagonal 1-RDM envelope cannot infer the interaction functional without extra two-body constraints. This is an adversarial marginal case, not evidence that full 1-RDM functionals are impossible.

For an inexpensive classical-plus-hopping lower bound, the diagonal interaction floor is exactly `0` (Fraction arithmetic); a row-sum hopping bound is exactly `4`, yielding `E >= -4`. It is valid but too loose by 2.7639 Ha on this control. The gap is a decisive falsifier for promoting this SCE-like floor plus operator norm into a useful molecular certificate without a sharper kinetic correction.

## Route status

* 45 quantum recovery / approximate Markov: **conditional/deferred**. Recovery requires a genuine compatible state and quantitative conditional mutual-information or gap assumptions; relaxed overlapping marginals are insufficient. No PSD recovery claim is made.
* 46 quantum belief propagation: **deferred prerequisite**. Needs a certified finite-temperature correlation-decay/gap theorem and a computable partition-function enclosure.
* 47 conditional factorization / entropy bounds: **tested diagnostically**. The marginal counterexample shows local occupation data alone do not control pair energy; additional overlap and entropy assumptions are required.
* 48 finite-temperature enclosures: **conditional**. For Hilbert dimension `D`, a certified free energy `F_beta` implies `F_beta <= E0 <= F_beta + log(D)/beta`; reaching 0.0016 Ha requires large beta and a hard partition-function enclosure. No such enclosure was claimed here.
* 49 certified 1-RDM envelopes: **tested diagnostically**; diagonal 1-RDM is provably too weak on the ring. Full 1-RDM plus a certified interaction envelope remains open.
* 50 SCE multi-marginal transport: **conditional/deferred**. Classical diagonal transport can provide a floor, but the hopping correction failed the tightness test above.
* 51 analytic Lieb–Oxford/Lieb–Thirring inequalities: **deferred prerequisite**. Constants and basis/sector translation must be certified before use.
* 52 generalized Pauli constraints/quasipinning: **conditional**. They constrain 1-RDM spectra but do not determine the two-body interaction; pairing with a certified functional is required.

## Novel hypothesis and falsifier

Hypothesis: on localized molecular fixtures, a *full* 1-RDM envelope augmented by a small number of certified pair-density channels could reduce certificate cost relative to a generic SOS dictionary. Decisive falsifier: on H4/H6, after optimizing the same channel count, the best valid lower interval remains wider than the current baseline at equal verification memory/time. The present ring result already falsifies the diagonal-only version.

## Refinement replay

Audit correction: an earlier draft incorrectly used `sqrt(n_i(1-n_j))` as a coherence bound; an asymmetric one-particle state disproves it. That line is superseded and must not be treated as a result. `refinement_checks.py` now uses the valid bounds `|<c_i^†c_j>| <= min(sqrt(n_i n_j),sqrt((1-n_i)(1-n_j)))`, with an asymmetric regression check (`n_i=.1,n_j=.9` allows `.3`). At `n_i=1/2`, the four-edge ring obtains the exact rational bound `E >= -4`; this is structurally stronger but numerically unchanged, so useful progress requires genuine 2-RDM positivity and overlap constraints.

The three-qubit check compares GHZ with the incoherent `|000>,|111>` mixture. Both are PSD and have identical classical `XZ` and `ZY` overlaps, but `XXX` expectations are `1` and `0`; this demonstrates nonuniqueness, not incompatibility. The sharper incompatibility is Bell_AB together with Bell_BC and common maximally mixed B: purity of Bell_AB forces C to factor, contradicting Bell_BC. An approximate-Markov hypothesis and number/parity preservation remain necessary.
