# Matching certificates beyond the original Hubbard target

The boundary-transfer compiler now produces both upper and lower energy
certificates for three specified nearest-neighbor models,

\[
H(U,t,V)=-t\sum_{i,\sigma}(c^\dagger_{i\sigma}c_{i+1,\sigma}+h.c.)
+U\sum_iD_i+V\sum_i(n_i-1)(n_{i+1}-1).
\]

The following intervals apply to half-filled open chains of one million
sites. Displayed endpoints are rounded outward; exact rational endpoints
are stored in the replay receipt.

| U | t | V | Ground-energy lower/site | Ground-energy upper/site |
|---:|---:|---:|---:|---:|
| 4 | 1 | 1/2 | -0.66320650000 | -0.61067634705 |
| 4 | 1 | -1/2 | -0.56730850000 | -0.52407461381 |
| 3 | 2/3 | 1/2 | -0.42042183334 | -0.38877917564 |

The interval widths are approximately 0.052530152949, 0.043233886186,
and 0.031642657689 per site, respectively. This extends the earlier
upper-only target calculations with matching local positivity certificates.
It does not establish transfer to arbitrary interactions or geometries.

## What changes in the upper compiler

The physical H8 source state remains the original, independently validated
U4,t1 block. The accepting replay always reconstructs that source before
evaluating a new target. The gate parameters are retuned for each target:

| U,t,V | a | b |
|---|---:|---:|
| 4,1,1/2 | 0.374704 | 0.177582 |
| 4,1,-1/2 | 0.339851 | 0.179145 |
| 3,2/3,1/2 | 0.353530 | 0.174082 |

The gate is F=I-a*h+b*h^2, with h the unit-hopping contact operator. The
target hopping t multiplies the Hamiltonian, independently of this gate
parameterization. The energy compiler recomputes every onsite, hopping,
and density contribution for the specified U,t,V.

Density terms are assigned to the same bond partition as hopping: middle
bonds belong to a block insertion, the outer bonds belong to the boundary
insertions, and the three bonds around a filtered cut belong to its
four-site dressed patch. That patch contains F-dagger H_patch F. No scalar
remote-observable invariance is assumed.

The 16-dimensional norm transfer and 48-dimensional energy transfer are
unchanged in size. At 160-bit precision the three certified energy-density
enclosure widths are at most 1.265e-37. As before, these endpoints enclose
the trial-state energy; only the upper endpoint supplies a variational
ground-energy upper bound.

The resulting upper bounds improve the preceding linear-filter target
bounds by approximately 0.00983234, 0.00869068, and 0.00617535 per site.
Retuning a,b also improves on transferring the original quadratic gate
unchanged. Comparing the new upper endpoints with the unchanged gate's
lower trial-energy endpoints verifies that the trial expectations improve,
not merely their numerical upper estimates. Global optimality of these
parameters is unproved.

## Why the lower certificates match these targets

For a length-L window use the centered local operator

\[
K_L=\sum_{i=0}^{L-1}u_i(D_i-\tfrac12n_i+\tfrac12)
-\sum_{j=0}^{L-2}t_jT_j
+\sum_{j=0}^{L-2}v_j(n_j-1)(n_{j+1}-1),
\]

with exact profile sums

\[
\sum_i u_i=(L-1)U,\qquad
\sum_jt_j=(L-1)t,\qquad
\sum_jv_j=(L-1)V.
\]

Each replay checks K_L>=ell*I in every local Fock sector. For L=6 this
covers all 4,096 states, split into 94 reflected spin-number blocks with
maximum matrix dimension 200. Checking only a half-filled local sector
would not suffice: a window inside a half-filled chain can exchange charge
with its environment.

Summing translated windows around a periodic chain and dividing by L-1
gives

\[
\frac1{L-1}\sum_xK_L^{(x)}
=H_{\rm periodic}-\frac U2(\widehat N-N).
\]

The centering term vanishes at global half filling. Thus the periodic
energy density is bounded below by ell/(L-1). Removing the closing bond
changes the operator by norm at most 2t+abs(V), yielding

\[
\frac{E_0^{\rm open}}N\ge
\frac\ell{L-1}-\frac{2t+|V|}{N}.
\]

This argument applies to either sign of V. The verifier checks the target
profile means before accepting a pairing of upper and lower certificates.

The three six-site lower values are -3.31602, -2.83653, and -2.10210.
Their profile shapes are scaled from the earlier weighted window; these
certificates do not assert that those shapes are optimal for the new
targets. The density profile equals V times the old unit-mean hopping
profile, so its five entries sum to 5V exactly.

## Verification and failed proposals

Direct CAR tests compare norms and the full physical Hamiltonian energy
on one, two, and three four-site blocks. The source includes charge
fluctuations, and the independent reference computes the density energy
directly from occupation labels rather than using the new symbolic
density-term generator.

Additional tests check the density diagonal on all two-site Fock states,
profile symmetry refusals, and target-window mismatch refusals. A density-
only local model has an exact minimum of -1; attempting to certify -1/2
is rejected for either sign of the interaction. This ensures the new
interaction cannot silently disappear from the lower proof.

The first implementation exposed a hopping-parameter shadowing bug during
regression. It was fixed before acceptance. Numerical proposals from that
version are retained only as a rejected diagnostic, with their invalidity
explicitly recorded. A test-reference omission of the second site's onsite
term was also corrected. Neither failing run is counted as validation.

```sh
OPENBLAS_NUM_THREADS=1 python -S results/marginal_graded_hubbard8/discovery/density_transfer.py
```

The fresh three-target replay completed in 156.867 seconds. Each target
includes an exact 24-site comparison, a million-site comparison with the
independent earlier linear-filter implementation, a retuned quadratic
upper, and its own fresh all-Fock lower proof. All 12 source/input hashes
matched. The original U4,t1,V0 certificate also replayed unchanged after
the extension, in 55.575 seconds.

All 31 focused tests passed in 22.919 seconds. All 461 regression tests
passed in 335.097 seconds. Validation is recorded in
`results/marginal_final_validation.json`.

## An exact ceiling and an extension obstruction

For U4,t1,V1/2, the six-site weighting family now also has an exact
limitation certificate. A physical six-particle ring state, specified by
400 integer amplitudes, is exactly invariant under fermionic cyclic
translation. Its centered onsite, hopping and density expectations are
uniform. Therefore every six-site profile with sums 20,5,5/2 has the same
expectation on this state:

\[
\langle K_6\rangle=\frac56\langle H_{\rm ring,6}\rangle.
\]

An all-Fock inequality K6>=ell*I must also hold on this particular state.
Consequently the best possible periodic lower density from this profile
family obeys the exact, outward-rounded ceiling

\[
\sup_{\rm profiles}\frac\ell5\le -0.66302383.
\]

The witness's exact value is approximately -0.6630238313738540.
The current certified periodic lower density is -0.663204. Thus all
remaining profile tuning can recover at most about 0.000180168627 per
site, under 0.35% of the present interval. This ceiling imposes no
positivity or reflection restriction on the profile entries themselves;
it follows only from their sums and physical operator positivity.

The witness also has exactly equal five-site reductions on its first and
last five sites. Their purity is approximately 0.3130168036318084, strictly
less than one. The equality is checked entry by entry with integer
arithmetic, so every Hermitian even five-site boundary correction
A_left-A_right has zero expectation. The same ceiling therefore survives
arbitrary such telescoping corrections.

Nevertheless, the pure six-site marginal cannot repeat on two overlapping
windows of a seven-site state. A pure marginal on sites 1 through 6
forces factorization with site 7. Its marginal on sites 2 through 7
would then be rho_(2..6) tensor sigma_7, whose purity is strictly less
than one. It cannot equal the required pure six-site marginal. This is
an explicit local-consistency witness excluded by global extendibility.

```sh
python -S results/marginal_graded_hubbard8/discovery/density_window_ceiling.py
```

This replay validates the integer ring state, all six signed translations,
uniform observables, the energy identity, exact overlap equality, and
mixed-reduction purity. The finite-ring trial energy is used solely as
a ceiling on the local lower-certificate family; it is not treated as a
bulk ground-energy bound. Its four source/input hashes are recorded in
`density_transfer/window_ceiling_replay.json`.

An auxiliary profile-search diagnostic had omitted the reflection-orbit
Gram normalization. Its overstated improvement was rejected. An independent
Gram-normalized recomputation gives approximately -0.6631938638557365
for its rational proposal, a small gain over the inherited numerical
profile minimum. It is not an accepted new lower certificate. The exact
ring-state ceiling establishes the limitation without trusting that search.

The quantitative extendibility target has now been achieved for this
witness. Exact four- and five-projector inequalities cross the family
ceiling; the refined million-site lower bound is approximately
-0.646165607326505, with the same upper endpoint. The full derivation and
469-test validation are in [the six-site projector report](marginal_six_site_projector.md).
Remaining tasks include charged-sector constraints, improved block states
and contact operators, broader interaction ranges and geometries, and
accuracy-versus-cost guarantees. General marginal representability,
thermal and dynamic prediction, and materials synthesis remain unproved.
