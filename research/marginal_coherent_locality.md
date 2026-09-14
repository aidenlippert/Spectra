# Local charge metrics and coherent occupation proofs

The H6 charge metric now needs only **16 factors coupling sites at separation at most two**. Exact replay certifies the nearby localized Hamiltonian's ionic complement at **−6.264 Ha**, with an achieved rounded-metric floor of approximately **−6.263678275261903 Ha**. Both complete electronic energy proofs pass: the original interval has width **7.416980322924775e−12 Ha** with17 fresh response directions, and the localized hopping perturbation has width **7.227269910284279e−12 Ha** with18 fresh directions.

A new coherent occupation proof certifies the same gap with the exact unrounded product metric, without constructing a sector matrix or calling a determinant action. Its strongest tested tree has364 completed occupation leaves and certifies16 additional ionic configurations in shared branches. It nevertheless conditions transition bounds on all380 ionic configurations. This is a finite reduction in leaf checks and a different proof representation, not a resolution of enumeration or scalability.

## Restricting the metric to local pairs

Starting from unit factors and the actual380-by-380 complement matrix, a numerical scan limits the maximum spatial separation R in the one-site/two-site charge family:

| R | Factors | Achieved numerical row bound, Ha |
|---:|---:|---:|
| 1 | 12 | −6.296622773981219 |
| 2 | 16 | −6.263678258526106 |
| 3 | 20 | −6.255706713827292 |
| 4 | 22 | −6.248386099413061 |
| 5 | 24 | −6.235591867703217 |

Only the R=2 proposal receives a new exact export here. The other scan entries remain numerical proposals, not optimum or acceptance claims. The R=5 result reproduces the preceding full-pair proposal. The extended recipe retains the exact old dictionary when no distance is specified, and explicitly rejects invalid distances.

For fixed R this metric family has a number of parameters linear in site count. Construction and proof cost have not been shown to share that scaling. The Hamiltonian still contains long-range transitions, the valence reference grows, and response and rational-precision costs remain uncontrolled. Two-site charge functions can contain products of four occupation operators; this is not a closed two-particle-marginal parametrization.

## Transfer with a weaker excitation threshold

The R=2 gap shifts from−6.264 to−6.304 Ha under the independently verified hopping norm1/25. This is below the previous perturbed excitation proposal, so the new proof uses an offset of1/10000 Ha above its physical Rayleigh value instead of1/100. Fresh response discovery still requires18 directions. The smaller Temple denominator widens the certified interval to7.227269910284279e−12 Ha, while satisfying the requested target.

The original witness is inherited and rechecked; the unperturbed response is freshly discovered. The perturbed witness is refined again for24 steps on its actual Hamiltonian, and its response is freshly discovered. The tested perturbation is the same hopping between the first two localized orbitals as the previous experiment, not a new molecular system.

## Combine transition amplitudes before taking absolute values

For a normal-ordered CAR monomial, split creation and annihilation sets into changed orbitals C,A and common spectator orbitals S. The compiler checks the exact symbolic identity

\[
\text{monomial}=\pm T_{C,A}\prod_{i\in S}n_i.
\]

All terms with the same C,A share the same bare transition and fermionic phase. Their signed coefficients can therefore be combined into an amplitude polynomial. For the admitted one-/two-body Hamiltonian this polynomial is constant or affine in spectator occupations. H6's1818 Hamiltonian terms compile to1140 nonzero transition groups. Their absolute values reproduce the physical transition magnitudes, including cancellation between density-assisted and ordinary hopping.

`CoherentCharge` evaluates these amplitudes and exact positive metric ratios using occupation rules. It handles the complete valence projector implicitly: fixed spin populations are propagated, the number of allowed spin completions is computed combinatorially, and a site DP subtracts the singly occupied valence completions. Transitions whose entire source or target branch lies in the valence space contribute nothing to QHQ.

The branch proof applies the weighted row inequality

\[
H_{ss}-\sum_{t\in Q,t\ne s}|H_{ts}|\,u(t)/u(s)\ge\gamma.
\]

Every branch bound is conservative over its allowed configurations. Complete binary splits and exact counts establish coverage. No scalar gap or projected matrix is imported from the preceding proof.

## A checked charge/spin bound for diagonal energy

Let q_i=n_(i,alpha)+n_(i,beta)−1 and z_i=n_(i,alpha)−n_(i,beta). The verifier reconstructs the diagonal Hamiltonian exactly as

\[
V=c+\ell^Tq+\sum_i(U_i/2)q_i^2
+\sum_{i<j}C_{ij}q_iq_j+\sum_{i<j}S_{ij}z_iz_j.
\]

Since z_i²=1−q_i², each spin term obeys

\[
S_{ij}z_iz_j\ge-|S_{ij}|+
\tfrac12|S_{ij}|(q_i^2+q_j^2).
\]

This gives V≥c′+ell^Tq+q^TMq. Exact rational LDL verifies M−lambda I positive semidefinite at **lambda=44661/250000**. For any fixed rational anchor v,

\[
q^T(M-\lambda I)q\ge2v^T(M-\lambda I)q-v^T(M-\lambda I)v.
\]

The remaining separable site costs are minimized by a DP enforcing both spin populations and at least one doublon. Two anchors—zero and the charges fixed by the branch—give valid lower bounds whose maximum is used. At the root this certifies diagonal energy at least **−5.45637498498686 Ha**. This is a diagonal bound, not a spectral bound for the full Hamiltonian.

The charge bound alone does not prune any additional H6 configurations at the requested threshold. Its exact algebra and domain checks are retained, and the failed compression outcome is recorded.

## Preserve mutually exclusive transition conditions

Taking a separate worst-case penalty for every transition can charge one configuration for mutually incompatible source occupations. The next refinement multiplies each penalty by its source indicator before summing, forming a multilinear occupation polynomial on the free bits.

With k_alpha,k_beta occupied free modes, exactly

\[
\binom{k_\alpha}{a}\binom{k_\beta}{b}
\]

monomials of spin bidegree(a,b) evaluate to one. Summing that many smallest coefficients, including implicit zero coefficients, gives an exact lower bound without enumerating free assignments. This retains some cancellation between conditional penalties.

| Coverage method | Tree nodes | Completed occupation leaves | Q configurations certified in shared branches |
|---|---:|---:|---:|
| Coherent amplitudes and independent ratio bounds | 799 | 380 | 0 |
| Add the checked charge/spin diagonal bound | 799 | 380 | 0 |
| Also retain source-indicator polynomials | 783 | 364 | 16 |

The final proof still encounters all380 fully conditioned transition-source patterns. This counter prevents interpreting the absence of determinant-action calls as absence of configuration-dependent work. The tree also retains20 valence-only leaves. The accepted leaf polynomials have at most four nonzero coefficients, so the observed pruning is local and small.

## Complete energy integration and verification scope

The energy verifier accepts the coherent-tree proof both directly and in a transferred reference. It rebuilds the proof from the caller's actual Hamiltonian and threshold and requires P to contain the entire valence manifold. A different P, an overclaimed threshold, or a hidden Hamiltonian inside the proof recipe is rejected.

Replacing the matrix-based complement recipes in both R=2 energy certificates with the coherent tree preserves both interval widths and response dimensions. The response and witness are reused for this substitution and all final premises are rechecked. Full energy replay still processes400 configurations through its response and upper-witness operations.

Focused tests compare grouped rows against actual CAR actions, check partial-occupation inequalities, validate the charge/spin identity and tangent bounds, exercise exact source-cardinality algebra, and reject incomplete or incorrectly bound proofs. A64-mode number-operator fixture is covered by one symbolic branch with determinant actions forbidden; that trivial family does not establish interacting molecular scaling.

The full marginal regression suite passes **319 tests in444.067 seconds**. Nine new standard-library replays match every corresponding saved replay field: the local metric gap, its two matrix-based energy intervals and perturbed upper witness, three coherent gap-tree variants, and both energy intervals after substituting the coherent tree. The tree recipe contains no hidden Hamiltonian or explicit sector-state list; the energy verifier supplies the actual Hamiltonian and checks the complete valence reference.

The next unresolved step is to preserve the joint dependence of spectator amplitudes and metric ratios inside branch inequalities. Current independent maxima still force most branches to individual configurations. The saved `frontier_branch_obstruction.json` gives an explicit branch where the bound fails even though both exact child rows pass; it diagnoses the present relaxation and is not a no-go theorem for all occupation certificates.


Follow-up: [Joint occupation algebra](marginal_joint_occupation.md) repairs the saved two-mode branch and verifies a335-node simplex proof. The new proof still evaluates400 physical endpoints; global compression remains open.
