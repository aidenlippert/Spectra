# Resource-bounded discoverability of physical law

This benchmark studies mechanism discovery under an explicitly bounded,
operational sensor. For three qubits, the candidate dictionary is every
nonidentity Pauli string of weight at most two: (3\binom31+9\binom32=36)
terms. A probe is a Pauli observable (Q) and Pauli-biased preparation

\[
\rho_S=(I+S)/2^3.
\]

For candidate Hamiltonian term (P), the ideal short-time response is the
coefficient of (S) in (i[P,Q]). The learner receives the resulting feature
matrix and measured responses, never the hidden support labels. The simulator
adds uniformly bounded deterministic noise, so the result is an ideal calibrated
linear-response sensor rather than full state or process tomography. Identity
terms are excluded because a scalar Hamiltonian offset is unobservable in this
sensor.

The experiment has three sequential stages. Stage 1 learns a public base
mechanism. Stage 2 introduces a previously absent noncommuting term and tests
warm-start OMP, scratch OMP, retrieval-only refitting of the old support, frozen
stage-1 predictions, and full-dictionary least squares. Stage 3 adds a second absent interaction. Every stage has a separate changed-coefficient context; held-out measurement observables have a different Pauli weight, so the evaluation probes are disjoint from acquisition probes. All methods receive the same sensor records where
possible; feature computations, records, and fitting work are counted. A
weight-three out-of-dictionary term supplies the misspecification case: a valid
system must report residual failure or abstain rather than force a nearest
in-dictionary explanation.

“Discovery” is support recovery plus predictive transfer under new probes,
coefficients, and combinations. Retrieval may match warm-start OMP and that is
a legitimate result; no compounding capability is inferred merely from cache
reuse. A small fixed Pauli vocabulary also cannot establish open-ended law
discovery. Any claim beyond this closed scope requires an explicit vocabulary
extension and an identifiability argument.

For a finite family (M_1,\ldots,M_K), let (P_i^e) denote the classical
outcome distribution of an allowed experiment (e). If one experiment has a
known pairwise separation

\[
\min_{i\ne j}\operatorname{TV}(P_i^e,P_j^e)\geq\Delta,
\]

then pairwise tests can identify the model with logarithmic dependence on
(K/\delta) and quadratic dependence on (1/\Delta); for bounded event
statistics, a Hoeffding-plus-union-bound sufficient scale is
(n=O(\Delta^{-2}\log(K/\delta))). This is an upper bound under the stated
experiment and does not imply that TV alone gives a universal lower bound with
the same constants: KL/Hellinger behavior and the outcome model matter. If all
allowed experiments have zero separation for a pair, no finite data can
distinguish that pair; Le Cam’s two-point method supplies the corresponding
minimax obstruction. Adaptive policies must count experiment-selection cost and
use the best budgeted separation actually achieved.

Primary sources: Evans, Harper & Flammia, [Scalable Bayesian Hamiltonian
Learning](https://arxiv.org/abs/1912.07636), for adaptive Hamiltonian learning
and conditioning/sample bounds; Zhang et al., [Identifiability Guarantees for
Causal Disentanglement from Soft Interventions](https://arxiv.org/abs/2307.06250),
for intervention-dependent identifiability and equivalence classes; and
Polyanskiy & Wu, [Dualizing Le Cam’s
Method](https://arxiv.org/abs/1902.05616), for two-point minimax lower-bound
logic.
