# Fractional/overlapping anticommuting covers: primary-source notes

## Question and proposed certificate

For a real Hermitian Pauli expansion H = sum_i c_i P_i, choose pairwise
anticommuting groups G and rational coefficients x_(G,i) satisfying exact
reconstruction c_i = sum_(G containing i) x_(G,i). Then

\[
\|H\|_\infty\leq\sum_G\sqrt{\sum_{i\in G}x_{G,i}^2}.
\]

Proof: each group sum squares to `(sum_i x_(G,i)^2) I`, and the triangle
inequality combines groups. Minimizing this bound subject to reconstruction
is a second-order cone program. Keeping the original partition feasible while
adding overlapping decompositions cannot worsen the mathematical optimum.
A finite optimizer may still cost more or return a worse candidate; the checker
must validate the actual submitted split. This is an elementary derivation,
with no claim of priority or demonstrated computational headroom.

## Closest primary literature

1. **Wu, Sun, Huang, Yuan, “Overlapped grouping measurement: A unified framework for measuring quantum states,” 2021.** [arXiv:2105.13091](https://arxiv.org/abs/2105.13091), with published version linked from the arXiv record.  The paper explicitly allows *overlapped groups* of compatible Pauli measurements and introduces coefficient splitting/optimization to reduce estimation cost.  This is the closest precedent for the reconstruction constraint and overlapping-group optimization.  Its objective is measurement variance/shots (state-dependent covariance and estimator design), so it does **not** by itself prove an operator-norm upper bound of \(\sum_G\|x_G\|_2\).

2. **Huang, Kueng, Preskill, “Predicting many properties of a quantum system from very few measurements,” 2020.** [Nature Physics](https://www.nature.com/articles/s41567-020-0932-7), [arXiv:2002.08953](https://arxiv.org/abs/2002.08953).  Classical-shadow analysis gives a distinct state-estimation/sample-complexity framework for many observables.  It is useful context for why measurement grouping and operator-norm certification should be kept separate; it does not supply the proposed decomposition certificate.

3. **Sarkar and van den Berg, “On sets of commuting and anticommuting Paulis,” 2019.** [arXiv:1909.08123](https://arxiv.org/abs/1909.08123).  This is a primary structural study of maximal commuting and pairwise-anticommuting Pauli sets, including efficient generation of maximal anticommuting sets.  It supports the combinatorial side of selecting candidate groups, but does not optimize split coefficients or state the SOCP certificate.

4. **de Gois, Hansenne and Gühne, “Uncertainty relations from graph theory,” Phys. Rev. A 107, 062211 (2023).** [primary paper](https://arxiv.org/abs/2207.02197).  For dichotomic observables, the paper relates sums of squared expectations/variances to graph parameters, including the Lovász number of the anticommutativity graph.  This is a stronger global graph/SDP viewpoint than a clique cover in some instances, but it bounds expectation/uncertainty quantities and is not automatically an upper bound on \(\|\sum_i c_iP_i\|_\infty\) for arbitrary coefficients.

5. **Lovász, “On the Shannon capacity of a graph,” IEEE TIT 25 (1979).** [IEEE DOI](https://doi.org/10.1109/TIT.1979.1055985).  This is the primary source for the \(\vartheta\) function.  Its SDP and graph sandwich theorem motivate Lovász-theta relaxations on (anti)commutation graphs, but applying \(\vartheta\) to a Pauli operator-norm bound requires an additional operator-specific argument; the graph parameter alone is not the proposed weighted norm decomposition.

6. **McNulty, “A Graph-Theoretic Approach to Quantum Measurement Incompatibility,” 2025.** [arXiv:2511.15954](https://arxiv.org/abs/2511.15954).  This recent primary preprint develops anti-commutativity graphs and bounds incompatibility robustness using Lovász number, clique number, and fractional chromatic number.  It is relevant evidence that fractional graph covers are established in a neighboring measurement-incompatibility objective.  It does not establish the coefficient-reconstruction SOCP or identify its value with a Hamiltonian operator norm.

## Assessment

The *ingredients* are established: exact coefficient splitting over overlapping measurement groups (Wu et al.), anticommuting Pauli clique structure (Sarkar--van den Berg), and graph/fractional-cover/\(\vartheta\) bounds for uncertainty or incompatibility (Hansenne et al.; McNulty).  I did not find a primary paper that names or proves the exact weighted operator-norm certificate
\(\inf_{Ax=c}\sum_G\|x_G\|_2\) for overlapping pairwise-anticommuting Pauli groups.  The elementary proof above is valid independently of a priority claim. The limited source search does not establish novelty, and measurement-shot grouping alone does not prove an operator-norm advantage.

The baseline is strong in the sense needed for safe use: every feasible split gives a rigorous upper bound, and the SOCP can be solved without constructing the full (R^2) matrix.  Claims of *tightness*, superiority to the partition bound on a target instance, or comparison with a Lovász-theta/SDP bound still require explicit numerical or analytical validation.  Candidate group enumeration can itself be expensive; using a selected family of cliques keeps the certificate scalable but makes the result family-dependent.

## Scope/assumptions to carry into implementation

- (P_i) are Hermitian Pauli strings and (c_i\in\mathbb R); identity terms are allowed but cannot belong to a nontrivial anticommuting group.
- Each (G) is pairwise anticommuting.  Commuting measurement compatibility is a different relation from anticommutation and must not be substituted silently.
- Reconstruction is exact coefficient-wise.  Approximate or probabilistic coefficient allocation would require a residual term.
- The certificate is an operator-norm bound obtained by triangle inequality; measurement-shot or variance improvements do not imply operator-norm improvements.
