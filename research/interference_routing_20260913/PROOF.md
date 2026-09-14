# Signed grounding and local fermionic routing

The attachment's routing implication is correct. The following derivations extend its sufficient test to diagonal corrections and make shared routing costs explicit. They use established congruence and signed-Laplacian algebra; no theorem of novel general molecular solvability is claimed. [Primary antecedents](SOURCES.md) include [Zelazo–Bürger](https://arxiv.org/abs/1408.2187) and [Chen et al.](https://eeqiu.people.ust.hk/wp-content/uploads/2021/09/Characterizing-the-Positive-Semidefiniteness-of-Signed-Laplacians-via-Effective-Resistances.pdf).

**Put the diagonal correction inside the routing problem.** Let H be a finite real symmetric matrix on the entire specified physical sector, let every phi_x be nonzero, and set D=diag(phi), ell_x=(H phi)_x/phi_x, w_xy=-H_xy phi_x phi_y. Then

\[
M_b=D(H-bI)D=L(w)+\operatorname{diag}(q),\qquad
q_x=\phi_x^2(\ell_x-b).
\]

Introduce a vertex g with an edge (x,g) of signed weight q_x for each x. Its augmented signed Laplacian is

\[
\widetilde L_b=
\begin{pmatrix}M_b&-q\\-q^T&\sum_xq_x\end{pmatrix}
=[I,-\mathbf1]^T M_b[I,-\mathbf1].
\]

The equality uses M_b 1=q. The rectangular map [I,-1] is surjective, and M_b is also a principal submatrix of the augmented matrix. Hence

\[
\widetilde L_b\succeq0\iff M_b\succeq0\iff H\succeq bI.
\]

Negative q_x are retained as negative edges. They are never declared harmless. Positive q_x can compensate for an indefinite L(w). This removes the separate requirements L(w)>=0 and ell_x>=b, while preserving the entire positivity obligation. The construction itself is an exact reformulation, not a complexity reduction.

**A tree turns globally interacting routes into a single matrix.** Choose an oriented spanning tree T of the augmented graph, or of an ungrounded connected component for the attachment's separated test. The implementation requires its edges to have positive weight and refuses if such a tree is unavailable. Let B_T be its incidence matrix. For each signed edge e, let f_e be its unique oriented tree path, so B_T f_e=b_e. Define

\[
C=\sum_e w_e f_e f_e^T.
\]

Then L=B_T C B_T^T. Since B_T has full column rank, L>=0 iff C>=0. Every cross term of every path remains in C. The verifier recomputes C exactly and performs rational LDL, rejecting negative pivots and nonzero coupling through a zero pivot. This basis change is equivalent to testing the full component's positivity; its dimension is the component size minus one. In particular, storing a short tree does not make a large component cheap to check.

The ungrounded molecular mode uses this test for L(w), sets its lower endpoint to min ell, and takes the minimum over **all** physical connected components. The grounded mode verifies an independently supplied rational b for each component. Both evaluate a separate nonzero integer trial vector's Rayleigh quotient exactly for the upper endpoint. Zeros are allowed in the upper vector but refused in D.

**The smallest simple-graph shared-route counterexample has four vertices.** Take positive edges 01,12,23 of weight 1 and negative edges 02,13 of weight -1/2. Each negative edge individually saturates its length-two path's resistance bound. Together their capacity is

\[
C=I-\tfrac12(1,1,0)^T(1,1,0)
     -\tfrac12(0,1,1)^T(0,1,1).
\]

For r=(1,2,1), Cr=-r/2. The node vector x=(2,1,-1,-2) has those tree differences and gives x^T Lx=-3, or Rayleigh quotient -3/10. Thus both isolated tests pass while the combined operator is indefinite. Three vertices cannot furnish this particular shortcut's counterexample: a connected positive simple subgraph already requires two of the three available edges, leaving at most one distinct negative edge. This counterexample does not refute the attachment's matrix inequality, which correctly includes joint effects.

The smallest mixed-slack example is H=[[1,-2],[-2,4]], phi=(1,1), b=0. Its local energies are (-1,2), yet H=(1,-2)^T(1,-2)>=0. Grounded routing accepts it exactly. A separate positive-slack example uses the attachment's H_t at t=3/5 plus I/5: its signed L is indefinite, but the combined H is PSD. Both are regression tests.

**A positive orbital graph cannot simply replace the configuration graph.** Consider four binary occupation modes and K_pq=I-SWAP_pq, with the parity-dressed CAR realization defined below. Put weight -1 on edges 01 and 23 and +1 on 02,03,12,13. The one-particle Hamiltonian is exactly uu^T for u=(1,1,-1,-1), hence PSD. In the two-particle sector ordered by occupation masks (3,5,9,6,10,12), however,

\[
v=(0,1,-1,-1,1,0),\qquad Hv=-2v.
\]

Equivalently, v is the product of antisymmetric pairs on 01 and 23. Writing H=sum_(p<q) K_pq-2(K01+K23), the total-spin-zero pair product has complete-graph energy 6 and subtracted energy 8. The direct integer matrix identity in [the receipt](/Users/aidenlippert/Documents/Spectra/results/interference_routing_20260913/sector_counterexample.json) proves the result without invoking spin representation theory.

This is the smallest number of modes for this failure within the occupation-exchange family: for at most three modes the only nontrivial sectors are one particle and its bit-complement image, which have identical exchange matrices. Four modes introduce an inequivalent two-particle sector. Thus exact one-particle routing is insufficient for a general many-particle exchange operator. This does not contradict the attachment, which explicitly starts in the complete physical sector, or the valid one-body implication A>=0 => dGamma(A)>=0. K contains density terms and parity-dressed hopping; it is not dGamma of the mode-graph Laplacian. The long edge's short parity program can contain degree-six CAR words, so this example is not asserted to be a two-body Coulomb molecule. Independent CAR expansion verifies every sector of the four-mode counterexample.

Congruence also requires the metric term. Already in one dimension, H=1, T=2, b=2 gives THT-bI=2>=0 while H-bI=-1. The correct expression T(H-bI)T=-4 rejects the false bound.

**A local CAR identity compresses whole families of routes.** For ordered spinless fermionic modes, define

\[
K_{pq}=n_p+n_q-2n_pn_q
-\left[a_p^\dagger\prod_{p<j<q}(1-2n_j)a_q+\mathrm{h.c.}\right].
\]

The parity product cancels the occupation-basis hopping sign. In that basis K_pq is I minus the operator that swaps the two occupation bits. It vanishes on 00 and 11 and is [[1,-1],[-1,1]] on 01,10. It is therefore PSD and annihilates the uniform amplitude in each fixed-N sector. These are mode-occupation pseudospins: no physical-site single-occupancy projection is assumed.

For three consecutive modes, restrict a K12+b K23-t K13 to the one-particle or two-particle sector. Each is a three-vertex path plus a negative chord. The empty and full sectors vanish. In either nontrivial sector the tree capacity is

\[
\begin{pmatrix}a-t&-t\\-t&b-t\end{pmatrix}.
\]

For a,b,t>=0 and a+b>0, positivity is equivalent to t<=ab/(a+b). For a=b=0 only t=0 passes. An exact boundary parameterization is

\[
a=t(1+\rho),\quad b=t(1+1/\rho),\quad \rho>0.
\]

Its capacity has nonnegative diagonal and zero determinant. Embed these even local CAR operators in an arbitrary longer chain. If a_i is charged to nearest edge i and b_i to nearest edge i+1, require

\[
b_{i-1}+a_i\le J_i
\]

on every edge, taking missing endpoint charges as zero. Summing the certified triangles and the unused positive K terms proves

\[
H_0=\sum_i J_iK_{i,i+1}-\sum_i t_iK_{i,i+2}\succeq0.
\]

The uniform fixed-N state is annihilated exactly. Overlapping templates cannot spend the same positive capacity twice. The next-neighbor positive off-diagonal entries form sign-frustrated triangles for nontrivial particle sectors, so a diagonal sign gauge cannot make all these couplings nonpositive.

For uniform J=1, t=1/4 and rho=1, this proves positivity with endpoint slack. Under the exact occupation-pseudospin identification, K_pq=2(1/4 I-S_p dot S_q). Thus this is the known frustrated ferromagnetic J1-J2 chain, not newly discovered solvable physics. [Richter et al.](https://arxiv.org/abs/0811.3549) discuss the corresponding one-quarter threshold. The implemented contribution is a rational routing compiler with explicit shared budgets and an inhomogeneous finite-grammar search.

**Discovery, accuracy, and bit costs remain separate.** The search picks rho_i from {1,1/2,2,1/4,4}. A dynamic program needs only the preceding ratio to check an edge budget. With q choices its arithmetic work is O(M q^2), storage O(M q) for backtracking. It is exhaustive only within that finite grammar; a refusal proves nothing about arbitrary positive real ratios.

Add the declared diagonal potential V=c+sum h_i n_i+sum v_i n_i n_(i+1). Then

\[
\min_{|s|=N} V(s)\le E_0(H_0+V)
\le \binom{M}{N}^{-1}\sum_{|s|=N}V(s).
\]

The upper bound is a rigorous uniform trial-state expectation because H_0 annihilates that state. The existing fixed-N chain DP computes both quantities exactly in O(MN) arithmetic work for this fixed range. It retains a particle counter and a constant-length occupation tail. The witness has O(M) rational entries; bounded-range products and sums have polynomial bit growth in M and input bit length. Timings do not bound arbitrary huge rational inputs; the recorded bit metric covers DP accumulators, not every temporary allocation.

Cheap evaluation does not guarantee a narrow interval. At half filling and V=delta sum n_i n_(i+1), delta>0, the minimum is zero and the uniform expectation is delta N(N-1)/M. Therefore this particular certificate has exact width delta(M/4-1/2). No choice of routing ratios changes this width. The remaining obstacle is the supplied upper/local-potential witness, even after sign frustration is handled.
