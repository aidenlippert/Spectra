# Three-patch projector correction

Suppose exact local checks establish, on three four-site patches,

\[
 A_i\succeq g_i(I-P_i)-\epsilon_i I,
 \qquad P_i=|u_i\rangle\langle u_i|,
 \quad g_i\ge0.
\]

The projectors are embedded into the full ladder by tensoring with identity on
the sites outside patch (i). Then

\[
 \sum_i A_i\succeq \sum_i(g_i-\epsilon_i)I-\sum_i g_iP_i.
\]

This is a global lower bound without constructing the eight-site Hilbert
space. Let (c_{ij}) be a certified upper bound for the norm of the overlap
operator between the ranges of the embedded projectors. For rank-one patches,
the nonzero spectrum of (sum_i g_iP_i) is bounded by the (3\times3)
comparison matrix with diagonal (g_i) and off-diagonal

\[
 \sqrt{g_i g_j}\,c_{ij}.
\]

Thus

\[
 \sum_i A_i\succeq
 \left(\sum_i(g_i-\epsilon_i)-\lambda_{\max}(G)\right)I.
\]

For disjoint outer patches, (c_{02}=0). For adjacent patches sharing one
rung, (c_{12}) is the operator norm of the contraction of the two local
bipartite vectors. It can be certified from the (16\times16) contraction
matrix (K_{ac}=\sum_b\overline{u_{ab}}v_{bc}), with outward rational norm
rounding. This is the only overlap calculation; no global determinant list is
needed.

To avoid unverified square roots, use weighted Young inequalities. For each
edge choose a positive rational (s_{ij}) and replace

\[
2c_{ij}\sqrt{g_i g_j}
\le c_{ij}\left(s_{ij}g_i+s_{ij}^{-1}g_j\right).
\]

The resulting rational row-sum/Gershgorin bound is a certified upper bound on
\(\lambda_{\max}(G)\). Optimizing the (s_{ij}) is numerical proposal only;
the checker verifies the final rational inequalities. The exact local
36-dimensional PSD checks establish the patch inequalities, while this small
comparison bound composes them globally.

This is an established projection-overlap lemma applied constructively, not a
new many-body breakthrough. It becomes useful only if the correction remains
positive after the ε terms and if the patch overlap graph has bounded degree.
The adversarial failure is a chain of nearly parallel local projectors: then
\(\lambda_{\max}(G)\) grows with the number of patches even when every local
gap (g_i) is large. Local gaps alone therefore do not prove a thermodynamic
global gap.
