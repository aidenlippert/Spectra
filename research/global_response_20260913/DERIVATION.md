# Global response and terminal transport

All statements concern the exact rational input Hamiltonian in its specified
fixed-particle sector. They do not bound molecular modeling error.

Write `M = H − bI`, with orthogonal retained and eliminated projectors `P,Q`.
In these subspaces, `M = [[A,B],[B*,D]]`. Exact orbital and local CAR
certificates establish `δI ≤ D ≤ MI` and `||B|| ≤ g`.

For H6 the simultaneous eliminated space is the union of double occupation
of the last two spatial orbitals. Its orthogonal decomposition is
`Q = Q1 + R2`, where `R2 = (I−Q1)Q2`; its complement has 532 labels
combinatorially. No list of those labels is needed to certify the response.
The other geometries and H8 use the certified single-orbital fallback.

If the two diagonal eliminated blocks have upper bounds `M1,M2`, and the
internal off-diagonal coupling has squared norm at most `ν`, then

\[
M_Q\le\frac{M_1+M_2+\sqrt{(M_1-M_2)^2+4\nu}}2.
\]

Dropping this internal coupling would be unsound. For the coupling from the
final retained space, the two output ranges are orthogonal, so

\[
\|[B_1;B_2]\|^2
=\|B_1^*B_1+B_2^*B_2\|
\le g_1^2+g_2^2.
\]

Both input bounds are restricted to the same final `P`. The second is an
upper bound on `R2 H P` obtained by restriction of `Q2 H (I−Q2)`.
All scalar square roots and endpoints are rounded outward rationally.

**Response and complete lower condition.** Set `c=(M+δ)/2`, `a=(M−δ)/2`,
`z0=c/a`, and use the polynomial

\[
r_k(d)=\frac{T_k((c-d)/a)}{T_k(z_0)},\qquad
p_k(d)=\frac{1-r_k(d)}{d}.
\]

The numerator vanishes at zero, so `p_k` is a polynomial of degree `k−1`.
With `X=p_k(D)B*` and the tall embedding `V=[I;−X]`,

\[
K=V^*MV=A-BX-X^*B^*+X^*DX
  =A-BF_k(D)B^*,\quad
F_k(d)=\frac{1-r_k(d)^2}{d}.
\]

The exact Schur complement is

\[
S=A-BD^{-1}B^*=K-R^*D^{-1}R,\quad R=B^*-DX.
\]

Since `||R|| ≤ g/|T_k(z0)|`,

\[
S\succeq K-\eta I,\qquad
\eta=\frac{g^2}{\delta T_k(z_0)^2}.
\]

Therefore `K−ηI ≥ 0` together with `D>0` proves `H≥bI`. A certified response
alone does not establish that terminal inequality.

**What the inherited terminal proof establishes.** If a separately checked
global SOS certificate proves `H≥L` with `γ=L−b>0`, then

\[
K=V^*MV\succeq\gamma V^*V\succeq\gamma I.
\]

Consequently `K−ηI ≥ (γ−η)I`. This is a direct positive-congruence argument;
no determinant matrix or inverse is needed. The inherited certificate is
essential to this proof. Transporting it is not independent discovery of a
small terminal factorization, and cannot be credited as new bound improvement.

For the old nested construction, `L1=K1−η1I`. The second `K2` is formed from
`L1`, so the actual terminal is `K2−η2I`, and its transported margin is
`γ−η1−η2`. Subtracting `η1` again from `K2` would double-count it.

The terminal-aware variant selects `η_budget=min(γ/2,0.001 Ha)` and takes the
first Chebyshev order meeting that budget. H6 needs order 87, giving 175
base-H actions. This rule uses the inherited terminal margin. The separately
frozen transfer rule continues to use the uniform `0.000001 Ha` budget.

**Congruence metrics.** For an invertible, nonunitary square `T`, a negative
error measured in the transformed identity metric obeys

\[
T^*(H-bI)T\succeq-\varepsilon I
\Longrightarrow H\succeq bI-\varepsilon T^{-*}T^{-1}
\succeq(b-\varepsilon\|T^{-1}\|^2)I.
\]

The relative condition `T*(H−bI)T ≥ −εT*T` instead gives `H≥(b−ε)I`.
The tests include an exact diagonal counterexample to omitting the inverse
norm. The positive terminal-transport argument above uses `V*V≥I` and does
not make that omission.

**Equivalent factoring and execution.** The scalar identity

\[
T_k(x)^2=\tfrac12(T_{2k}(x)+1)
\quad\Longrightarrow\quad
F_k(d)=\frac{T_{2k}(z_0)}{T_{2k}(z_0)+1}\,p_{2k}(d)
\]

is an exact operator identity because all powers within it are of the same
`D`. The surrounding noncommuting projectors and Hamiltonian factors keep
their order. It gives a single recurrence for `F_k`, but does not reduce
the physical nested count `(2k1+1)(2k2+1)=12,667`. The simultaneous response
instead has `2k+1=237` actions at the uniform budget. Actual per-RHS counters
verify these counts, without treating a block of RHS vectors as one action.

**Collective coupling rule.** Given nonnegative operator-block norm bounds
`n_ij` and positive row/column weights `p_i,q_j`, the weighted Schur test is

\[
\|B\|^2\le
\left(\max_i\sum_j n_{ij}q_j/p_i\right)
\left(\max_j\sum_i n_{ij}p_i/q_j\right).
\]

This follows by weighted Cauchy–Schwarz and summing the block-vector norms.
The bounded search checks at most `3^7=2187` powers-of-two weight choices.
It improves the old internal coupling majorant, but is not a terminal
positivity proof. Reusing that map for a different exceptional complement
would be invalid; that earlier attempt was rejected.

**Compact upper controls.** A Slater determinant with occupied columns `C`
has norm `det(C*C)` and density projector `R=C(C*C)^{-1}C*`. Disjoint rational
orbital rotations make the Gram matrix diagonal. Wick contractions evaluate
every one- and two-body Hamiltonian term without a many-body vector.
For CH2 the common-orbital, nested alpha/beta occupations certify pure spin;
spin projection alone would not suffice.

The response lift uses `χ=QH|HF>` with `Q|HF>=0` and

\[
U(\alpha)=\frac{e_{HF}-2\alpha\|\chi\|^2+
\alpha^2\langle\chi,H\chi\rangle}
{1+\alpha^2\|\chi\|^2}.
\]

The rational grid includes zero. Exact replay checks its moments and positive
norm. This implementation materializes sparse excited determinants during
the contractions; their count remains in the ledger. It never constructs
the complete fixed-N basis, but must not be described as using no many-body
labels. Independent small explicit CAR calculations test both upper formulas.

**Relation to existing work.** Elimination and spectral estimates are
established Feshbach–Schur machinery; the new measured result here is the
particular jointly eliminated response and its execution cost.
[Dusson, Sigal and Stamm](https://arxiv.org/abs/2105.02058).
Noncommutative rational realizations offer a relevant representation
precedent, without implying cheap inversion or molecular positivity.
[Volčič](https://arxiv.org/abs/1505.07472),
[Porat and Vinnikov](https://arxiv.org/abs/1905.11304).
Renormalization-based lower bounds are a relevant comparator for compressed
consistency constraints; this campaign did not run that external method,
so it makes no competitive runtime claim against it.
[Kull, Schuch, Dive and Navascués](https://arxiv.org/abs/2212.03014).
