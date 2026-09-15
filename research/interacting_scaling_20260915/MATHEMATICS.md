# Interacting local and collective operator certificates

## Exact orbital changes and the earlier restricted family

Let R be a real rational spatial matrix with R^T R=I. Apply the same R to alpha
and beta orbitals. The CAR-preserving substitution is

\[
 a_{p\sigma}\longmapsto\sum_qR_{pq}a_{q\sigma},\qquad
 a^\dagger_{p\sigma}\longmapsto\sum_qR_{pq}a^\dagger_{q\sigma}.
\]

It is implemented by a number-preserving unitary on Fock space, without
constructing that space. It preserves N, N_alpha, all spin generators and
S^2. Normal-ordered degree and charge are preserved. Complete fixed-degree
operator spaces and complete degree-bounded ideal multiplier spaces are
therefore invariant. Dictionary changes act by invertible linear maps;
Gram matrices transform by congruence, preserving positivity.

The old molecular dictionary is not the unrestricted complete cubic space:
its pure triples are confined to two index groups, and its Gram matrices and
number multipliers use binary spatial parities. An orbital rotation crossing
those groups sends an allowed triple to a combination including forbidden
supports. A diagonal parity D becomes R^T D R, which is generally not diagonal.
Neither old index groups nor old bit masks may be inherited as local labels.
The exact coordinate tests exhibit this counterexample and check inverse
transport, particle number and S^2 preservation.

The small coordinate control uses complete cubic operators, no spatial
pruning, and complete two-body number multipliers on both sides. This makes
the positive/ideal cone invariant. Native coefficient-L1 residual penalties
are still basis dependent: the positive-cone equivalence alone does not prove
identical optima for the residual-penalized certificate family. An exactly
transported residual retains its original norm allowance by unitarity; using
the newly expanded coefficient norm is a different sufficient bound. This
distinction is explicit in any native-coordinate comparison.

For the rounded rotated Hamiltonian H_r, exact coefficient reconstruction gives
\(\|H_r-\Gamma(R)^\dagger H\Gamma(R)\|\le\eta_R\). Thus a completely
verified interval [L_r,U_r] implies [L_r-eta_R,U_r+eta_R] for the original
Hamiltonian. The original transform and both actual proof witnesses are checked;
a saved success flag is insufficient. Nuclear constants, if added, shift both
endpoints equally.

## Direct local and collective construction

Keep all global quadratic operator channels. Add cubic operator dictionaries
on overlapping spatial clusters, closed under both spin ladder actions. Also
keep a collective cubic dictionary containing all words supported on at most
two spatial orbitals. This is one coherent Gram family per spin channel:
products of words on distant pairs remain available.

With fixed local cluster width k and s spatial orbitals, a chain cover has
O(s k^6) local Gram entries. The collective pair-supported cubic list has
O(s^2) coordinates, hence O(s^4) Gram entries. Global quadratic channels also
have O(s^4) entries. These are representation counts, not runtime or accuracy
theorems. Increasing k or adding collective supports can erase that saving.

Each new cluster is added alongside existing blocks. Setting its Gram matrix
to zero recovers the previous family. Identical earlier ideal terms remain
available. The mathematical maximum certified lower bound cannot decrease
under these inclusions, although bounded numerical runs can return worse points.

The collective dictionary can also include every cubic word on a declared
three-, four- or five-orbital window. These words share a Gram block with the
earlier pair-supported words; this preserves their cross terms with other
windows and distant pairs. It is different from merely adding independent
window energies. For bounded width this adds O(s) words to the O(s^2) collective
list. It does not establish that bounded width suffices for a prescribed energy.

When a collective basis changes from V_old to V_new, the initializer solves
V_new A=V_old after expressing both in the same physical-word coordinates.
The equality is checked with integer arithmetic using the exact supported
rational denominators. It initializes Q_new=A Q_old A^T. The associated ideal
coefficients are matched by their actual polynomials. This preserves the old
operator identity mathematically. Floating-point transport is a proposal and
must pass the usual exact factor export and replay before any bound is accepted.

The number ideal uses global (N-N_target), multiplied by all global one-body
Hermitian operators, all long-range occupation products and local two-body
Hermitian operators on the declared clusters. It does not impose a fragment
charge. Alpha-number, scalar S^2 and spin-ladder freedoms are retained with
their original singlet-only interpretation. These multiplier restrictions are
part of the declared family and may be expanded; they are not equivalent to
the old complete global quartic multiplier space.

Generate coefficient rows directly from H, declared Gram products and ideal
products, then close under exact spin averaging and Hermitian conjugation.
No omitted Hamiltonian term is silently zeroed. The original checker forms
the complete CAR remainder and bounds it by its exact spin-averaged coefficient
norm. Sextic cancellation is enforced collectively across all blocks and ideals.

The MPS supplies numerical moments to the optimizer. It does not restrict the
unknown physical ground state. Accepted energy gain, construction cost and
complete replay determine which additions are retained.

## Coupled fragments

For a disjoint spatial partition, split the exact polynomial into terms whose
support lies in one fragment and terms touching several fragments:

\[
 H(\lambda)=H_{\rm within}+\lambda V_{\rm between},\quad 0\le\lambda\le1.
\]

The coefficient split and Hermiticity are checked exactly. Lambda=1 reconstructs
the original Hamiltonian, including hopping and charge-transfer terms. Only
global particle number is fixed. Intermediate lambdas are diagnostic models.
Any continuation state or retained proof is recorded and charged as continuation.

At lambda=0, a charge-aware local bound H_i>=ell_i+mu(N_i-n_i) adds to sum ell_i
on total N=sum n_i. A common supporting slope must cover every local charge
sector. At nonzero coupling, that argument alone gives no tight result; all
coupling must be represented or bounded. Exact cancellation of shared ideal
multipliers does not cancel the sum of local remainder allowances.

The separate disconnected control constructs each declared fragment Hamiltonian
on every local charge and spin-projection block. Exact local factor checks give
L_i(q)<=E_i(q)<=U_i(q). With all fragment charges permitted, min-plus convolution
over the total number gives

\[
 L_0=\min_{\sum q_i=N}\sum_i L_i(q_i),\qquad
 U_0=\min_{\sum q_i=N}\sum_i U_i(q_i).
\]

The lower follows from the block decomposition under commuting local number
operators. The upper is attained by a tensor product of the checked local
Rayleigh-quotient witnesses with the minimizing total charge. Dynamic programming
over the partial total charge avoids listing global determinant states or all
charge assignments. For a two-spatial-orbital fragment it enumerates 16 local
Fock labels, over all charges, rather than a full-system particle sector.

For any exact slope mu, set ell_i(mu)=min_q[L_i(q)-mu q]. Then
H_i>=ell_i(mu)+mu N_i, yielding the additional composable lower
sum_i ell_i(mu)+mu N. The checker evaluates candidate rational envelope
breakpoints and retains the best exact supporting value. The convolution bound
can be stronger when local charge energies are nonconvex. Both constructions
retain all local lower-bound losses. Neither applies when an interfragment term
is present; the control expressly refuses that case.

For the separate frozen-state response control, let chi_0 be the normalized
integer MPS constructed at zero coupling, and evaluate its **exact Rayleigh
quotient** e_ref on the fully coupled rotated Hamiltonian H_r. This is not merely
an upper estimate of that particular state's energy. If the original-model
ground interval is [L,U] and the verified rotation allowance is eta_R, then

\[
 \langle\Gamma(R)\chi_0,H\Gamma(R)\chi_0\rangle-E_0(H)
 \in[e_{\rm ref}-\eta_R-U,\ e_{\rm ref}+\eta_R-L].
\]

The state tensors, charges and normalization are held fixed; only the
Hamiltonian binding changes. This bounds the energy cost of retaining the old
state. It includes every type of variational readjustment and does not isolate
entanglement energy or prove that a particular chosen SOS block is necessary.

## Exact singlet trace and forced null spaces

For s spatial orbitals and N=2n, define

\[
 d_0={s\choose n}^2-{s\choose n+1}{s\choose n-1},\qquad
 \tau_0(A)=\frac{\operatorname{Tr}_{N,M=0}\mathcal T(A)
 -\operatorname{Tr}_{N,M=1}\mathcal T(A)}{d_0}.
\]

Each integer-spin S>=1 multiplet contributes once to each magnetic trace and
cancels. The S=0 component occurs only in the first. Haar averaging is required
for a general operator; commuting with S_z is not sufficient.

For a canonical word with identical creation/annihilation index sets of size k,
containing a alpha modes and k-a beta modes, its unnormalized magnetic trace is

\[
 (-1)^{k(k-1)/2}{s-a\choose N_\alpha-a}
 {s-k+a\choose N_\beta-k+a}.
\]

All other normal-ordered words have zero magnetic trace. Binomials outside their
range are zero. This evaluates bounded-degree traces without particle-sector
enumeration. Independent validation uses a 15-determinant spin projector only
as a separately charged test.

Tau_0 is normalized, positive, and vanishes on the singlet ideals. But its Gram
matrices can have forced null vectors. For example S_+ annihilates a singlet,
so tau_0(A^dagger S_+)=0 for every A. A mixture with tau_0 cannot correct a
negative value or a nonzero row in such a null direction at weight less than one.

A valid dual must satisfy normalization, Hermiticity, all declared singlet
ideal equations, positivity on every declared operator Gram, and the dual
coefficient-box conditions for the actual residual rule. Its objective is an
upper bound on the lower certificates attainable by that family. A physical
trace is a repair ingredient, not an accepted sharp family obstruction.

### Implemented affine repair

The repair exports the actual rational operator coordinates, ideal polynomials
and twirled coefficient residual rule of one declared family. Write its Gram
functional as G_C(y), with y(1)=1, y(I_j)=0 and |y_w|<=1. The last box is the
dual of the coefficient-L1 remainder penalty, including Hermitian pairing.

For each exact null vector v of G_C(tau_0), impose G_C(y)v=0 before mixing.
The trace is faithful on the singlet subspace, so such a polynomial combination
annihilates every singlet. These equations select a physically compatible face;
they need not be consequences of a previously restricted collection of ideals.
Adding them to the repair search cannot invalidate a resulting feasible dual.
Exact sparse affine elimination enforces those equations, the ideals,
normalization and any trace-saturated coefficient coordinates. Then search
y_t=(1-t)y_repaired+t tau_0 over a declared rational grid. Every accepted
candidate undergoes fresh exact Gram PSD, coefficient-box and ideal checks.

If D=y_t(T(H)), weak duality gives b-||R||_1 <= D for every certificate in
this exported singlet family. The physical-H lower endpoint from that family
is therefore at most D-delta_spin, before any orbital allowance. This bounds
the attainable lower certificate; it is not a lower bound on the physical
ground energy. A family limitation relative to a frozen upper U would require
U-(D-delta_spin)>0.0016 Ha, with all other conventions matched. Merely finding
an exact dual with a large objective, or failing to repair one, proves no such
limitation. Full fixed-N traces cannot replace tau_0 in these singlet equations.

## Exact spin-pattern reuse and nonsinglet completion

For a normal-ordered word w, relabel its distinct spatial indices in increasing
order to 0,...,k-1, preserving each alpha/beta label. A spatial relabeling P
commutes with the global spin representation. Therefore T(Pw)=P T(w).
Caching the exact rational twirl of this bounded-degree pattern reuses algebra
without dropping coefficients. Order-preserving relabeling avoids an additional
fermionic permutation sign. The H10 comparison rebuilds all coefficient maps
and finds bit-identical outputs before claiming a preparation speed improvement.
The independent accepting checker continues to reconstruct the original words.

For spin-invariant even-N H, every S>=1 multiplet has an M_S=1 member. A lower
bound on the entire M_S=1 sector therefore bounds all nonsinglets. Construct
that screen directly from quadratic operators first, with exact occupation-
trace initialization. If insufficient, add local cubic windows and then
coherent pair-supported cubics. Use N_alpha-N/2-1 as the magnetic ideal, and
do not use the singlet S^2 or ladder ideals. A magnetic-screen lower L_1 and
singlet lower L_0 yield min(L_0,L_1)-delta_spin for the original spin-defective
operator. This screen is part of discovery and final exact replay, not a free
external assumption.

See [primary sources](SOURCES.md) for established embedding and coarse-graining
precedents. None of these identities proves favorable scaling for molecular
Hamiltonians.

## Complete accepted interval

Write H_s=T(H_r) and let delta_spin bound ||H_r-H_s||. Suppose the exact
singlet and magnetic replays give L_0=b_0-epsilon_0 and L_1=b_1-epsilon_1
for H_s, where each epsilon includes the actual reconstructed remainder and
any charged factor removal or rounding. Let q be the exactly evaluated
Rayleigh quotient of the retained, normalized fixed-N MPS on H_r. Then

\[
 E_0(H)\in
 \left[
 \min(L_0,L_1)-\delta_{\rm spin}-\eta_R,
 \ q+\eta_R
 \right].
\]

The spin defect is charged once on the lower side. It is unnecessary on the
upper side because q is evaluated on H_r itself, rather than H_s. The orbital
allowance is charged at both endpoints. All comparisons use the resulting
exact interval width, not a numerical primal/dual residual or a solver stopping
flag. Intermediate coupling models use their own exact H(lambda); transfer to
the original molecular input is applied at the fully coupled endpoint.
