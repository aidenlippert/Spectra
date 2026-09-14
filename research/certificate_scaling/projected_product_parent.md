# Direct recognition of an interacting path parent

The implemented sufficient condition recognizes the original rational Hamiltonian
as a sum of positive local terms and supplies a nonzero common null state. Both
bounds are exact, with no occupation-sector enumeration. It is a restricted
parent family, not a claim that molecular Hamiltonians satisfy the condition.

After ordering the modes along a connected hopping path, define

\[
B_i=q_{i+1}a_i(1-n_{i+1})-q_i a_{i+1}(1-n_i),\qquad
H=cI+\sum_{i=0}^{M-2}g_i B_i^\dagger B_i,
\]

where all q_i and g_i are positive rationals. The exact CAR expansion is

\[
g_iB_i^\dagger B_i=\alpha_i n_i+\beta_i n_{i+1}
-(\alpha_i+\beta_i)n_i n_{i+1}
-t_i(a_i^\dagger a_{i+1}+a_{i+1}^\dagger a_i),
\]

with alpha_i=g_i q_(i+1)^2, beta_i=g_i q_i^2 and t_i=g_i q_i q_(i+1).
In the repository's ascending canonical annihilator convention, n_i n_j has
coefficient **minus one** on a_i† a_j† a_i a_j. The implementation reconstructs
actual CAR factor products; it does not trust a manually copied quartic sign.

## Recognition from H alone

Read onsite coefficients A_i and negative nearest-neighbor hoppings -t_i.
Set alpha_0=A_0 and recover deterministically

\[
\beta_i=t_i^2/\alpha_i,\qquad
\alpha_{i+1}=A_{i+1}-\beta_i.
\]

Require positive t_i and alpha_i, the final equation A_(M-1)=beta_(M-2),
and the observed density coefficient -(alpha_i+beta_i) on each edge.
Starting from q_0=1, recover q_(i+1)=q_i alpha_i/t_i and
 g_i=t_i/(q_i q_(i+1)). Finally expand the reconstructed sum g_i B_i†B_i
and compare it with the entire canonical original H. Any leftover term refuses.
No q, matching, factor, spectral solution, or upper witness is supplied.
An unknown chemical potential term mixed into every onsite coefficient is
outside this recognizer's condition; it is not silently removed.

For q=(1,2,1,3), g=(2,5,7), the exact data are
alpha=(8,5,63), beta=(2,20,7), t=(4,10,21), onsite=(8,7,83,7), and
physical density coefficients=(-10,-25,-70). The recurrence recovers q up to
its irrelevant common scale. A separate three-mode example q=(1,2,3),
g=(5,7) has onsite=(20,68,28), hoppings=(-10,-42), and densities=(-25,-91).

## Inferring path order and hopping signs

The public recognizer now infers the hopping graph from all nonzero one-body
offdiagonal terms. It requires exactly one connected path, traverses it from
the lowest-labeled endpoint, and fixes real hopping signs with a recursively
determined vertex gauge of plus or minus one. It applies that signed mode
permutation to EVERY original polynomial term, canonicalizes, and checks the
inverse transformation exactly. Then it runs the recurrence and full CAR
reconstruction above. No favorable orbital ordering or signs need be supplied.

A signed permutation of fermionic modes preserves CAR and the spectrum.
The null state is transported by its inverse, including the induced occupation
permutation signs. An independent four-mode N=2 test constructs those wedge
signs directly and checks H psi = c psi after a nontrivial path permutation and
sign gauge. Branches, cycles, disconnected graphs and non-Hermitian inputs
refuse. This permits relabeling a path; it does not turn a branching tree into one.

## Exact upper bound and non-Gaussian state

In every sector 0<=N<=M, the unnormalized state with occupation amplitudes
psi(S)=product(q_j for j in S), |S|=N, is nonzero and annihilated by every B_i.
If exactly one edge endpoint is occupied, the two contributions have equal
fermion signs: adjacent modes with the edge emptied have the same occupied
prefix. Their amplitudes cancel. Empty and doubly occupied edges contribute
zero because of the occupation projectors. Consequently E0=c in each sector.

The normalization is the elementary symmetric polynomial e_N(q_0²,...,q_(M-1)²),
computed with an exact descending dynamic program. At M=4,N=2, the amplitudes
violate the Slater two-form Pluecker relation: p01*p23-p02*p13+p03*p12 equals
the positive product of all four q values. Thus this fixed-N pure state is
non-Gaussian. This is an interacting control beyond hidden quadratic examples.

## Cost and executed checks

Graph recognition, signed relabeling and CAR reconstruction take O(T+M) arithmetic operations for T
input terms at fixed degree four. The optional explicit normalization takes
O(MN) rational operations and O(N) DP entries. There are M-1 factors, each with
four monomials. Arithmetic counts are not unit-cost bit complexity: input,
recovered parameter, and normalization bit lengths must also be counted.
The continued-fraction recurrence can be expressed using tridiagonal
continuants, so its rational bit growth is polynomial in M and input bit size;
products defining q and the elementary symmetric DP also have polynomial bit
bounds. This does not assert constant per-operation runtime.

`projected_product_parent.py` provides a python -S CLI. Six focused tests
cover independent occupation-bit action at N=2, rational parameters, every
allowed particle sector's recognition, unsupported terms, signs, boundaries,
invalid inputs, and the actual CLI. The Fock action is test-only. Production
recognition enumerates zero occupation states.

The executed nonuniform ladder M=4,8,16,32,64,128,256 is in
`results/certificate_scaling/projected_product_parent/ladder.json`. Every row
has exact interval width zero. At M=256,N=128, recognition plus normalization
took 0.0913 seconds in this run, with 255 factors, 1,020 factor monomials,
32,768 DP updates and a maximum recorded DP rational size of 508 bits.
These timings exclude fixture generation. All four frozen molecular fixtures
H4/H6/H8/H10 were tested and refused this condition.

The parent construction is related to the established stochastic-matrix-form /
Rokhsar-Kivelson approach; no novelty claim is made for that broad mechanism.
See [Castelnovo, Chamon, Mudry and Pujol (2005)](https://arxiv.org/abs/cond-mat/0502068).
The remaining research question is which broader chemistry-relevant structures
admit similarly direct discovery, with controlled error and total cost.

## Counterexample to an unchanged tree extension

For the four-mode star with edges (0,1),(0,2),(0,3), set q=g=1 and use
the same pair-projector factors on every edge. In the N=2 sector their sum
of squares has determinant 4 and characteristic polynomial
(lambda-2)^2 (lambda^2-4lambda+1)^2. Its minimum eigenvalue is 2-sqrt(3)>0:
there is no common null vector. A root-run independent occupation-bit test
checks determinant 4 and the exact matrix identity (H-2I)^3=3(H-2I).
This star has treewidth one. The counterexample refutes extending THIS
common-null-state construction based on treewidth alone; it is not a proof
that every certificate method fails on trees. Adjacent-mode fermion signs
are an essential assumption of the implemented path theorem.

The current inferred-path ladder is separately preserved under
`results/certificate_scaling/projected_product_parent/inferred_path/`.
At M=256 it takes 0.1457 seconds including graph recognition, signed
permutation and its exact round trip; the preceding 0.0913-second figure is
for the earlier natural-order core. Both have exact zero width. All frozen
molecular fixtures also refuse the expanded inferred-path condition.
