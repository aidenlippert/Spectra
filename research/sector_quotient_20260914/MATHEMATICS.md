# Restoring fixed-number cubic completion

The user's pasted request specifies this experiment. The attached derivation is
supporting mathematics; its conclusions were tested, not adopted as an authority.
All molecular statements in the report refer to the supplied frozen H8 Hamiltonian.

## Coefficient convention and exact reduction

For increasing k-tuples I,J, use
`E_IJ = a_i1† ... a_ik† a_jk ... a_j1`, without a factorial. A coefficient
matrix is indexed by k-orbital tuples, not N-particle determinants. The existing
CAR code stores annihilators in increasing order, so conversion multiplies
coefficients by `(-1)^(k(k-1)/2)`.

The insertion lift satisfies `hat(L_k V)=(Nhat-k)Vhat_k`. At k=2 this gives
`hat(L_2 V)=(N-2)Vhat_2` on the physical sector. If an aggregate positive
factorization expands as `S=S_low+hat(W)_3`, admissible exact completion is
`W=L_2 V`. The identity then reads

`H-b = S - (Nhat-N)Vhat_2 + lower-degree ideals`.

The unchanged verifier computes `R=H-b-S-(Nhat-N)X`, so **X4=-Vhat_2**.
Adding X4 while still requiring W=0 cannot help: L2 is injective for m>=5.

The signed contraction C is adjoint to L in the full coefficient Frobenius
inner product, including both off-diagonal entries. With K=C3(W), T=C2(K),
set `t=tr(T)/(3(m-2))`, `B=(T-t I)/(2(m-3))`, and
`V=(K-L1(B))/(m-4)`. Then `Z=W-L2(V)` has C3(Z)=0 and is orthogonal to
the lift image. The accepting equality is **Z=0**, not C3(W)=0.

The exact implementation checks the identities on all 325 pair-matrix units
for m=5,6 and sparse m=12,16 cases. An independently written occupation-bit
action checks a six-mode example on all 64 Fock inputs. Those 64 enumerated
inputs are an explicitly isolated algebra test; no molecular constructor or
checker enumerates N-electron determinants. The disjoint-triple counterexample
has zero contraction and squared norm 2 and remains nonzero. The molecular
integration test also confirms a disjoint-spatial triple constraint survives.

For signed factors `W=sum_t w_t v_t v_t^T`, remove each orbital from each triple
vector with its exterior sign to construct K on the pair coefficient space.
`||Z||_F^2=sum_tu w_t w_u (v_t^T v_u)^2-<V,K>` is checked exactly. This avoids
constructing W when signed factors are supplied. This campaign does **not**
claim a short signed W factorization for arbitrary molecular mixed cubic squares.

## Actual strong H8 diagnostic

The strong inherited proof is used only by `diagnostic.py`. No strong factors,
eigenvectors, or multipliers are read by the compact constructor. Its actual
expanded sextic W and quartic X4 satisfy `W+L2(X4)=-R6` exactly, including the
old proof's numerical-export residual. Projection also verifies
`V+X4=-projection(R6)` and `Z=-contraction_free(R6)` exactly. Consequently the
old proof uses a substantial lift component, but is not an exactly zero-residual
lift certificate. Its existing exact wedge witness pays the remainder.

That proof binds original H on the whole fixed-N sector and needs no separate
spin correction. The new compact construction below does use spin averaging,
so its distinct spin accounting must not be inherited implicitly from the old proof.

## Enlarged compact molecular cone

Keep the prior paired H8 dual48 maps V, selected through quartic MPS guiding
moments and prior compact coefficient duals. Each charge-compatible operator
dictionary receives an independent positive matrix Q. All entries and all cross
terms in those retained spans remain available. Partner matrices need not agree.

Use every degree-six coefficient equation and all symmetry-compatible quartic
number multipliers. This is the linear formulation
`S6+L2(X4)=0`; it is equivalent to collective Z=0 and corrected quartic matching.
No quadratic norm equality is inserted into an SDP. The model has 91,752 Gram
entries, versus 54,932 shared entries in the prior paired model. The separate
nonsinglet proof adds another 18,128 entries to either complete construction.

Equal partner Q matrices and X4=0 give the old paired construction as a feasible
special case. Odd adjoint cancellation transposes the partner Gram indices;
transposing leaves symmetric Gram evaluation unchanged. This is checked on every
actual molecular pair. The initial old numerical point matches the enlarged
equations to 4.29e-10 per coefficient. Exact exported factors still undergo full
residual accounting, including final rounding away from any nominal floating span.

The coefficient dual contains the additional degree-six variables. The guiding
MPS supplies no invented higher moments, and no MPS moment is imposed as the
unknown state's moment. The numerical probes used to check adjoints or condition
linear solves are random numerical vectors, not physical state moments.

### Number-dressed linear completion, after the unchanged-span test

The full H8 cubic dictionary had omitted linear words using the fixed-number
identities `Nhat a_i=(N-1)a_i` and `a_i† Nhat=N a_i†` on input N particles.
The selected compact cubic spans do not inherit this equivalence automatically.
An explicit probe found relative distances of 0.9483–0.9759 between the required
`Nhat a_i` directions and their projections into the retained spans.

The bounded follow-up appends exactly `Nhat a_i/4` for the 16 relevant modes,
and the corresponding adjoint columns, with quarter-integer coefficients. The
factor 1/4 is only a coordinate scale. Old columns remain unchanged; no rank
search, full-proof factor, or higher guiding moment chooses the new directions.
This raises the independent Gram-entry count from 91,752 to 96,904. The number
identities are checked as exact right-ideal identities at m=6,12,16.

This matters for interference with the other retained cubic operators too:
`C†(Nhat a_i-(N-1)a_i)=C†a_i(Nhat-N)`. The cross term is precisely a quartic
number multiplier, already restored in this campaign. Thus this is a completion
of the same fixed-number representation, not a separate physical model.

## Spin and complete lower-bound scope

Let T be the exact SU(2) average and Hs=T(H). Averaging preserves positivity,
particle number, and body degree, and commutes with L2. Thus the actual aggregate
operator for this singlet certificate is T(S); its sextic lift condition is
imposed collectively. A rational projector and a modular independence check
reduce 28,461 full coefficient rows to 8,533 independent invariant rows.
The projector has exact integer idempotence and rank equal to its exact trace.

The existing singlet checker allows quadratic spin multipliers and a Casimir
multiplier; this campaign does not silently extend those schemas. Its core H is
the exact auxiliary quartic polynomial prescribed by the sealed checker. The
quartic number multiplier is already admitted by the unchanged degree-three
CAR verifier.

The full lower endpoint is `min(L_singlet,L_MS1)-delta`, where the separate
MS1 certificate bounds every S>=1 multiplet of Hs. The input defect is
`delta=59/250000000000 Ha`, charged once. The upper endpoint is the same frozen
rational MPS upper used in the retained H8 comparison. An auxiliary singlet
number alone is never reported as a complete ground-energy interval.

## Numerical work is a proposal, not a proof

Matrix-free forward evaluation expands `V Q V^T` only inside the fixed cubic
operator dictionaries and applies sparse CAR maps. Its adjoint forms the exact
corresponding numerical congruence. It avoids storing the much larger dense
map from reduced Gram coordinates to all sextic coefficients. All 2,009,792
nonzeros of the sparse physical maps are counted as construction data.

Boundary-point optimization, numerical QR elimination, and fixed-bound
Douglas–Rachford feasibility are alternative numerical treatments of the same
enlarged cone. QR removes duplicate multiplier coordinates, not physical
equations; its floating defects remain subject to exact export/replay. Failed
convergence does not certify a cone obstruction.

The optional dense normal-matrix preconditioner is explicitly a large numerical
cost on coefficient equations. It is not a small proof, a determinant matrix,
or an accepting certificate. The report includes its allocation and construction.

`repair.py` can recover V from the actual rounded squares, replace X4 by its
exact projected value, and propose a residual wedge witness. It retains exactly
the remaining -Z in the averaged sextic residual. It never discards nonzero Z.
The sealed standard-library verifier independently re-expands the resulting
certificate and charges every residual and spin defect.

After the first measured 177.7-second repair, proposal-side spin projection was
accelerated using the fact that adjacent spin-index flips preserve canonical
order or create a repeated index and vanish. It agrees exactly with the sealed
projector on all 662 number-conserving six-mode matrix units and a sixteen-mode
sparse mixture. The accepting replay still uses the original sealed projector.

## Prior context, not a novelty claim

Fixed-number SOS and contraction-consistent decompositions are established
tools. This campaign tests a specific omission in Spectra's compact constructor.
Relevant primary sources checked for context:

- [Mazziotti, 2023: higher-order constraints and two-body unitary decomposition](https://arxiv.org/abs/2304.08570).
- [Rubin, Low and DePrince, 2026: weighted SOS and particle/spin constraints](https://arxiv.org/abs/2602.05069).
- [Lackner et al.: contraction-consistent two-particle reconstruction](https://arxiv.org/abs/1411.0495).

These sources do not validate this campaign's numerical results. The local
exact certificates and receipts provide that evidence. No general efficient
many-body solution or scaling theorem follows from this H8 experiment.
