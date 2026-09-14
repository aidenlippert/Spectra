# What is exact, and what remains a numerical hypothesis

The new accepting schema is `mechanism_spin_subspace_v1`. Its expansion is
the preceding completed-spin checker with a new schema name and a local
baseline checker whose explicit cap is sixteen spin orbitals. Hamiltonian
and tail bindings, rational factors, number-ideal validation, exact adjoints,
and the degree-four residual gate are retained. Frozen historical modules
were preserved. The enlarged cap is tested against a sixteen-mode diagonal
Hamiltonian with a known ground energy, plus malformed inputs and the old
H6 certificate.

For every accepted case, exact expansion gives

    Hret = b I + S + (N-n)D + R,

where S is a sum of Hermitian squares, including both B^dagger B and
B B^dagger built from the same exact B. On the fixed-N sector the ideal
term vanishes. The existing coefficient norm argument gives

    E0(Hret) >= b - ||R||_coefficient,1.

The separately checked collective tail gives Hret+ell I <= H, so
L=b-||R||_1+ell is a molecular lower. Exact rational trial-state evaluation
supplies U. No eigensolver status, approximate rank or proposed rule is an
acceptance condition for L <= E0(H) <= U.

## A genuine identity inside the candidate rule

For a spatial pattern L_k, let Q_k be its spin-summed density and let
T_ki=(L_k a)_i. The CAR imply

    [Q_k,a_i] = -T_ki,
    [T_ki,Q_k] = (L_k^2 a)_i.

Thus, for Hret=h_tilde+(1/2)sum_k lambda_k Q_k^2,

    [Hret,a_i] = [h_tilde,a_i]
                  -sum_k lambda_k Q_k (L_k a)_i
                  -(1/2)sum_k lambda_k (L_k^2 a)_i.

The alpha=1, spin-summed portion of the rule is the cubic contraction in
this identity. `identity.py` checks every annihilation mode of the training
fixture by exact polynomial subtraction. The alpha=1/2 version and the
separate spin components are candidate extensions, not identities equating
them with the commutator.

The selected directions are normalized and rounded at 10^10. The literal
identity above refers to the exact unrounded alpha=1 operators; it is not
silently asserted for rounded vectors. Each rounded vector remains an
explicit exact combination of allowed generators and receives the ordinary
residual check when used in an energy certificate.

## Structure measurements are not energy guarantees

If C contains exact generator coefficients converted to floating point and
F,D are the accepted root and direction factors, the physical coefficient
matrix of proof operators is W=C D^T F^T. With M the uniform fixed-N trace
matrix, Z=M^(1/2)W measures their aggregate squared operator size. Singular
values of Z and projections onto a proposed template subspace describe this
particular proof's size distribution. They do not bound the energy loss
after deleting directions. Small omitted operators may be essential to
coefficient cancellations. Likewise, an approximately shared spin-flip Gram
does not establish an exact coefficient identity between two fitted blocks.

The recipe uses two scalar exponents selected solely from those training
diagnostics. Its coefficients are recomputed from each case's own spatial
patterns. Numerical duplicate removal is logged, and no claim of complete
family coverage follows from it. Held-out solves use full coupling or only
separate nonnegative contributions within that fixed recipe; they receive
the same exact verifier.

## Interpreting the remaining H6 interval

Let P <= E0(H) <= U be the independently audited physical interval and
C be the accepted full-family dual ceiling. The exact relations are

    E0(H)-L in [P-L,U-L],
    E0(H)-L_family >= P-C,
    U-L_family >= U-C.

Here P-C=1.458470 mHa and U-C=1.513465 mHa. These do not rule out the
1.6 mHa target. C-L=1.901956 mHa is an unresolved gap, not a measured amount
that better selection can necessarily recover; C need not be the tightest
family ceiling.

The non-overlapping identity in `interval_diagnosis.json` is

    U-L = (U-E0(H)) + (E0(H)-E0(Hret)-ell) + (E0(Hret)-b) + ||R||_1.

The first two terms are bounded respectively by 0.054994 and 0.108331 mHa;
the exact residual is 0.009140 mHa. Hence E0(Hret)-b is at least 3.242955
mHa for this certificate. The scalar b itself is not asserted to be a valid
lower before its residual is paid. The separate 0.006569 mHa loss relative
to the reported floating objective is a diagnostic, not an additional
penalty to add again.
