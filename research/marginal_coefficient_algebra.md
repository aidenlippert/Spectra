# Coefficient-space SOS proposer: audit and fair comparison

The coefficient-space formulation should use a fixed canonical CAR word basis.
For each dictionary word `w_q`, compute its normal-ordered coefficient vector
`v_q` (body ranks 0–3). A PSD Gram matrix `Q` contributes the coefficient
vector of `sum_{q,r} Q_qr w_q†w_r`; a number-sector gauge contributes the
linear image of `(Nhat-N)X`, with `X` restricted to number-conserving words of
degree at most four. The equality map is therefore a sparse rational linear
map from `(b,Q,X)` to the target Hamiltonian coefficients.

The SDP proposal is:

`maximize b subject to Q >= 0 and A(b,Q,X)=coeff(H)`.

Keep `X` Hermitian in the real coefficient basis. Never impose PSD on the
surviving two-body coefficient matrix: positivity belongs to `Q`; the target
operator is obtained only after CAR expansion. The number gauge is algebraic
bookkeeping on the fixed sector, not a positive term.

Comparison with the earlier 20-dimensional fixed-sector solver must separate
three effects. The sector solver permits arbitrary finite-sector factor rows
and can therefore exploit many-body identities hidden from a degree-bounded
coefficient dictionary. The coefficient proposer is a restricted, portable
certificate family and may produce a weaker lower bound. Conversely, its
certificate size scales with canonical body-word counts (roughly
`binom(M,k)^2`) rather than `binom(M,N)^2`, and its verification is independent
of sector diagonalization.

The map has substantial gauge redundancy: CAR identities, Hermitian conjugacy,
particle-number contractions, and null dictionary combinations create dependent
rows and nonunique `(Q,X)`. Build a sparse exact row basis before numerical SDP,
or use QR/SVD only to identify a candidate rank and then recheck every proposed
identity with exact rational CAR arithmetic. Near-zero eigenvalues of `Q` must
be treated as uncertified until rationalized and PSD-checked; rounded `X` must
be included in the residual norm.

The strongest fair experiment fixes the same word dictionary, denominator,
solver tolerance, and Hamiltonian coefficient convention across baseline DQG,
T1/T2 cubic blocks, and the coefficient-space lift. Report proposal objective,
exact replay residual, and certified lower bound separately. A failed exact
replay is rejection, never a widened error bar. No conclusion about generic
representability follows from a successful small-M run.
