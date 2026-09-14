# M10 blocker audit: what is and is not solved

The ten-mode result does not yet prove the optimum of the cubic SOS relaxation. `solve_coefficients` solves only the certificate (dual) side: it maximizes `b` subject to an exact coefficient identity represented in floating point, then exports a rounded rational certificate. The exact checker makes each exported lower bound sound, but it cannot certify that SCS found the best `b`.

The decisive next artifact is a feasible *moment functional* for the same finite dictionary. It is the SDP dual of the coefficient identity. If `y` assigns a value to every independent canonical coefficient, then feasibility requires the affine CAR contraction/number constraints and positive semidefiniteness of each moment matrix

\[
  M_B(y)_{uv}=y(v^\dagger u)\succeq0,
\]

for every Gram block, with `y(1)=1`. Its objective is `y(H)`. A rational interval enclosure of such a `y`, with Cholesky/LDL certificates for every block, proves a matching upper bound on the relaxation optimum. Matching it to the exported lower bound would close the M10 optimization blocker. In CVXPY the equality constraint's dual value is the natural floating-point seed; it must be exported and independently checked, rather than treated as a proof.

The observed failure of pair-hopping additions has a precise algebraic explanation for the tested antisymmetric construction. For any operator `A`,

\[
 (A-A^\dagger)^\dagger(A-A^\dagger)
 =A^\dagger A+AA^\dagger-A^\dagger A^\dagger-AA.
\]

If `A` is a pair annihilator, the last two terms are pair creation/annihilation quartics and may carry coherence. They do not vanish merely from nilpotency. They vanish only when the particular `A` is itself a single monomial (where `A^2=(A^\dagger)^2=0`) or when an explicitly verified CAR cancellation is present. Thus “pair hopping is diagonal” is not a general obstruction; it applies only to the single-transition operators actually tested. The correct experiment is to inspect the residual dual moment matrix's negative eigendirections and generate the corresponding linear combinations of pair operators, then verify whether their quartic square has a nonzero Hamiltonian projection.

The rank/kernel route is therefore concrete: obtain `y`, find a near-zero or negative eigendirection in each moment block, and add that operator family. If no such negative direction exists after rational certification, the current cubic relaxation optimum is proven and the remaining gap is entirely an upper-state issue. A theorem that all charge-neutral quartic additions are ineffective would require proving the relevant projection vanishes for every coefficient matrix; the existing data do not establish that.
