# Certified directional response between molecular reference blocks

The connected H4 square precision failure is closed on both exported stress fixtures. A second-order resolvent bound preserves directional coupling between small blocks and bounds only the remaining response. At hopping strength 1/1000, the rational energy interval improves from **1.2116563809476015e-5 to 2.5780947601473916e-10 Hartree** without changing P32 or its inherited reference partition. Selecting blocks directly from the perturbed matrix improves the width further to **1.7180947601473918e-10 Hartree**.

These are finite, exact certificates for exported rational 70-state Hamiltonians. They do not remove determinant enumeration, supply a general marginal-cone boundary, or establish efficient molecular scaling.

## The response inequality

Partition the actual Q block as H_QQ=C+R, with C block diagonal and R containing every interblock entry. Let W=H_QP and take the exactly computed symmetric row-sum bound eta>=||R||_2. Every reference block must satisfy C_j>=gamma_j I with gamma_j>b+eta, including blocks with zero direct P coupling.

Define

\[
A=C+R-bI,\qquad G=(C-bI)^{-1},\qquad
T=(C-(b+\eta)I)^{-1}.
\]

The exact resolvent identity and inverse order give

\[
A^{-1}=G-GRG+GRA^{-1}RG
\preceq G-GRG+GRT RG.
\]

Writing X=GW and Y=RX, the computed self-energy upper bound is

\[
\Sigma_1=W^TX-X^TRX+Y^TTY.
\]

Exact positivity of H_PP-bI-Sigma_1 proves H-bI positive. Unlike a scalar norm shift applied to the entire response, this expression retains the sign and direction of the first interblock correction. It is a sufficient bound, not an assertion that the truncated expression is the exact connected inverse or always dominates every other surrogate.

## Exact polynomial evaluation

Each small block's recurrence is discovered on all identity columns. This proves f_j(C_j)=0 on the entire block, including directions reached through R that were absent from the original W. The inverse at shift z is evaluated as -q_z(C_j)/f_j(z), where q_z(x)=(f_j(x)-f_j(z))/(x-z). Replay checks f_j(z) nonzero and multiplies (C_j-zI) by the proposed inverse to verify identity exactly.

The implementation never solves the full connected 38-state Q inverse in this response path. It still constructs the full 70-state Hamiltonian, performs a 32-dimensional retained positivity check, and uses a preexisting integer upper witness proposed from full numerical diagonalization. Polynomial evaluation does not erase those costs.

## Controlled results

All rows preserve the actual perturbed Hamiltonian, P32, and the independently evaluated integer Rayleigh upper witness for the corresponding strength.

| Hopping strength | Norm-shift response width (Ha) | Second order, inherited blocks | Second order, automatic blocks |
|---|---:|---:|---:|
| 1/1000 | 1.2116563809476015e-5 | 2.5780947601473916e-10 | 1.7180947601473918e-10 |
| 1/10000 | 1.027610466794303e-8 | 1.0010466794302913e-10 | 1.0010466794302913e-10 |

Inherited block dimensions and full recurrence degrees are 6,6,7,6,1,7,4,1. Automatic dimensions and degrees are 5,6,7,6,4,7,3. The automatic method starts from singleton Q states, processes actual edges by decreasing exact absolute weight with deterministic index tie breaks, and merges only when the combined size is at most seven. It has no partition-quality guarantee. Replay reconstructs this selection from the actual Hamiltonian and rejects a mismatched serialized partition.

Automatic interblock norm bounds are 1000000587/250000000000 and 100000587/250000000000, approximately 0.004000002348 and 0.000400002348. The inherited bounds were 0.007 and 0.0007. This closes dependence on a supplied original Q partition for these two tests; the retained P space still comes from the preceding experiment.

Certificates are approximately 50 KB. Independent standard-library replays take about 2.45–3.37 seconds in the recorded runs. Export, including proposal and replay, took about 7.7–10.3 seconds. Timings are measurements, not complexity bounds.

## Conditional higher-order error theorem

The implemented response is r=1 of a valid hierarchy. For any positive integer r, define

\[
U_r=\sum_{n=0}^{2r-1}(-GR)^nG+(GR)^rT(RG)^r.
\]

With K=G^(1/2) R G^(1/2), the odd finite geometric identity gives an exact positive remainder K^r(I+K)^(-1)K^r. Transforming back yields

\[
U_r-A^{-1}=(GR)^r(T-A^{-1})(RG)^r\succeq0.
\]

Let delta be any certified lower bound for lambda_min(C)-b with delta>eta. Since -eta I<=R<=eta I, inverse order and the spectral calculus of D=C-bI imply

\[
0\preceq T-A^{-1}
\preceq(D-\eta I)^{-1}-(D+\eta I)^{-1}
\preceq\frac{2\eta}{\delta^2-\eta^2}I.
\]

Consequently,

\[
0\preceq U_r-A^{-1}
\preceq\frac{2\eta}{\delta^2-\eta^2}
\left(\frac{\eta}{\delta}\right)^{2r}I.
\]

The self-energy error is bounded by the same scalar times any certified kappa>=||W||_2^2; the exact squared Frobenius norm is one available choice.

There is also an energy statement. Let u be a valid upper witness, delta_u=min_j gamma_j-u>eta, and

\[
\epsilon_r=\kappa\frac{2\eta}{\delta_u^2-\eta^2}
\left(\frac{\eta}{\delta_u}\right)^{2r}.
\]

This uniformly bounds self-energy error for every b<=u. For b<E0, minimizing the quadratic form of H-bI over Q gives S_exact(b)>=(E0-b)I in Euclidean P coordinates. Thus b=E0-epsilon_r-tau, for any tau>0, leaves the approximate Schur matrix at least tau I. The supremum of admissible lower endpoints is therefore at least E0-epsilon_r. The achievable interval width is at most the upper-witness error u-E0 plus epsilon_r, in the limiting supremum sense.

This proves geometric decay of the response conservatism inside the certified gap regime. It does not bound the upper-witness error, prove monotonicity of the numerical proposer, guarantee the greedy partition's gap, or control construction and rational bit costs as size or order grows. Only r=1 is currently implemented and exported in this module.

## Verification and remaining work

Four new exports have independent `python -S` replays under `results/marginal_molecular_gershgorin/square_connected_1_{1000,10000}_{second_order,automatic_second_order}`. Tests compare the bound to a direct inverse on a noncommuting rational example, recover exact response at R=0, reject uncoupled low states and insufficient shifted gaps, verify polynomial inverses, corrupt recurrences, patch direct solves to fail during replay, and reject false bounds and selection recipes.

The remaining central task is to discover and apply a compact retained space and a sufficiently gapped reference directly from an implicit many-body Hamiltonian, with controlled accuracy and resource growth. These finite experiments identify a useful response algebra; they do not solve that construction problem. General representability, the asymmetric quartic relaxation optimum, larger molecular transfer, and physical errors beyond the exported finite-basis electronic Hamiltonians remain separate open questions.
