# Charge-sector spin polynomials and exact projected positivity

The H6 complement proof now shares metric work across spin occupations. The accepted proof covers140 neutral ionic charge patterns using112 coefficient-count inequalities and28 exact3×3 projected positivity checks. It certifies **QHQ≥−6.264Q Ha** without spin-endpoint evaluation or determinant actions in construction/replay. Explicit charge-pattern enumeration remains.

## Freeze charge, retain spin algebra

At each spatial site let q_i=n_iα+n_iβ−1. Fix a neutral ionic charge pattern q∈{−1,0,1}^m with at least one doublon. At a doublon both occupations are1, at a holon both are0, and at a singly occupied site

\[
n_{i\alpha}=x_i,\quad n_{i\beta}=1-x_i,\quad x_i^2=x_i,
\qquad\sum_i x_i=r/2,
\]

where r counts the single sites. Each source event and coherent amplitude is now a Boolean spin polynomial. The target-valence exclusion is determined by q+Δ, and the positive metric ratio is a rational scalar shared across every spin occupation with this q. The compiler verifies the sign of each active grouped amplitude over its allowed source event before replacing its absolute value. Sign-changing amplitudes are rejected; they are not silently linearized.

Metric ratios are cached by Δ within each charge pattern. The complete row polynomial is formed before a lower bound is applied. Squarefree monomials with degree greater than r/2 vanish on this Boolean slice and are discarded. This identity would not hold on an unrestricted continuous hyperplane.

H6 charge accounting is exact:

| Doublons | Charge patterns | Spin occupations per pattern | Total occupations |
|---|---:|---:|---:|
| 1 | 30 | 6 | 180 |
| 2 | 90 | 2 | 180 |
| 3 | 20 | 1 | 20 |
| Total | 140 | — | 380 |

The20 last patterns fully determine the occupations. Zero spin-endpoint calls does not mean zero enumeration of all kinds.

## The first bound diagnoses the remaining difficulty

For each degree d, precisely binomial(r/2,d) squarefree monomials are active. Summing that many smallest coefficients, including implicit zeros, gives a rigorous lower bound. This certifies112 charge patterns. Its worst bound is **−6.651783167300551 Ha**, and28 one-doublon patterns fail the desired threshold. The failed bound is retained as `cardinality_probe.json`.

All failing patterns have four single sites. Their spin slice has exactly two alpha occupations, so every remaining polynomial is quadratic. The next gate uses that quadratic structure directly.

## A projected quadratic certificate

Write

\[
p(x)=p_0+\sum_i p_i x_i+\sum_{i<j}p_{ij}x_ix_j.
\]

Set x_i=(1+z_i)/2, giving p=c+lᵀz+zᵀMz with

\[
c=p_0+\tfrac12\sum_i p_i+\tfrac14\sum_{i<j}p_{ij},\quad
l_i=\tfrac12p_i+\tfrac14\sum_{j\ne i}p_{ij},\quad M_{ij}=p_{ij}/8.
\]

M has zero diagonal. The gate checks that every l_i is equal, so lᵀz=0 on the balanced subspace. It does not assume this coefficient identity from approximate symmetry.

For r=4 use the mutually orthogonal balanced vectors

\[
v_1=(1,1,-1,-1),\quad v_2=(1,-1,1,-1),\quad v_3=(1,-1,-1,1)
\]

and the diagonal correction

\[
d=-\frac14\sum_{a<b}(v_a^TMv_b)(v_a\odot v_b).
\]

Its entries sum to zero. Thus zᵀdiag(d)z=0 on Boolean z_i²=1, while M+diag(d) is diagonal in this balanced basis. For a requested bound γ, set λ=(γ−c)/r. With B having columns e_i−e_last, the verifier checks by exact rational LDL that

\[
B^T\big(M+\operatorname{diag}(d)-\lambda I\big)B\succ0.
\]

This bounds the diagonally corrected quadratic on the balanced sphere. The diagonal correction vanishes on physical Boolean spins, proving p≥γ for every physical spin assignment. The H6 matrices have dimension3. The checker uses strict positive pivots; a boundary case with zero pivots may be refused even when a semidefinite certificate exists.

For other even r, the implementation uses the zero diagonal correction and the same projected PSD test as a sufficient condition. Higher-degree polynomials and unequal linear coefficients are refused by this gate. The special four-spin exactness does not establish an exact solver for arbitrary spin slices.

## Arithmetic cost and independent checks

The proof compiles620 polynomial monomials across140 charge patterns. Its shared metric calculation uses **178164 scalar exponentiations**, compared with **862272 metric scalar endpoint evaluations** in the preceding simplex proof. This compares that counted operation only; it is not a total-runtime speedup claim. The new compiler also processes35550 transition groups, caches8910 charge-change kernels, and performs28 small PSD checks.

An independent validation computes all380 physical CAR rows and compares each exact rational row with its compiled polynomial. Every row agrees. This exhaustive check is labeled validation-only; the proof generator and standard-library replay do not use it. The same receipt records the28 PSD gates, their exact pivots and diagonal corrections, and independent row minima.

Focused tests also cover density-assisted and pair hopping, fixed-sign rejection, neutral charge accounting, coefficient bounds, four-spin extrema, a larger sufficient PSD gate, malformed supports, false thresholds, incomplete budgets, and energy Hamiltonian/reference binding.

## Direct transfer to the perturbed Hamiltonian

The unchanged unperturbed metric also certifies the actual perturbed spin-symmetric Hamiltonian directly. A first direct check passes at−6.304 Ha; a stronger check passes at **−6.264 Ha**, recovering the0.04 Ha sacrificed by the previous perturbation-norm argument. Neither direct check imports a reference gap or hopping norm.

This is the previously studied strength1/50 hopping between the first two localized orbitals. It establishes direct transfer to that changed Hamiltonian, not to a new molecule or size. The two inherited energy intervals also remain valid after substituting the charge-spin gap proof. Fresh response discovery with the stronger direct gap and excitation offset1/100 uses **17 directions**, down from18. The independently replayed nearby spin-symmetric interval has width **3.331192213308058e−14 Ha**. After the checked Hamiltonian approximation error is included, the input-Hamiltonian interval has width **3.9293895949815965e−12 Ha**, down from7.227269910284279e−12 Ha. The upper witness is inherited and rechecked; the response is newly discovered. Full-sector response work remains.

## Unresolved work

The spin dimension is certified collectively here, but charge patterns still number

\[
\sum_{d=1}^{m/2}\binom{m}{d}\binom{m-d}{d}.
\]

This grows combinatorially. Further progress requires a charge-level inequality or elimination scheme that shares work between patterns, while preserving the metric products and ionic exclusion. Larger spin slices may retain higher-degree terms or fail the sufficient PSD test. Sign stability can also fail for other Hamiltonians.

The full energy witness and response still use the finite spin sector. Reference growth, larger-system transfer, the full marginal cone, and finite-temperature/dynamic/synthesis predictions remain unproved. The present result is an exact finite operator certificate with shared spin algebra, not a general chemistry compiler.


Verification: **332 marginal regression tests passed in273.719 seconds**. Six new standard-library gap/energy replays match every saved replay field. Independent exhaustive validation matches all380 physical Q rows exactly. Detailed receipts are recorded in `results/marginal_h6/charge_spin/progress.json`.


Follow-up: [Shared charge-tail obstruction](marginal_charge_tail.md) exactly brackets a failed uniform-tail envelope and excludes two simpler metric families at the required threshold. The accepted charge-spin certificates remain valid.
