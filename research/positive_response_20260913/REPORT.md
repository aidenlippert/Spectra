**A molecular sector gap, and the limit of a common-zero response**

This pass supplies one concrete prerequisite for collective elimination on
the frozen H6 Hamiltonian. A sector containing 210 six-electron states has
an independently certified positive gap at the chosen energy target.
The proof uses a small one-body bound and the existing collective tail
certificate. It enumerates no states and computes no many-body inverse.
The response into that sector, and the lower bound on the remaining
effective Hamiltonian, are still open. No new full molecular ground-energy
interval is certified in this pass.

The second result explains a limitation of the synthetic response controls:
the common zero of the retained molecular density squares lies far above
the molecular ground energy. The third result sharpens the previous
consistency defect into a spin-rotation-invariant T1 separator.

All molecular statements concern the frozen rational finite-basis
Hamiltonians. The
[protocol](/Users/aidenlippert/Documents/Spectra/research/positive_response_20260913/PROTOCOL.md)
and
[derivations](/Users/aidenlippert/Documents/Spectra/research/positive_response_20260913/DERIVATION.md)
give the identities and scope of the tests.

**The usable molecular gap.** Let

\[
Q=n_{5\uparrow}n_{5\downarrow},\qquad P=I-Q,
\]

with zero-based spatial labels and interleaved spin orbitals, so Q fixes
spin orbitals 10 and 11 occupied. On the full N=6 sector, Q has dimension
\(\binom{10}{4}=210\); P has dimension 714. These are dimension counts,
not generated bases. Use the frozen physical upper
\(U\approx-6.333058626233001\) Ha and fix

\[
b=U-0.0016\ \mathrm{Ha}\approx-6.334658626233001\ \mathrm{Ha}.
\]

The exported certificate proves

\[
Q(H-bI)Q\succeq\delta Q,\qquad
\delta=
\frac{1912741854356760239502987364891}
{10000000000139879611500000000000}\ \mathrm{Ha}
\approx0.1912741854330005\ \mathrm{Ha}.
\]

Here is the entire bounding mechanism. Write the retained Hamiltonian as
\(H_{\rm ret}=c+\mathrm d\Gamma(t)+R\), where
\(R=\frac12\sum_i\lambda_i Q(L_i)^2\succeq0\), and reuse the exact
lower tail shift \(\ell=-3.632\times10^{-9}\) Ha. Compressing to Q leaves
\(c+2t_{55}+\mathrm d\Gamma(t_{\ne5})\), with four remaining electrons.
The sum of its four lowest diagonal spin-orbital energies, minus four
times the maximum absolute off-diagonal spatial row sum, bounds it below.
Adding \(\ell\) gives

\[
QHQ\succeq-\frac{7679230551}{1250000000}\,Q
=-6.1433844408\,Q\quad\text{(Ha)}.
\]

Subtracting b gives the displayed gap. The exact checker also verifies
spin independence of t and replays the tail certificate. This proof uses
only \(R\succeq0\); the quantitative Lie-bracket bound below is unnecessary.

The
[standalone certificate](/Users/aidenlippert/Documents/Spectra/results/positive_response_20260913/sector_gap_certificate.json)
is **430 bytes**, bound to **74,534 bytes** of existing fixture and tail
data. Its
[independent replay](/Users/aidenlippert/Documents/Spectra/results/positive_response_20260913/sector_gap_replay.json)
took 0.305 seconds including data loading in the recorded run. It imports
no numerical packages and needs no Lie-bracket reconstruction. The
certificate proves a gap for D; it does not assert \(H\succeq bI\).

All six analogous double-occupancy choices were checked. Their sufficient
lower bounds for shifted D, using \(R\succeq0\), are:

| Fixed doubly occupied spatial orbital | Lower bound on D (Ha) | Strict gap established? |
|---:|---:|:---|
| 0 | −0.663965133 | No |
| 1 | −0.697295367 | No |
| 2 | −0.663965133 | No |
| 3 | −0.438765006 | No |
| 4 | −0.087392613 | No |
| 5 | +0.191274185 | Yes |

The five negative sufficient bounds are inconclusive about their actual
gaps. No eliminated-sector eigenvalues were computed.

**Why a shared zero misses the molecular low energy.** For N=s electrons
in s spatial orbitals, the existing exact orbital/spin Casimir identity is

\[
C_{\rm orb}=\sum_{pq}(E_{pq}-\delta_{pq}N/s)^\dagger
                         (E_{pq}-\delta_{pq}N/s)
=2[S_{\max}(S_{\max}+1)-S^2],\quad S_{\max}=s/2.
\]

The retained patterns generate the full traceless spatial Lie algebra:
dimension 35 for H6 and 63 for H8. Thus a vector killed by every retained
density square is killed by all traceless spatial densities. The Casimir
identity identifies their common kernel exactly as the maximum-spin
multiplet, with dimension s+1. Its retained energy is
\(c+\operatorname{tr}t\). The original Hamiltonian's energy on this kernel
is bounded below by that value plus the certified lower tail shift.

| Fixture | Patterns | Common-kernel dimension | Retained kernel energy (Ha) | Physical upper (Ha) | Kernel lower minus physical upper (Ha) |
|---|---:|---:|---:|---:|---:|
| H6 | 10 | 7 | −5.775989062 | −6.333058626 | +0.557069561 |
| H8 | 14 | 9 | −8.487857873 | −9.252151640 | +0.764293759 |

Consequently, the ground state cannot lie in the common kernel of these
retained squares. A molecular construction must handle their nonzero
interaction energy and correlated fluctuations. This result leaves open
shifted squares, different positive decompositions, and compact responses
with a frustrated bulk. Full Lie closure alone does not rule out a short
noncommutative operator program.

The Casimir identity was already derived in the project's
[collective-tail proof](/Users/aidenlippert/Documents/Spectra/research/molecular_collective_20260913/PROOF.md).
This pass verifies it again by exact CAR expansion and applies it to the
kernel and energy question. H8 is an existing fixture held out from the
new H6 calculation. Its result tests the kernel diagnosis; it supplies no
new molecular response or T1 transfer result.

We also tried an explicit quantitative interaction gap. Normalized Lie
brackets propagate bounds
\(\|Q(A_j)\psi\|^2\le k_j\langle\psi,R\psi\rangle\).
Exact reconstruction of a traceless matrix basis then gives
\(R\succeq C_{\rm orb}/K\). On a specified doubly occupied orbital,
\(C_{\rm orb}\succeq12I\) for H6. The accepted value was

\[
K\approx4.44348451\times10^{16}\ \mathrm{Ha}^{-1},\qquad
12/K\approx2.70058329\times10^{-16}\ \mathrm{Ha}.
\]

This positive gap is quantitatively useless here. Construction took
0.601 seconds and final independent replay 0.369 seconds. The certificate
contains a 35-dimensional reconstruction, reaches 706-bit reconstruction
fractions, and occupies 215,118 bytes. Those costs and the weak constant
are preserved in the
[gap receipt](/Users/aidenlippert/Documents/Spectra/results/positive_response_20260913/gap_replay.json).
No additional optimization of this bound was attempted.

**A nonzero response error can remain in a positive block.** For an
already certified reference response X0 with zero reference residual,
\(K_0\succeq0\), and shifted reference block \(D_{\rm ref}\succ0\), adding positive diagonal blocks
V0 and V1 gives the exact congruence

\[
T^\dagger(H-bI)T
=\operatorname{diag}(K_0,D_{\rm ref})+\operatorname{diag}(V_0,0)
+[-X_0,I]^\dagger V_1[-X_0,I]\succeq0,
\quad T=\begin{pmatrix}I&0\\-X_0&I\end{pmatrix}.
\]

The actual residual is \(E=-V_1X_0\). It can be large and noncommuting
while the full displayed block remains positive. These are ordinary
congruence and positivity identities; rigorous elimination through Schur
maps is established mathematics.
[Dusson, Sigal and Stamm](https://arxiv.org/abs/2105.02058).

The concrete control adds the periodic conditional-exchange bath
\(L=\sum_i n_{i+2}(I-\mathrm{Swap}_{i,i+1})\) to the two sectors of the
previous homogeneous central-spin model, with nonnegative strengths
\(\lambda_0,\lambda_1\). Each term is a positive square and kills a
symmetric Dicke vector. The actual D is nonscalar; its minimum is still
the reference gap d. All three small controls prove \([L,J_-]\ne0\),
\(E\ne0\), and failure of \(K\succeq E^\dagger E/d\) on an explicit
retained vector. The full congruence succeeds.

| Bath spins; excitations | Added strengths λ0; λ1 | Enumerated control states | Negative expectation of K−E†E/d | Exact ground-interval width |
|---|---|---:|---:|---:|
| 5; 3 | 1/3; 10 | 20 | −49.5413 | 1.915×10⁻¹⁵ |
| 7; 3 | 2; 17 | 56 | −231.080 | 1.040×10⁻¹⁴ |
| 9; 4 | 5/7; 100 | 210 | −1982.49 | 1.136×10⁻¹⁴ |

These are synthetic model energy units. Separate numerical diagonalizations
agree with the analytic intervals within a 10⁻¹⁰ diagnostic tolerance;
numerics accept no certificate. Two further formula-specified cases with
64 and 1024 bath spins use no state enumeration. Their certificate inputs
occupy 370 and 383 bytes; their widths are 3.10×10⁻¹⁴ and 9.12×10⁻¹³.
This restricted model has a shared zero for its added interactions.
The inequality \(H\succeq H_0\) already explains its lower bound. The
control tests retention of a positive residual block when the separate
error estimate fails.

The project already contains noncommuting frustration-free controls in
[earlier structural routes](/Users/aidenlippert/Documents/Spectra/research/certificate_scaling/low_rank_structural_routes.md).
Conditional-exchange frustration-free spin models also have an established
literature, including the Fredkin chain; the precise periodic control here
is not identified with that model.
[Adhikari and Beach](https://arxiv.org/abs/1805.00532).

**The missing consistency condition survives spin completion.** Starting
from the frozen 12-term three-removal operator C, apply the adjoint spin
Casimir without fitting any coefficients:

\[
C_{1/2}=\frac{15C/4-J_{\rm ad}^2(C)}{3},\qquad
C_+=[S_+,C_{1/2}].
\]

The projected operator still has 12 terms on the same four spatial
orbitals. Its coefficients, over denominator 300 and in the frozen triple
order, are

```
-9, -3, -18, 61, 115, -9, -176, 60, -14, -128, 39, 142
```

The exact spin-one-half covariance identities hold. Complete the doublet:

\[
P_{\rm spin}=\{C_{1/2}^\dagger,C_{1/2}\}
             +\{C_+^\dagger,C_+\}\succeq0.
\]

Exact CAR expansion proves that P_spin commutes with all three total-spin
generators. Sixth-degree terms cancel, leaving 181 constant, quadratic
and quartic terms. Every used moment is explicitly present in the frozen
functional; no higher moments or unrecorded values are supplied. Its value is

\[
y(P_{\rm spin})=
-\frac{2831239196789434247458020091}
{825000000000000000000000000000}
\approx-0.00343180508701750.
\]

This is a **dimensionless positivity violation**, not an energy improvement.
The previous full-family feasibility audit remains hash-bound and accepted,
so this spin-scalar T1 condition still excludes its witness. T1 is an
established representability condition; the instance projection and
diagnosis are the results of this pass.
[Mazziotti](https://arxiv.org/abs/1207.0541).

The negative value is independent of the spin axis. Mixing the functional
with the uniform N=6 trace still violates this condition for trace fraction
\(0\le t<0.00932030038220734\); the exact threshold is recorded in the
[molecular replay](/Users/aidenlippert/Documents/Spectra/results/positive_response_20260913/molecular_replay.json).
Those mixtures inherit old-family feasibility by convexity. No old-family
feasibility claim is made for arbitrary spin-rotated functionals.

Independent fermionic action shows both doublet components and their
adjoints annihilate the entire maximum-spin multiplet. Generating its seven
vectors uses 64 amplitudes. A closed-shell singlet gives the positive
expectation 1/500, equal to four squared norms. Thus the separator probes
lower-spin correlations absent from the common-zero bulk. Its coefficients
still originate in the earlier fitted separator; this is not independent
molecular transfer.

**Verification and accounting.** All accepting replays use Python -S and
exact rational arithmetic, with no NumPy, SciPy, CVXPY or PySCF imports.
Twenty focused tests pass, including altered reconstruction coefficients,
false norm bounds, forward bracket references, nonpositive bulk strengths,
wrong input bindings, overstated sector gaps and altered spin projections.
The original checker remains unchanged and rejects the new model schema.
All 718 files in the preceding manifest match their frozen hashes.

The final combined molecular/spin replay took 4.788 seconds. The
[cost ledger](/Users/aidenlippert/Documents/Spectra/results/positive_response_20260913/cost_ledger.json)
records 25.125 seconds of summed measured stage times, including repeated
development checks. This sum is not calendar end-to-end time. Untimed
editing, preliminary inspections, import overhead outside module timers,
and an initial diagnostic refusal are disclosed. That refusal came from
using diagonal spread alone to diagnose nonscalarity: the M=5 block has a
constant diagonal and nonzero off-diagonal entries. The action calculation
exposed the issue; the corrected diagnostic and a focused regression test
cover it.

Each complete exact response control run generates 286 basis labels;
the separate numerical controls generate another 286 and store 47,636
matrix entries. Repetitions and test controls are charged in the ledger.
The molecular energy comparisons replay inherited FCI-derived upper
witnesses with 200 H6 and 1000 H8 amplitudes, totaling 3,095,600 word-state
checks per complete molecular run. Their earlier FCI discovery cost is
inherited, not removed. The old full-family audit is reused from its
frozen receipt. No new molecular SDP, FCI, or fixture generation occurred.

From the workspace root, the independent replays are:

```sh
/opt/homebrew/Caskroom/miniconda/base/bin/python -S -m research.positive_response_20260913.sector_gap
/opt/homebrew/Caskroom/miniconda/base/bin/python -S -m research.positive_response_20260913.coercivity --replay
/opt/homebrew/Caskroom/miniconda/base/bin/python -S -m research.positive_response_20260913.block_response --replay
/opt/homebrew/Caskroom/miniconda/base/bin/python -S -m research.positive_response_20260913.molecular_diagnostic --replay
/opt/homebrew/Caskroom/miniconda/base/bin/python -S -m unittest research.positive_response_20260913.test_rules -v
```

**The next mathematical target is now specific.** Use
\(Q=n_{5\uparrow}n_{5\downarrow}\), its independently certified gap, and
\(B=PHQ\) to construct a short response X. Derive
\(E=B^\dagger-DX\) and the retained remainder as operator programs before
expansion. Test whether their positivity can be certified while retaining
the lower-spin structure identified above. Charge program size, expanded
verification work and error together. A small D-gap certificate supplies
the first gate; useful compression still requires a compact response and
a controlled retained remainder.
