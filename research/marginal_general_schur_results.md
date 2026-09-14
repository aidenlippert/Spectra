# Transfer after breaking pair charges: exact coupling-resolvent certificates

The Schur method now certifies two connected ten-mode Hamiltonians that break **all five individual pair-charge conservation laws** and **global left/right exchange**. One fixture also adds a pair-transfer interaction and diagonal perturbations. At strength 1/1000, the exact interval widths are **5.603353530955592e-7** and **5.996685860823592e-7**.

The new Hamiltonian need not preserve the reference symmetries. The method uses a known reference gap, computes the perturbation's action exactly, and automatically finds a short recurrence for that action under the reference Hamiltonian. Its accuracy remains perturbative around the fixed matched reference. This is not a solver for unrestricted chemistry with uniform accuracy.

## Measured transfer

All models have M=10,N=5 and use the same interaction units as the preceding matched tests. Start from H0 at hopping t=1/5. At strength s, add centered matched-edge perturbations s(i-2) and cross-pair hoppings

\[
-s\sum_{i=0}^4(1+i/5)
\left(a_i^\dagger a_{5+(i+1)\bmod5}+\mathrm{h.c.}\right).
\]

Together with the matched edges, these form a connected ten-site cycle. The mixed fixture additionally includes (s/2)(a0† a1† a7 a6 + h.c.), (s/7)n0 n6, and (s/11)n2. The latter terms make the perturbation within the retained symmetric space nonzero; the verifier handles it explicitly.

| Fixture / strength | Exact lower | Exact upper | Width |
|---|---:|---:|---:|
| Connected cycle, 1/1000 | 3.280992989097722 | 3.280993549433075 | **5.603353530955592e-7** |
| Mixed interaction, 1/1000 | 3.281080640428784 | 3.2810812400973703 | **5.996685860823592e-7** |
| Connected cycle, 1/100 | 3.2787389251822607 | 3.279700047047185 | 9.611218649243363e-4 |
| Mixed interaction, 1/100 | 3.279492976277002 | 3.280577503850644 | 1.0845275736414738e-3 |

The larger perturbation exposes the accuracy loss of this six-dimensional construction. A subsequent [enlarged coupling-space experiment](marginal_enlarged_schur_results.md) reduces these two larger-perturbation widths to 9.197e-8 and 8.325e-8 using 32 and 42 retained dimensions. Uniform accuracy as the norm bound approaches the reference gap remains unproved.

The scalar-gap baseline at strength 1/1000 gave widths 1.3327e-5 and 1.3438e-5. At strength 1/100 its widths were 0.01963 and 0.02580. Grouping Hermitian monomials improves the rigorous norm estimate; resolving the reference response improves it further. All baseline and refined artifacts remain in distinct directories under `results/marginal_general_schur/`.

## The exact sufficient condition

Let Z embed symmetric amplitudes f(k) into the canonical five-electron Fock basis. Its columns have 32 total nonzero entries, with fermionic signs from reordering pair-ordered creation operators into canonical mode order. The verifier checks exactly that

\[
D=Z^TZ=\operatorname{diag}\binom5k,
\qquad H_0Z=Z H_{\rm collective}.
\]

Set P=ZD^-1 Z^T and Q=I-P. The previously verified reference sector decomposition gives QH0Q >= cQ, where c is obtained from the d=1,2 Jacobi brackets. For the actual Hamiltonian H=H0+delta, compute

\[
A=D^{-1}Z^T\delta Z,
\qquad W=\delta Z-ZA=Q\delta Z.
\]

The calculation retains nonzero P delta P. It does not assume the perturbation preserves P, Q, pair charges, or exchange parity.

For a rigorous eta >= ||delta|| and z=b+eta<c,

\[
Q(H-bI)Q\succeq Q(H_0-zI)Q\succ0.
\]

Inverse order and the Schur complement then give the sufficient lower-bound test

\[
Z^THZ-bD-W^T[Q(H_0-zI)Q]^{-1}W\succ0.
\]

This is a six-dimensional rational positivity check. It treats the actual perturbation in the retained space and bounds its remaining effect through the reference complement.

## Automatically compressing the reference response

The code forms W,H0W,H0²W,... and finds an exact rational linear dependence

\[
f(H_0)W=0.
\]

It substitutes the discovered coefficients back into the original columns to verify the relation. No guessed spin labels for the perturbed Hamiltonian enter this step.

Writing g_z(t)=(f(t)-f(z))/(t-z), the response is

\[
W^T[Q(H_0-zI)Q]^{-1}W
=-\frac{1}{f(z)}\sum_k g_k(z)M_k,
\qquad M_k=W^TH_0^kW.
\]

The sign follows from (H0-zI)g_z(H0)W=-f(z)W. The verifier checks f(z) is nonzero. Exact LDL additionally verifies symmetry and positivity of the final matrix.

The cycle fixture has recurrence degree **four**:

```
f(t) = t^4 - 20 t^3 + (738/5)t^2 - (11892/25)t + 352809/625.
```

The mixed fixture has recurrence degree **six**. Both are discovered by exact elimination, rather than supplied as model-specific formulas. The small degrees reflect the few reference eigenvalues reached by the coupling. They do not imply a small recurrence for arbitrary reference Hamiltonians.

## Norm certificate

The basic estimate sums absolute canonical coefficients. The refined estimate groups each monomial with its Hermitian conjugate. A nondiagonal normal-ordered number-conserving monomial M is a partial isometry with M²=0; its domain and range are disjoint occupation subspaces. Consequently ||cM+c*M†|| <= |c|. Common creation/annihilation indices merely add occupation projectors. A diagonal monomial is a signed occupation projector and obeys the same coefficient bound.

After grouping, eta is 13s for the cycle fixture and (2115/154)s for the mixed fixture (exactly 423/30800 at s=1/1000). Tests exhaust all one- and two-body transition patterns on four modes and check that (M+M†)² is a diagonal zero/one projector. These are rigorous upper bounds, not measured operator norms.

## Replay, cost, and remaining limits

`experiments/marginal_general_schur.py` accepts a real, Hermitian, number-conserving degree-at-most-four Hamiltonian in its certificate. Replay reconstructs the reference embedding, norm bound, coupling recurrence, and positivity condition from that Hamiltonian. It independently evaluates the integer upper witness through exact CAR action. It also reports which pair charges are broken and whether global exchange remains invariant.

All nine saved baseline/refined artifacts passed a separate `python -S` replay. The four final certificates have directories ending in `_paired_resolvent`. For example:

```sh
python -S -m experiments.marginal_general_schur --verify results/marginal_general_schur/mixed_1_1000_paired_resolvent/certificate.json
```

The refined generation runs took approximately 0.22–0.35 seconds locally, including numerical proposal of the upper witness and exact replay. These are individual observations, not controlled scaling benchmarks. Lower construction explicitly enumerates 32 reference configurations and the configurations reached by their couplings; the upper witness enumerates the full 252-state fixed sector. The coupling reaches 104 basis states for the cycle and 118 for the mixed fixture.

This removes preservation of pair charges and exchange symmetry as requirements on the tested new Hamiltonians. It retains a specially structured reference, a certified reference gap, and a perturbation norm estimate. General reference discovery, scalable construction without Fock-state enumeration, larger chemical systems, and uniform accuracy at larger perturbations remain unresolved. None of these bounds determines the optimum of the quartic marginal relaxation.
