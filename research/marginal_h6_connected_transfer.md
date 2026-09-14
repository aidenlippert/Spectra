# Connected H6 transfer using reference factors

The spin-preserving hopping perturbation that joins the H6 complement into one 368-dimensional component now has a complete global energy certificate. At strengths1/1000,1/100, and1/50, the Schur interval widths are **2.0568984929976446e-10**, **1.5457570667278142e-10**, and **1.5408622850547893e-10 Hartree**. The proof retains the original 200- and168-dimensional factors and regenerates the physical upper and response on the connected Hamiltonian. It never factors the joined368-state block. The1/50 case requires the sharper square norm below. Subsequent [excitation-count certificates](marginal_spin_temple.md) tighten the last two intervals further.

This is transfer from the directly constructed H6 reference. It still evaluates explicit spin-zero configurations and does not establish scalable reference discovery or general chemistry.

## Reference gap, actual response

Let H0 be the exactly spin-symmetric reference Hamiltonian, Hs the spin-symmetric perturbed Hamiltonian, and H the actual perturbed rational export. The checker treats the two changes separately.

First it verifies complete Q coverage and fixed-point factors for the reference:

\[
QH_0Q\succeq\gamma_0Q.
\]

It then computes a rigorous bound eta on ||Hs-H0|| and requires

\[
\gamma\le\gamma_0-\eta.
\]

This gives QHsQ>=gamma Q for the same retained P. The current Hs may connect the reference blocks. Replay never claims that the reference block decomposition is also a decomposition of Hs.

The response matrices, residuals, and upper proposal are all generated using Hs. Only the complementary spectral threshold is inherited through the norm shift. The optimized residual Schur inequality can therefore recover much higher accuracy than an energy interval shifted directly by eta; eta does not impose an O(eta) floor on the final energy width.

Finally, exact coefficient comparison supplies epsilon>=||H-Hs||. The checked Hs lower is reduced by epsilon, while the upper is evaluated directly on H. For these perturbations epsilon remains5.4e-11 Ha.

Both H0 and Hs are independently validated as Hermitian, number conserving, and exactly SU(2)-invariant on the same mode/particle sector. The reference factor coverage is recomputed using the actual P. A supplied gap, norm, or factor receipt is never trusted without replay.

## The connected perturbation

For spatial orbitals0 and1, add

\[
\Delta H=t\sum_{\sigma=\alpha,\beta}
(a^\dagger_{0\sigma}a_{1\sigma}+a^\dagger_{1\sigma}a_{0\sigma}).
\]

This preserves spin symmetry but joins Q into one component of dimension368. The coefficient-norm gap bound is eta=4|t|. At t=1/1000, gamma0-eta lies about0.06416 Ha above the old evaluated upper, so the reference gap remains usable.

The constructor starts upper refinement from the reference's physical witness, then runs24 sparse Hamiltonian-action directions on the new Hs. Its resulting witness has400 nonzero amplitudes, compared with200 before the perturbation. It constructs a fresh response basis rather than assuming the old response remains adequate.

The strength1/1000 result uses28 response directions, reaches Hs width1.5168984929975216e-10 Ha, and transfers to original-H width2.0568984929976446e-10 Ha. The recorded transfer phase takes about341.81 seconds, excluding initial source-certificate verification and preparation. It includes repeated exact replays and is not a comparative performance benchmark.

The tenfold stronger case also passes using the coefficient-sum norm: eta=0.04 Ha and gamma=-6.30489910104 Ha. It uses31 response directions, reaches Hs width1.0057570667277666e-10 Ha, and transfers to original-H width1.5457570667278142e-10 Ha. Its recorded phase takes about410.22 seconds. The two final certificates occupy2,237,340 and2,283,057 bytes. Both use a connected368-state Q action graph and the same200/168 reference factors. These two H6 intervals use the coefficient-sum bound; the1/50 extension uses the sharper square certificate.

Final replay evaluates368 reference-H0 source states,400 current-Hs source states, and400 original-H source states, all drawn from400 distinct determinants. Source-certificate verification is accounted for separately:400 source-Hs and200 source-original states. These phase categories overlap; they must not be summed as a count of unique configurations. Repeated replays also repeat some operator/source evaluations.

## A sharper algebraic norm certificate

The coefficient sum is conservative for this hopping. For one spin, write K=a_p†a_q+a_q†a_p. The CAR identities imply

\[
I\pm K=\frac12\left[
(a_p\pm a_q)^\dagger(a_p\pm a_q)
+(a_p^\dagger\mp a_q^\dagger)^\dagger(a_p^\dagger\mp a_q^\dagger)
\right].
\]

Adding the two spins yields **four positive squares for each inequality**

\[
2|t|I\pm\Delta H\succeq0.
\]

The new norm checker accepts general bounded weighted CAR square recipes only after exact expansion proves both identities eta I +/- DeltaH = sum positive squares. Weights must be nonnegative, and no unverified residual is discarded.

For t=1/50, the exported bound is eta=1/25. A normalized physical vector proportional to amplitudes[1,-1,1,1] on determinants[243,246,249,252] attains Rayleigh value1/25 in the12-mode,6-electron sector. Thus the global Fock-space and fixed-N operator norms are exactly1/25. This witness intersects P and by itself does not certify the Q-compressed norm. A subsequent [four-state Q-supported witness](marginal_factor_width.md) on[627,630,633,636] proves the compressed norm also equals1/25 for the current P.

This sharper norm changes the usable gap range. At t=1/50, the coefficient-sum shift places the threshold0.011840474828199237 Ha below the old evaluated upper, failing that sufficient gate. The square-certified shift places it0.028159525171800762 Ha above the upper. The completed t=1/50 transfer uses31 fresh response directions and an independently refined24-direction upper witness. Its Hs width is1.0008622850545999e-10 Ha; original-H width is1.5408622850547893e-10 Ha. The certificate occupies2,288,690 bytes and its recorded transfer phase takes404.28 seconds. Exact replay again uses368 reference-H0,400 current-Hs, and400 original-H source states within400 distinct configurations. This resolves the gap-gate failure for this finite perturbation without a joined-Q factorization.

## Validation and remaining cost

Focused tests connect previously separate reference blocks, require the same-sector reference proof, reject incomplete coverage and optimistic gamma, and reject a lower endpoint that would pass on H0 but fails on the actual Hs. They also separate the original-H allowance from the reference gap allowance. A small transfer test prohibits new factor construction entirely, regenerates upper/response, and reaches the requested width using the old scalar factors.

The norm-square tests cover both signs of t, zero strength, negative weights, and mismatched identities. A small example fails the coefficient-sum gap gate but succeeds with the square-certified norm and independently replayed energy interval.

```sh
python -S -m experiments.marginal_spin_reduction \
  --verify results/marginal_h6/connected_spin/certificate.json

python -S -m experiments.marginal_spin_reduction \
  --verify results/marginal_h6/connected_spin_1_100/certificate.json

python -S -m experiments.marginal_operator_norm \
  --verify results/marginal_h6/hopping_norm_squares/certificate.json

python -S -m experiments.marginal_spin_reduction \
  --verify results/marginal_h6/connected_spin_squares_1_50/certificate.json
```

The joined-block factorization blocker is resolved on the accepted transfer. The larger unresolved issue is avoiding explicit growth of the spin-zero action support and discovering useful references and block decompositions as the system changes. A retained norm-shifted reference is conditional on a usable gap; it is not a universal solution to physical-marginal representability.
