# Reviewed local construction and collective correction

All energies below use Hubbard hopping t=1, not molecular hartrees. The model
is the fully coupled open 2-by-4 repulsive Hubbard ladder, U=8, with four up
and four down electrons. No claim of a general many-body solution or novelty
of projection/SDP inequalities is made.

## Exact overlapping representation

Number the four two-site rungs A,B,C,D. Let

\[
 h_r=T_{r,\mathrm{vertical}}+8(d_{r,1}+d_{r,2})-4(N_r-2).
\]

Let v_rs contain both leg hoppings between neighboring rungs. Define

\[
 h_0=h_A+\tfrac12h_B+v_{AB},\quad
 h_1=\tfrac12h_B+\tfrac12h_C+v_{BC},\quad
 h_2=\tfrac12h_C+h_D+v_{CD}.
\]

Their sum is H-4(N-8), hence equals H in the declared global charge sector.
For arbitrary real symmetric, spin-number-preserving one-rung matrices J_0,J_1,
set

\[
 A_0=h_0+J_0^{B}-\ell_0 I,\quad
 A_1=h_1-J_0^{B}+J_1^{C}-\ell_1 I,\quad
 A_2=h_2-J_1^{C}-\ell_2 I.
\]

If all A_i are positive semidefinite on **every** local occupation sector,
then E_0 >= L=sum ell_i. The two boundary operators cancel exactly. No local
particle count is fixed. Each two-rung Fock space has 256 states, decomposed
into 25 (N_up,N_down) blocks, with maximum dimension 36. The construction and
checker never form the 4,900-state global balanced sector.

The constructor writes this exact local model directly to a conic solver.
Solver outputs propose J and ell; rationalization followed by exact
fraction-free LDL decides acceptance. The public cone convention is
[Clarabel's scaled upper triangle](https://clarabel.org/stable/api_cone_types/).

## An exact ceiling for the entire boundary-operator family

Suppose positive local density matrices rho_0,rho_1,rho_2 have trace one and
agree on each shared rung:

\[
 \operatorname{Tr}_A\rho_0=\operatorname{Tr}_C\rho_1,\qquad
 \operatorname{Tr}_B\rho_1=\operatorname{Tr}_D\rho_2.
\]

For every accepted local certificate of the above form,

\[
 0\le\sum_i\operatorname{Tr}(\rho_i A_i)
   =\sum_i\operatorname{Tr}(\rho_i h_i)-\sum_i\ell_i.
\]

Thus D=sum Tr(rho_i h_i) is an upper ceiling on the lower bounds attainable
by that complete boundary family. The matrices need not extend to a physical
global state, so D is **not** a physical ground-energy upper bound.

The exact dual repair first fixes trace by adding (1-trace)I/256. Shared
one-rung marginals are replaced by their average. If the required changes
in a patch's left and right reductions are delta_L,delta_R (both traceless),
add delta_L tensor I/16 + I/16 tensor delta_R. This enforces overlap equations
exactly without changing the trace. Finally mix every patch with I/256 using
a common rational weight. Consistency survives. Exact PSD checks, not the
size of the mixture or a numerical eigenvalue, decide acceptance.

This is an established local-marginal/SDP mechanism; see
[Eisert (2026)](https://doi.org/10.1103/cz6k-y46r) and
[Fawzi, Fawzi, Scalet (2024)](https://arxiv.org/abs/2311.18706).
The new evidence here is the replayed bound for this fixed interacting case,
not discovery of the general SDP principle.

## Carrying a collective compatibility penalty

Choose integer local vectors u_i of definite spin populations, with nonzero
norms n_i. Put P_i=u_i u_i^T/n_i on its patch. The checker proves

\[
 A_i\succeq g_i(I-P_i),\qquad g_i>0,
\]

by exact local PSD tests. No vector is assumed to be an exact local eigenvector.
Its rational projector and the stated inequality are checked directly.

The outer patches AB and CD are disjoint and their projectors are even
fermionic operators. Hence P_0 and P_2 commute. With g=min(g_0,g_2), h=g_1,

\[
 A_0+A_1+A_2\succeq g(I-P_0P_2)+h(I-P_1).
\]

The product projector P=P_0P_2 is rank one on the whole four-rung space, but
is represented as a product of two local vectors. Its squared overlap with
the middle projector is computable locally:

\[
 p=\|PP_1\|^2
 =\frac{\langle u_1|\rho_B^{(u_0)}\otimes\rho_C^{(u_2)}|u_1\rangle}
 {n_0n_1n_2},
\]

where the reductions on the numerator use unnormalized integer vectors.
Only 16-by-16 one-rung reductions and local vector contractions are needed.

The standard two-projection bound is

\[
 g(I-P)+h(I-Q)\succeq
 \frac{g+h-\sqrt{(g-h)^2+4gh\|PQ\|^2}}2 I.
\]

One proof uses the synthesis operator W=[sqrt(g)V_P,sqrt(h)V_Q] for
isometries into the two projector ranges. Its Gram has diagonal blocks gI,hI
and off-diagonal norm at most sqrt(gh)||PQ||. A 2-by-2 scalar quadratic form
bounds its largest eigenvalue; subtract from (g+h)I.

The checker does not need an irrational square root. It accepts rational gamma
only after verifying

\[
 0\le\gamma\le\min(g,h),\qquad
 (g-\gamma)(h-\gamma)\ge ghp.
\]

The full ground-energy lower bound is then L+gamma. This inequality carries
information about incompatible overlapping local subspaces that is absent from
a sum of separately positive patch matrices. Its proof has global validity,
but its accepting data and operations are local.

## Actual accepted values

The zero-boundary lower is -4.4925060891t. Optimized local boundary operators
raise it to

\[
 L=-4.375173654t.
\]

The independently repaired local-marginal dual gives

\[
 D=-4.3751532484263764t.
\]

Therefore the optimum of the declared local family is enclosed in a width of
0.0000204055736236t. More optimization of this same family cannot remove its
large physical error.

The exact local spectral profiles have gaps

\[
 (g_0,g_1,g_2)=(0.30865585,\ 0.5208039841,\ 0.30865585)t.
\]

They give p approximately 0.064278254 and accept

\[
 \gamma=0.2678136631t,\qquad L+\gamma=-4.1073599909t.
\]

In particular,

\[
 E_0-D\ge0.2677932575263764t.
\]

This is an exact separation from the old local family, independent of a global
ground-state reference. It is not merely a better solution within that family.
The exact rationals and both independent witnesses are in the result artifacts.

## Exact information loss in the single-gap compression

Consider the distinct coarsened Hamiltonian

\[
 F=\sum_{i=0}^2 g_i(I-P_i).
\]

It satisfies all three retained profile inequalities with equality. Therefore
any larger universal correction inferred from those inequalities alone must
also lower-bound F. Choose the normalized state
|a>_A tensor u_1/sqrt(n_1) tensor |d>_D, with the two local occupation labels
chosen to give global N_up=N_down=4. Its middle penalty is zero. The other
two expectations are exactly

\[
 q_0(a)=\frac{\sum_c(\sum_b u_0(a,b)u_1(b,c))^2}{n_0n_1},\qquad
 q_2(d)=\frac{\sum_b(\sum_c u_1(b,c)u_2(c,d))^2}{n_1n_2}.
\]

Thus

\[
 \lambda_{\min}(F)\le g_0(1-q_0(a))+g_2(1-q_2(d)).
\]

The exact locally enumerated choice a=d=6 gives 0.5678716126607801t. Even an
optimal global analysis using only these fixed coarse profiles therefore
cannot guarantee a correction greater than that number. With the frozen
physical upper and baseline L, this information alone leaves at least
0.7813792356485363t of interval width. This is an information limit across
operators satisfying the same inequalities; it does not bound a method that
uses the original Hamiltonian's additional excitation structure. The trial
energy here is an upper for F, **not** an upper for the physical H.

## Limits and costs

The corrected lower is still more than 1t below the independently certified
variational upper near -3.025922806t. The 0.001t diagnostic interval target is
not met. This is a Hubbard-model test, not a molecular accuracy milestone.

A numerical screen retaining ranks 1,2,4,6,8,12,16 in this same grouped-projector
rule found no gain beyond rank one. It is a bounded search, not an obstruction
to all larger ranks or richer use of the local spectra. Higher-rank overlap
Gram dimensions grow as r^2 and contraction entries as 256 r^3 here; at r=16
the numerical overlap contraction already contains 1,048,576 entries.

The accepting lower checker enumerates 256 local occupations, not the global
sector. A separate regression reconstructs the Hamiltonian on all 4,900 global
configurations to catch errors in model binding. A separately charged global
reference supplies a physical upper. Those uses of enumeration are explicitly
not part of a claimed enumeration-free complete pipeline.

This pass does not construct X,K,R for the earlier relative-residual theorem,
or prove that induced operators stay compact after repeated reductions. It
supplies a different, exactly checked local-to-global correction and diagnoses
its current failure quantitatively. Neither that improvement nor a small
fixed-size calculation establishes a field-level many-body breakthrough.
