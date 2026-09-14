# V9 approximation constructions

These procedures use established numerical ideas supplied to the experiment.
Their exact certificates do not imply learned procedures or hardware validity.

## Fixed-space Galerkin and whole-trajectory enrichment

For a retained Pauli coordinate set B containing the initial support, let P_B
be coordinate projection and define

    c0=O0,    c_(k+1)=P_B G c_k/(k+1),
    Otilde_m(t)=Σ_(k=0)^m c_k t^k.

Set q_k=(I−P_B)Gc_k. The residual coefficients are exactly

    r_k=−q_k     for k<m,
    r_m=−Gc_m.

In particular, omitting the last q_m from basis-selection scores loses a real
part of the residual. The checker is independent of these formulas: it
recomputes the entire polynomial derivative and generator action.

A valid, inexpensive l1 diagnostic separates the terminal inside component:

    E_l1 = Σ_(k=0)^m [T^(k+1)/(k+1)] ||q_k||_(Pauli,l1)
           + [T^(m+1)/(m+1)] ||P_B Gc_m||_(Pauli,l1).

This is an upper bound, not an infeasibility test. The proposer increases order
on one fixed B until E_l1 meets tolerance, the inside tail is small, or the
order budget is exhausted. It then tries the ordinary norm partitions. If the
certificate is still too wide, it scores each omitted label by its accumulated
absolute weighted coefficients, adds a bounded batch, and restarts. All failed
rounds and checks remain charged.

The BFS control uses the same recurrence, tail criterion and norm builder, but
adds complete breadth-first layers. Thus the comparison separates the policy
for choosing a basis from redundant witness reconstruction or a coarse degree
grid. Neither policy is learned in V9.

## Direct derivative collocation in a Krylov subspace

Build an untrusted numerical Krylov basis V and projected generator A, with
initial projected coordinate y0. Choose m nodes s_i in (0,1), let l_j be their
Lagrange polynomials, and define

    Q_ij = ∫_0^(s_i) l_j(s) ds.

Solve the stage equations

    U_i = y0 + T Σ_j Q_ij A U_j,

then form

    p(Ts) = y0 + T Σ_j [∫_0^s l_j(v)dv] A U_j.

This constructs a degree-m polynomial directly from the ODE. It does not
require m+1 independently computed matrix-exponential node values. In the
implementation, node calculations, orthogonalization and the dense stage
solve use floating point. The Pauli coefficients in normalized time are
rounded to a common rational grid, converted exactly into powers of physical
time, and exported with the exact requested initial observable.

None of those proposal calculations is trusted. The original checker verifies

    ||O0−Otilde(0)||∞ + ∫_0^T ||Otilde'−G Otilde||∞ ds <= eps

by exact rational Pauli arithmetic, including numerical basis/solve/rounding
errors. No eigenvalue-only bound for a nonnormal generator is substituted.
For projected dimension k, the direct stage system has dimension mk, not m;
its dense work and storage scale cubically and quadratically in mk. V9 caps
mk at384 and records the solve size.

The matched projected-Taylor arm uses the same Krylov builder and rational
grid. The fraction-free full Taylor baseline remains a separate, stronger
comparison when the Krylov setup and residual conversion cost more than they
save. The relevant error-control literature is recorded in
[polynomial_sources.md](polynomial_sources.md), including
[Hochbruck–Lubich](https://doi.org/10.1137/S0036142995280572) and
[Jawecki's defect-based study](https://arxiv.org/abs/2001.11922).

## What Bernstein subdivision would prove

For operator-valued Bernstein coefficients C_i on an interval of length h,
nonnegative basis functions and their exact integrals give

    ∫ ||R(t)||∞ dt <= h/(m+1) Σ_i U(C_i).

Dyadic de Casteljau subdivision constructs child coefficients as rational
convex combinations of parent coefficients. Exact norms give a nonincreasing
integrated coefficient envelope under refinement, by convexity and preservation
of the integrated weights. Independently chosen upper witnesses need not be
monotone; retain the parent certificate as an option. The scalar residual
1−2t has envelope1 before subdivision and1/2 after one split when integrated
coefficient weights are used. A maximum-per-interval bound remains1; it must
not be mistaken for the sharper integral envelope.

This extension has not been implemented as a V9 checker or counted as a
headroom win. It increases coefficient and norm-check work and leaves the
original generator replay cost intact.

## Incremental basis updates: an exact identity, limited measured opportunity

For B contained in B', the new Taylor coefficient c'_k equals the old c_k
through every prefix whose raw derivative has no component in B' minus B.
Writing delta_k=c'_k−c_k gives, while the old recurrence is available,

    delta_(k+1) = [P_B' G delta_k + (P_B'−P_B)G c_k]/(k+1).

Raw derivatives of an unchanged coefficient are unchanged. Their residual
classification can change at the first newly retained component, so a raw
cache hit does not automatically justify reusing a norm witness. A final
checker still recomputes the submitted residual independently.

The development diagnostic validates both prefix equality and
G c'_k = G c_k + G delta_k exactly on four actual trajectories. Prefix reuse
alone removes only about6–19% of generator multiply-adds. Difference updates
remove more generator multiply-adds but introduce enough map merging and
rational subtraction that no full-cost improvement follows from the counts.
The saved raw-derivative representations also consume memory. This is a
conventional incremental-computation opportunity, not an acquired m1 or a
causal reduction in m2 acquisition cost.
