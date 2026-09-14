# v9 polynomial and truncation baselines: bounded source review

## Scope and correction to v8

The v8 proposal is best treated as a proposal generator. Rational nodes are a
valid implementation convenience; an algebraic certificate for the usual
Chebyshev nodes is not needed if the exported rational polynomial is checked by
the exact residual checker. Likewise, independent enclosures of every node
value are not a validity requirement when the checker recomputes the
coefficients/residual from the submitted polynomial. The construction should
collocate the ODE derivative, rather than first paying for `m+1` expensive
exponential node values and interpolating them.

## Sound collocation certificate

Let the Heisenberg ODE be `Y'(t)=G Y(t)`, `Y(0)=O`, and let `P` be a continuous
degree-`m` polynomial on `[0,T]`. Define the exact residual

```text
R(t) = P'(t) - G P(t).
```

For a positive unital Lindblad semigroup `U(t)=exp(tG)`, the induced operator
norm is contractive: `||U(t)X||_∞ <= ||X||_∞`. Variation of constants therefore
gives the directly checkable bound

```text
||U(T)O - P(T)||_∞ <= ||P(0)-O||_∞ + ∫_0^T ||R(s)||_∞ ds.
```

The existing exact Pauli residual records can certify the right side by
integrating a positive pointwise norm majorant (monomial or Bernstein form).
Collocation equations are only a way to choose `P`; they are not part of the
certificate. This also explains the domain limitation: the argument requires
the stated semigroup contraction (or an explicitly certified growth factor).
It does not justify an error bound from a scalar spectral interval of a
nonnormal Lindblad generator.

A practical derivative-collocation construction chooses rational nodes `t_j`
and solves the linear equations

```text
P'(t_j) - G P(t_j) = 0   (or a prescribed defect target),
```

in a rational polynomial basis, with `P(0)` fixed or included as an initial
condition. The candidate is then converted to the checker’s power/Bernstein
representation and fully replayed. If the collocation system is solved only
approximately, its algebraic solve error is simply part of the exact residual
of the exported polynomial; no trust in floating-point node values is needed.
For a discontinuous Galerkin variant, each element has its own polynomial and
the bound must additionally charge the norm of every inter-element jump (or a
continuous reconstruction). Makridakis--Nochetto explicitly use a higher-order
reconstruction to restore continuity and derive an a posteriori residual
equation for dG(q)/Radau methods. [Makridakis--Nochetto, 2006](https://doi.org/10.1007/S00211-006-0013-6)

## Fair conventional baselines

* **Global collocation / cG.** A degree-`m` polynomial has `m+1` coefficients;
  derivative collocation requires `O(m)` generator applications per node in a
  dense coefficient formulation, plus an `O(m^3)` exact/dense solve. The fair
  comparison charges basis conversion, failed orders, residual replay, and
  coefficient growth. Classical collocation error analysis relates solution
  error to interpolation/residual error for smooth ODE solutions, but the
  v9 certificate must use the exact residual bound above rather than an
  interpolation theorem. (A primary numerical-analysis treatment is Ascher,
  Mattheij and Russell, *Numerical Solution of Boundary Value Problems for
  Ordinary Differential Equations*, SIAM, whose collocation chapters are
  indexed at [SIAM](https://epubs.siam.org/doi/book/10.1137/1.9781611975703).)
* **cG/dG time discretization.** Piecewise polynomial Galerkin methods are
  established residual-based baselines. dG requires jump/reconstruction terms;
  cG avoids jumps but couples the global solve. The Makridakis--Nochetto
  theorem applies to dissipative evolution operators under their angle and
  regularity hypotheses, which are stronger and different from finite Pauli
  operator-space assumptions. It is a baseline reference, not an automatic
  proof for every Lindblad representation.
* **Arnoldi/Krylov exponential action.** Hochbruck--Lubich prove convergence
  bounds for Arnoldi/Lanczos approximations to `exp(τA)v`, with dependence on
  numerical range/spectrum/pseudospectrum and matrix-vector products, and use
  them for time integration. [Hochbruck--Lubich, 1997](https://doi.org/10.1137/S0036142995280572)
  Defect-based work gives an integral representation and computable a posteriori
  bounds for Krylov propagators, while warning that quadrature defect estimates
  are not generally proven upper bounds. [Jawecki, 2020](https://arxiv.org/abs/2001.11922)
  Thus v9 should compare against the actual Arnoldi implementation and charge
  all orthogonalization, projected solves, and defect certification.
* **Adaptive quantum truncation.** Ensemble-rank truncation for Lindblad
  dynamics constructs an approximate Kraus map and truncates principal
  components; the paper reports `O(K N_H^2)` per step versus full
  matrix-matrix `O(N_H^3)` and suggests choosing rank by discarded weight.
  [McCaul, Jacobs and Bondar, 2020](https://arxiv.org/abs/2010.05399)
  This is a relevant computational baseline, but its displayed method has an
  `O(dt^2)` local approximation and discarded-weight heuristic; it is not by
  itself a rigorous operator-norm final-time certificate. It should be included
  only when its truncation error and all conversion costs are independently
  bounded under the same residual/norm contract.

## Expected cost versus Taylor

Taylor needs a recurrence `Y_{k+1}=G Y_k/(k+1)` and therefore roughly `m`
sparse generator applications, with exact support growth and the checker’s
residual work. Derivative collocation can avoid constructing `m+1` exponential
node values, but it introduces a coupled `O(m)`-unknown solve, repeated
generator applications at collocation equations, dense/fraction-free linear
algebra, and potentially severe coefficient growth after conversion to powers.
It is plausible only when the chosen polynomial basis keeps generator supports
small and the solve is materially cheaper than the Taylor recurrence at the
same certified bound. No such headroom follows from Chebyshev approximation
alone, especially for nonnormal dynamics. The proposed preflight should report
accepted bound, generator applications, solve/orthogonalization work,
coefficient bit lengths, residual replay, and refusals; it must not claim a
physical error from scalar spectral data.

## Theorem domain and decision boundary

The residual certificate is valid for a bounded finite-dimensional generator
whose Heisenberg semigroup is positive, unital, and operator-norm contractive,
with an exact or certified positive residual majorant. For a general
nonnormal matrix, a spectral radius/interval is insufficient; use a certified
semigroup bound or the residual argument. dG additionally needs jump terms or
continuous reconstruction. Krylov a posteriori formulas can depend on
dissipativity/Ritz geometry, and defect quadrature may be an estimate rather
than an upper bound. Quantum rank truncation requires a separately proved
discarded-component bound. These limitations keep the comparison fair while
leaving the conventional methods as strong baselines before any scientific
method acquisition is considered.


## Root source/cost correction

The arXiv 2001.11922 paper is single-authored by Tobias Jawecki; root verified that record directly. The related three-author exponential error-bound paper has a different identifier. A collocation solve is not merely O(m^3) independent of state dimension: for projected dimension k the direct stage system has dimension mk and a dense factorization cost scaling as O((mk)^3), with O((mk)^2) storage. V9 records these sizes. The Makridakis–Nochetto DOI did not resolve through root's browser in this pass, so its detailed theorem application is not independently established here. None of these external bounds replaces the exact local checker.
