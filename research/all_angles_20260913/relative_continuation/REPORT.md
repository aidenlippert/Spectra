# Routes 37–44 and 53–56: relative continuation experiment

## Result

The relative perturbation certificate is valid and gives a real finite test
improvement. With `Eref=0`, `P=diag(0,g)`, and rational `R` chosen so that
`R + (1/5)P + (1/10)I = vv^T`, exact rational arithmetic verifies the premise.
For `g=1/10`, `v=(1/10,1/10)`, the dense independent eigenvalue is
`E0=-0.0912310562562` is diagnostic only. The relative lower bound is `-0.1`; the
exact rational trial `e_0` gives upper `-0.09`, hence certified width `0.01`.
The full residual-norm lower bound is
`-0.114142135624` (width `0.02291107937`): 2.61x narrower.

Across `g=0.1,0.2,0.4,0.8,1.6`, the diagnostic eigenvalue gap stays near
`0.009`, while the norm width grows from `0.023` to `0.320`. Thus the benefit survives a
coupling sweep for this tensor-product-like block family. These are
abstract energy units, not Hartree.

The adverse gap-closing sweep (`g=10^-1,10^-2,10^-4,10^-8`, with
`v=(1/20,1/20)`) remains valid even as the parent gap tends to zero. The
certificate becomes less informative only in the expected sense that its
fixed epsilon dominates; it does not assume a spectral gap. This specifically
guards against silently relying on a gapped parent.

An interacting cyclic chain attempt (`L=4,6,8`) uses a domain-wall parent and
a rank-one perturbation coupling every computational basis state. It has an
exact rational lower and rational `e_0` trial upper, with certified width
`0.0025` at each size. Dense ground values remain diagnostics only.

`results/all_angles_20260913/relative_continuation/receipt.json` contains the
machine-readable receipt. `test_relative.py` is both the exact rational
premise checker and the independent NumPy dense validation.

## Route accounting

* 37 frustration-free/parent controls: **known control only**. The prior
  projected-product parent ladder is exact width zero but rejects molecular
  fixtures; this experiment is a perturbed parent, with nonzero width and a
  verified improvement over `||R||`.
* 38 stoquastic gauges: **not tested here**; requires sign-structure search.
* 39 integrable/Bethe controls: **not tested here**; no claim for molecules.
* 40 supersymmetric/stabilizer/matchgate controls: **prior supersymmetric
  chain route is a known extensive/control result; no repeat**.
* 41 sign-removing gauges: **not tested here**.
* 42 fixed-node restrictions: **not tested here**.
* 43 similarity/transcorrelated transforms: **not tested here**.
* 44 continuous unitary flows: **not tested here**.
* 53 continuation along geometry: **proposal only**. The same inequality is
  stable under parameter changes if the rational PSD slack remains certified;
  no geometry fixture was run.
* 54 reaction-energy differences: **not tested**; total-energy certificates do
  not automatically cancel errors in differences.
* 55 response/forces: **not tested**.
* 56 local-observable guarantees: **not tested**.

## Proof and limits

For any normalized `x`, `xᵀRx >= -eta xᵀPx - eps ||x||²`; hence
`xᵀ(H-Eref I)x >= (1-eta)xᵀPx-eps||x||² >= -eps||x||²`.
Therefore `H >= (Eref-eps)I`. The checker validates the PSD premise by exact
construction (`vvᵀ`) and the eigenvalue independently in floating point.

The two-dimensional comparison and cyclic chain do not establish chemistry relevance,
asymptotic discovery cost, or a universal relative bound. For a tensor sum of
blocks, epsilon and total-energy error add extensively, while energy density
can remain bounded; the receipt reports single-block total energies only.
