# Stabilizer reduction of the complete quartic block

For the ten-mode matched model, the neutral charge block of the complete
degree-four dictionary has dimension 186. Its stabilizer under the exact
matched-flavor action is the full `S5 x C2`, of order 240. I represented each
group element as a signed permutation matrix on the CAR words and formed the
integer group-algebra operator

\[
 R=\sum_{g\in G}(k_g+1)(U_g+U_g^T).
\]

The spectrum of `R` has eigenspace dimensions

```
11, 2, 6, 144, 2, 11, 10
```

and the dimensions sum to 186. A Gram matrix invariant under the stabilizer
commutes with `R`, so replacing the 186-dimensional PSD variable by PSD
variables on these eigenspaces is lossless at the exact symmetry level. The
largest multiplicity block is still 144; this particular random group-algebra
element therefore does not provide the desired compression by itself. Several
independent group-algebra elements, or a character/projector decomposition,
are required to split that multiplicity space.

The signed action was checked on every word in the block and every one of the
240 stabilizer elements; all images remained in the block. This computation is
purely coefficient-level and constructs no occupation-sector matrices.
