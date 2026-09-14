# Exact Pauli closure scaling boundary

The experiment uses independent \(Z_i\) and nearest-neighbor \(X_iX_{i+1}\) generators, seeds the closure with all quadratic generators, and applies exact Pauli commutators. The quadratic/free-fermion prediction \(n(2n-1)\) is confirmed for \(n=2\) through 6: sizes 6, 15, 28, 45, 66. This verifies known Lie-space structure; it is not a new theorem. The Jordan–Wigner mapping and free-fermion solvability are standard; see Lieb, Schultz, and Mattis, *Annals of Physics* 16 (1961), 407–466, DOI [10.1016/0003-4916(61)90115-4](https://doi.org/10.1016/0003-4916(61)90115-4).

Adding \(Z_0Z_1\) as an additional control while retaining the original quadratic seeds leaves \(n=2\) unchanged (the added term is central in this small generated algebra), then grows the closure to 30 terms at \(n=3\), 126 at \(n=4\), and 510 at \(n=5\). At \(n=6\), cap 512 is reached and the result is explicitly incomplete; omitted witnesses are recorded in `results/scaling_results.json`. No dimension is inferred from that capped run. This is empirical closure-boundary evidence, not a universality proof.

The runner records `complete`, zero-residual verification, operation count, wall time, cap, work cap, and omissions. It uses exact Pauli algebra; capped output is only a lower bound.
