# Symmetric matched-flavor variational upper bound

For `M=2m`, fixed `N=m`, use the normalized ansatz supported on configurations
with exactly one fermion in each matched pair `(i,i+m)`. Let `k` be the number
of left occupations and assign a common amplitude `a_k` to each of the
`binom(m,k)` flavor subsets. In a flavor-ordered creation convention, a gauge
choice makes every matched hop between the `k` and `k+1` sectors have the same
sign. Counting edges gives

`||psi||^2 = sum_k binom(m,k) a_k^2`,

`<V> = sum_k binom(m,k)[binom(k,2)+binom(m-k,2)]a_k^2`,

and

`<T> = -2t sum_{k=0}^{m-1} binom(m,k)(m-k)a_k a_{k+1}`.

The second formula counts the within-half density pairs. The third counts the
`m-k` available right-to-left matched hops from each configuration and includes
the Hermitian reverse hop. Equivalently, in normalized Dicke coordinates
`c_k=sqrt(binomial(m,k))*a_k`, the variational matrix is tridiagonal with

`D_kk=binom(k,2)+binom(m-k,2)`,
`D_{k,k+1}=-t*sqrt((k+1)(m-k))`.

The lowest eigenvector of this `(m+1)`-by-`(m+1)` matrix supplies a valid
Rayleigh-Ritz upper bound. It is an ansatz bound for every `m`; it does not
assert the global ground state, especially after asymmetric hopping breaks the
permutation symmetry.

Independent verification needs no Fock enumeration or NumPy: for integer
amplitudes, compute the two displayed integer sums using exact `binom`, then
compare `energy/norm` with the reported fraction. For a small sign audit,
construct each flavor-ordered determinant by applying creation operators and
evaluate one matched hop using the CAR parity `(-1)^(number of occupied modes
below the acted mode)`; after the stated gauge convention, all contributing
edges must agree with the negative hopping sign.
