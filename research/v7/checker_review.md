# V7 checker follow-up audit

The updated Bernstein conversion is algebraically correct. For power
coefficients `r_ell` on `[0,h]`, the implementation forms
`b_k = sum_{ell<=k} r_ell h^ell C(k,ell)/C(m,ell)` and charges `h/(m+1)` per
Bernstein coefficient. This is the standard exact conversion and integral
triangle envelope. Trimming trailing zero coefficients before choosing `m` is
also sound. Jump records remain full operator-norm witnesses with weight one.

The first draft of `tests/test_v7_dense_certificate.py` was invalid: it used
`U.conj().T` as the inverse of an eigensystem for a nonnormal damped generator,
and compared propagation from the target against a certificate whose initial
observable was empty. That draft is explicitly discarded. The repaired test
uses an independent scaling/squaring Taylor matrix exponential, constructs the
certificate with `{target: 1}`, and compares the dense endpoint against the
degree-12 Taylor candidate itself. It also cross-checks `gamma=0` against direct
unitary conjugation. Cases cover 2 and 3 qubits, `T=0.1,1`, and `gamma=0,0.2,2`;
all 12 repaired cases pass. This is only a numerical sanity check, not a proof.

Remaining shortcomings:

* The exact rational parser now rejects floats and bounds textual numerator and
  denominator lengths before constructing `Fraction`; the earlier parsing DoS
  issue is fixed. `expected_time` is authenticated against the sum of piece
  durations, fixing the implicit-horizon issue.
* Cost accounting is still explicitly structural/partial. It does not count
  every Fraction arithmetic, dictionary operation, Pauli multiplication, or
  bit operation. This is acceptable only if reports call these counters
  structural proxies and separately include the complete wall-clock/resource
  limit, as stated by the new design.
* The physical claim remains conditional: the checker verifies the declared
  Pauli Hamiltonian plus uniform local depolarization generator, but cannot
  establish that an external device has those couplings or noise.

Test log (from the populated workspace, with `PYTHONPATH=.`):

```
PYTHONPATH=. python tests/test_v7_dense_certificate.py
.
----------------------------------------------------------------------
Ran 1 test in 0.140s

OK
```


Root replay also passed `python -m unittest tests.test_v7_dense_certificate -v`
from the Documents/Spectra workspace; Python supports this namespace-package
import here. The earlier assertion that this invocation was unavailable was
incorrect.
