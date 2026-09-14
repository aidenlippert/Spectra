# Final correction status

The original observations below describe intermediate agent drafts and are
superseded where noted. Root replaced the locality script with actual CAR
word_product checks and matching primal/dual certificates; use
locality_dual_exact.json. The periodic bond checker now has a stdlib dagger,
computes the commutator of energy terms, and evaluates its determinant upper
through DeterminantOracle. Cluster code checks exact four-state local
matrices; its global-chain bound remains an analytic composition theorem,
not a fully emitted global SOS artifact. The new independent audits for
stoquastic_chain_certificate and fermionic_ratio_chain cover their actual
implementations and are separate from these early diagnostics.

---

# Structural verification audit

The structural claims have useful finite evidence, with several scope and
implementation corrections required before presenting them as general
theorems.

## `low_rank_structural_check.py`

The `build()` result is a valid exact zero certificate when run with ordinary
Python: `verify()` returns lower `0`, residual `0`, and four bond factors. The
all-up determinant is indeed a zero witness for the displayed annihilation
factors, and adjacent factors overlap in mode support. However:

* The file imports `dagger` from `experiments.marginal_coefficient`, which
  imports NumPy. Thus `python -S low_rank_structural_check.py` fails before
  verification with `ModuleNotFoundError: numpy`; the exact production
  verifier itself is stdlib-only, but this structural generator is not.
* `bond_words()` uses `(site + 1) % sites`, making a periodic ring, while the
  module description calls it a chain. The factor count and overlap claim are
  therefore ring results. If an open ladder is intended, remove the modulo;
  otherwise label it periodic.
* `verify()` checks the exact polynomial identity and norm-≤1 CAR monomial
  residual, but it does not independently evaluate the all-up witness. The
  witness claim is transparent and can be checked separately, but should not
  be described as an output of `verify()`.

## `cluster_structural_check.py`

The script's dynamic program and assertions correctly establish the stated
bookkeeping values for the toy local energies and edge-norm penalties. It does
not construct or replay the fermionic dimer Hamiltonians, prove the local
energies `b(0)=1,b(1)=0,b(2)=1`, or verify the operator norm equality for the
edge. The comment correctly states that local dimer occupancy is not
conserved. Consequently this is an exact accounting check conditional on the
declared local model, not an exact CAR certificate generator. Fixed coupling
has an extensive width; the scaled-coupling table is not fixed-accuracy
scaling in system size because the physical coupling changes as `1/L`.

## `locality_dual_counterexample.py`

The reduced one-particle argument and the displayed `min(0,1-|lambda|)` value
are mathematically clear, and the negative determinant is a valid finite
warning against replacing a full PSD pair block by individual monomial
squares. The script itself does not perform exact matrix algebra: `matrices()`
uses NumPy floats and is unused by `run()`, while `pseudofunctional_obstruction()`
returns hard-coded symbolic facts. Call this an analytic counterexample with
an executable data sweep, not an independently exact CAR replay. Also, the
claim concerns a deliberately restricted cone and cannot be generalized to
all local SOS or chordal Gram constructions.

## Complexity and scope

`verify()` is polynomial in the emitted certificate's number of factor rows,
word lengths, and rational bit lengths when the maximum factor word degree is
fixed (the current admitted factors have degree at most three). It is not a
polynomial-time guarantee for unrestricted degree, since CAR expansion can
grow combinatorially, nor in physical system size when certificate output or
coefficient bit length grows exponentially. Chordal/local block storage is
linear only under fixed block width, bounded overlap, and a supplied local
residual decomposition. The scripts do not establish fixed total energy error
as `N` grows: the cluster penalty is extensive at fixed coupling and fixed
cluster size. The four-ladder controls therefore support conditional finite
families, not a universal chemistry or representability theorem.
