# Target 7: a concrete stability obligation

The [May 2026 preprint](https://arxiv.org/html/2605.10943) by Balasubramanian, Davydova and Lin constructs a 3D CSS memory with exponential lifetime under its specified low-temperature thermal dynamics. It explicitly defers TQO-2, required for its proposed arbitrary-local-perturbation stability argument. This establishes neither a physical material nor operation at 300 K. We inspected version 1 and its [version record](https://arxiv.org/abs/2605.10943).

## Diagnostic, derived from stabilizer algebra

The TQO-2 framework compares local and global ground-space reductions. For a small region A and allowed enlarged region B, the local ground space on B and global ground space must induce the same reduced support on A. See [Bravyi, Hastings and Michalakis](https://arxiv.org/abs/1001.0344) for the geometric and Hamiltonian assumptions of the stability theorem. Matching this condition on a few regions does not establish it uniformly.

For positive-sign CSS codes, global stabilizers supported inside A form a binary subspace. Local stabilizers generated using checks wholly supported inside B form another. The supplied `experiments/stabilizer.py` computes both by GF(2) elimination and returns a missing constraint when containment fails. X and Z sectors are handled separately after checking commutation; this restriction avoids silently discarding sign-consistency obligations for general stabilizer presentations.

The manufactured example uses X1X2 and X2X3. With A=B={1,3}, the global product X1X3 is unavailable from buffer-local generators. Including site 2 in B restores it. The disconnected A in this example is a diagnostic test input, not a region satisfying every hypothesis of a published TQO-2 theorem. It does not reconstruct the recursive memory code.

A bounded next attack is to implement the actual recursive generators, specify valid interior/buffer geometry, and search junctions for a witness. A uniform positive result needs an induction proving local generation at every recursive scale. Finite successes cannot replace that induction. Spectral perturbation stability and thermal memory lifetime remain separate proof obligations.

## Quantifying an ambient-temperature proposal

The following arithmetic is illustrative, independent of any code construction. At 300 K, kBT is approximately 25.85 meV. An Arrhenius single-event estimate with attempt rate 10^12/s and a 10^9-second horizon requires a barrier around 25.85 meV × log(10^21), or 1.25 eV. This is not a memory theorem: correlated errors, the number of escape paths, entropy, disorder and bath coupling also matter.

An engineering claim should therefore supply an interaction Hamiltonian, a derivation of its effective protective terms, allowed disorder, a microscopic bath model, initialization and decoding costs, and a lifetime inequality in physical units. Weak effective couplings introduced by an implementation gadget must be charged before comparing with kBT.

A successful local-generation result would close one mathematical obligation. It would not automatically supply initialization, physical realization, thermal robustness or protected gates.
