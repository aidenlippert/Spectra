# Exact molecular operator construction

The final tensor_operator.py maps all 2,912 supplied CAR terms into products of exact local Jordan–Wigner matrices. Creation is encoded (1, mode); annihilation is (0, mode). In the occupation basis, creation is [[0,0],[1,0]] and annihilation is its transpose. No configuration basis is generated.

At each cut, a sparse coefficient matrix is indexed by (left state, local matrix) and the remaining suffix word. An exact minimum vertex cover of its nonzero pattern gives a factorization into the next operator bond. This is a structural factorization, not a claim of minimum algebraic rank. Coefficients remain rational. Re-expanding the final graph reproduces every original product coefficient exactly.

The accepted H8 construction has 3,916 symbolic local connections and bond widths [1,4,16,39,62,79,100,125,154,125,100,79,62,39,16,4,1]. Construction plus symbolic replay took 0.4499 seconds internally; the bounded process took 0.5373 seconds.

Independent small-system occupied-bit tests check operator action, including nilpotence, non-normal-ordered input, lower-degree terms and quartic terms. sparse_mpo_edges exposes (left,right,bra,ket,Fraction) entries.

An early implementation interpreted the input encoding incorrectly and reported eight Pauli words. That result is invalid and has been replaced. The retained old diagnostic is explicitly excluded in artifact_disposition.json.

This exact operator front end does not establish economical state evolution. The bond-64 action/compression bound remains too loose to certify the frozen weak-control task.
