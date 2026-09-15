"""Exact small checks for the local-overlap compression theorem.

Stdlib only.  The example uses Pauli matrices as a finite CAR/Jordan-Wigner
toy: it intentionally includes two noncommuting overlapping supports.  The
global 8-by-8 matrix is used only to check the toy identity; the advertised
certificate bound uses local norms and symbolic supports.
"""

from fractions import Fraction
from itertools import product


def mm(a, b):
    n, m, p = len(a), len(b[0]), len(b)
    return [[sum(a[i][k] * b[k][j] for k in range(p)) for j in range(m)] for i in range(n)]


def add(a, b):
    return [[a[i][j] + b[i][j] for j in range(len(a[0]))] for i in range(len(a))]


def kron(a, b):
    return [[a[i // len(b)][j // len(b[0])] * b[i % len(b)][j % len(b[0])] for j in range(len(a[0]) * len(b[0]))]
            for i in range(len(a) * len(b))]


def eye(n):
    return [[Fraction(int(i == j)) for j in range(n)] for i in range(n)]


I = [[Fraction(1), Fraction(0)], [Fraction(0), Fraction(1)]]
X = [[Fraction(0), Fraction(1)], [Fraction(1), Fraction(0)]]
Z = [[Fraction(1), Fraction(0)], [Fraction(0), Fraction(-1)]]


def op3(a, b, c):
    return kron(kron(a, b), c)


def assert_eq(a, b, label):
    assert a == b, label


def main():
    # Fixed-N warning: arbitrary multipliers of (n-N) do not vanish after
    # projection.  For one mode and N=0, P0 a n a^dagger P0 = P0 exactly.
    a = [[Fraction(0), Fraction(1)], [Fraction(0), Fraction(0)]]
    adag = [[Fraction(0), Fraction(0)], [Fraction(1), Fraction(0)]]
    number = mm(adag, a)
    bad = mm(mm(a, number), adag)
    assert bad == [[Fraction(1), Fraction(0)], [Fraction(0), Fraction(0)]]
    assert bad[0][0] != 0, "arbitrary multiplier must not be treated as sector-zero"

    # R01 = X_0 Z_1 and R12 = X_1 X_2 overlap at mode/cluster 1 and anticommute.
    r01 = op3(X, Z, I)
    r12 = op3(I, X, X)
    full = add(r01, r12)
    zero = [[Fraction(0) for _ in range(8)] for _ in range(8)]
    # Exact local norms are one.  The triangle certificate therefore prices eta=2.
    local_norms = [Fraction(1), Fraction(1)]
    eta = sum(local_norms)
    assert eta == 2

    # Anticommutation and R^2=2I are exact rational matrix identities.
    assert_eq(add(mm(r01, r12), mm(r12, r01)), zero, "overlap terms must anticommute")
    assert_eq(mm(full, full), [[2 * x for x in row] for row in eye(8)], "R^2=2I")
    # Hence ||R||=sqrt(2), while the executable local certificate reports eta=2.
    # Verify the squared norm inequality without floating point: 2 <= eta^2=4.
    assert Fraction(2) <= eta * eta

    # Fixed-N check: retain only computational-basis states with Hamming weight N.
    # Compression is safe: P_N R P_N is a principal submatrix, so its Frobenius
    # norm is at most the unrestricted operator bound eta.
    for n_particles in range(4):
        states = [bits for bits in product((0, 1), repeat=3) if sum(bits) == n_particles]
        indices = [bits[0] * 4 + bits[1] * 2 + bits[2] for bits in states]
        projected = [[full[i][j] for j in indices] for i in indices]
        frob_sq = sum(x * x for row in projected for x in row)
        assert frob_sq <= eta * eta, (n_particles, frob_sq)

    # Size/verification accounting is symbolic and never constructs a 2^m matrix.
    # Chain supports {i,i+1}: p clusters, |E|=p-1, degree <=2, local block size 4.
    for p in (2, 4, 8, 32, 128):
        edges = p - 1
        assert edges <= 2 * p
        described_coefficients = 2 * edges  # two rational Pauli coefficients/support
        global_dimension = 2 ** p
        assert described_coefficients < global_dimension
        print(f"p={p} edges={edges} local_coefficients={described_coefficients} "
              f"global_dim={global_dimension} eta_bound={edges}")

    print("PASS: noncommuting overlap, fixed-N compression, exact residual bound, symbolic size sweep")


if __name__ == "__main__":
    main()
