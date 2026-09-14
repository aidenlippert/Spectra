"""Exact toy check of the phase-estimation ground-overlap sampling barrier.

For |phi> = sqrt(p)|0> + sqrt(1-p)|1>, exact phase estimation returns E0
with probability p. Repeating until seeing E0 has geometric mean 1/p.
Amplitude amplification changes this to O(1/sqrt(p)) oracle uses, but needs
coherent reflections/state preparation; it does not make p disappear.
"""
from fractions import Fraction


def expected_trials(p: Fraction) -> Fraction:
    assert 0 < p <= 1
    return 1 / p


def success_after(n: int, p: Fraction) -> Fraction:
    """Probability of at least one exact ground outcome in n independent shots."""
    return 1 - (1 - p) ** n


def main() -> None:
    for p in (Fraction(1, 4), Fraction(1, 16), Fraction(1, 100)):
        n = 0
        while success_after(n, p) < Fraction(9, 10):
            n += 1
        print(f"p={p}: E[shots]={expected_trials(p)}, shots_90={n}")
    # A falsifiable identity: geometric expectation is exactly 1/p.
    assert expected_trials(Fraction(1, 16)) == 16
    assert success_after(16, Fraction(1, 16)) < Fraction(9, 10)


if __name__ == "__main__":
    main()
