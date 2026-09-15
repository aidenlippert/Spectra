from fractions import Fraction
import pytest

from experiments.v10_exp_readout import evaluate_modes


def test_exact_constant_and_zero_time():
    terms = [(Fraction(0), Fraction(0), 0, Fraction(3), Fraction(-2))]
    mid, err, cost = evaluate_modes(terms, Fraction(0), Fraction(1, 1000))
    assert mid == (Fraction(3), Fraction(-2))
    assert err == 0


def test_real_decay_is_enclosed():
    mid, err, _ = evaluate_modes([(Fraction(-1), 0, 0, 1, 0)], Fraction(1), Fraction(1, 100))
    # exp(-1) is between 0.36 and 0.37; the rational disk must contain it.
    assert mid[0] - err <= Fraction(37, 100)
    assert mid[0] + err >= Fraction(36, 100)
    assert err <= Fraction(1, 100)


def test_equal_exponents_cache_and_budget():
    rows = [(0, 1, 0, 1, 0), (0, 1, 1, 1, 0)]
    mid, err, cost = evaluate_modes(rows, Fraction(1, 2), Fraction(1, 100))
    assert cost["unique_exponents"] == 1
    assert cost["cache_hits"] == 0  # cache is per unique exponent
    assert err <= Fraction(1, 100)


def test_refusals():
    with pytest.raises(ValueError):
        evaluate_modes([(1, 0, 0, 1, 0)], 1, Fraction(1, 10))
    with pytest.raises(ValueError):
        evaluate_modes([(0, 0, 33, 1, 0)], 1, Fraction(1, 10))
