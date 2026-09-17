from fractions import Fraction as F
import pytest

from experiments.marginal_symbolic import encode
from research.rsi_discovery_20260916.v2.frontier_search import checked_basis


def test_candidate_artifact_round_trip_keeps_gaussian_rational_observables_exact():
    real = {((1, 0), (0, 0)): F(2, 7)}
    imag = {((1, 0), (0, 2)): F(-3, 11), ((1, 2), (0, 0)): F(3, 11)}
    artifact = [{"real": encode(real), "imaginary": encode(imag)}]
    assert checked_basis(artifact, 4) == [(real, imag)]


def test_candidate_artifact_cannot_smuggle_nonhermitian_correlations_to_checker():
    artifact = [{"real": encode({((1, 0), (0, 2)): F(1)}), "imaginary": []}]
    with pytest.raises(ValueError, match="Hermitian"):
        checked_basis(artifact, 4)
