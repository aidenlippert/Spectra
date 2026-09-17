from fractions import Fraction as F
from decimal import Decimal, localcontext
import pytest

from experiments.marginal_hunt_car import add, mono
from experiments.marginal_symbolic import number_shift, product
from research.interacting_scaling_20260915.dictionary import ideal_basis
from research.rsi_discovery_20260916.v2.algebra import oracle
from research.rsi_discovery_20260916.v2.certified import CertifiedTemplates
from research.rsi_discovery_20260916.v2.frontier import coordinates, drift, exp_upper, advance, upper_round, construct, validate_observables
from research.rsi_discovery_20260916.v2.ideal import IdealSpan, verify_witness
from research.rsi_discovery_20260916.v2.quotient import reduce_dictionary


def test_new_checked_operation_reuses_only_verified_templates():
    kernel = CertifiedTemplates(lambda p: oracle(p["left"], p["right"]))
    for offset in (0, 10, 1000000):
        l, r = ((0, offset),), ((1, offset),)
        assert kernel(l, r) == oracle(l, r)
    assert kernel.checks == 1 and kernel.hits == 2
    assert not kernel.receipt()["full_molecular_sector_enumerated"]


def test_corrupt_new_checker_primitive_cannot_self_authorize():
    kernel = CertifiedTemplates(lambda p: {(): -1})
    with pytest.raises(ValueError, match="exact witness"):
        kernel(((0, 3),), ((1, 3),))
    assert not kernel.templates


def test_dictionary_reduction_proves_every_cross_term():
    modes, particles = 6, 2
    columns = [(str(i), product(number_shift(modes, particles), p))
               for i, p in enumerate(ideal_basis(modes, [], max_body=2))]
    a, b = mono(((0, 0),)), mono(((0, 2),))
    operators = [a, b, add(a, product(a, number_shift(modes, particles)))]
    accepted = reduce_dictionary(operators, modes, particles, IdealSpan(columns))
    assert accepted["declared_raw_ideal_map_equivalence_proved"]
    assert accepted["reduced_size"] == 2 and len(accepted["witnesses"]) == 6
    rejected = reduce_dictionary(operators, modes, particles, IdealSpan([("N-n", number_shift(modes, particles))]))
    assert rejected["physical_singlet_reduction_proved"]
    assert not rejected["declared_raw_ideal_map_equivalence_proved"]


def test_observable_control_algebra_closes_with_correct_sign():
    d = ({((1, 0), (0, 0)): F(1), ((1, 2), (0, 2)): F(-1)}, {})
    w = ({((1, 0), (0, 2)): F(1), ((1, 2), (0, 0)): F(1)}, {})
    y = ({}, {((1, 0), (0, 2)): F(-1), ((1, 2), (0, 0)): F(1)})
    basis = [(mono(()), {}), d, w, y]
    c, remainder = coordinates(drift(d, w), basis)
    assert c == [0, 0, 0, -2] and not remainder
    c, remainder = coordinates(drift(w, d), basis)
    assert c == [0, 0, 0, 2] and not remainder


def test_dynamical_integrator_and_outward_envelope():
    value, error, growth = advance([[F(1)]], [F(1)], F(1, 2), terms=20)
    with localcontext() as context:
        context.prec = 60
        exact = Decimal('0.5').exp()
        lo = Decimal(value[0].numerator) / Decimal(value[0].denominator)
        hi = Decimal((value[0] + error).numerator) / Decimal((value[0] + error).denominator)
        assert lo <= exact <= hi
    assert upper_round(F(1, 3)) >= F(1, 3)
    assert upper_round(F(-1, 3)) >= F(-1, 3)
    with pytest.raises(ValueError):
        exp_upper(F(1000))


def test_frozen_control_task_refuses_changed_hamiltonian_or_state():
    with pytest.raises(ValueError, match="original H8"):
        construct({"modes": 16, "particles": 8, "hamiltonian": []}, {})


def test_independent_witness_entry_point_rejects_hidden_energy_span():
    number = mono(((1, 0), (0, 0)))
    columns = [("n", number), ("n_minus_one", add(number, {(): F(-1)}))]
    assert not verify_witness(mono(()), columns, {0: F(1), 1: F(-1)})


def test_observable_proposals_cannot_change_domain_or_use_nonhermitian_data():
    w = ({((1, 0), (0, 2)): F(1), ((1, 2), (0, 0)): F(1)}, {})
    assert validate_observables([w], 4) == [w]
    with pytest.raises(ValueError, match="Hermitian"):
        validate_observables([({((1, 0), (0, 2)): F(1)}, {})], 4)
    with pytest.raises(ValueError, match="original Hamiltonian"):
        validate_observables([w], 2)
    with pytest.raises(ValueError, match="exact rational"):
        validate_observables([({((1, 0), (0, 0)): 0.5}, {})], 4)
    with pytest.raises(ValueError, match="particle-neutral"):
        validate_observables([({((1, 0),): 1}, {})], 4)
