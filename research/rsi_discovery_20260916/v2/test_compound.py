from research.rsi_discovery_20260916.v2.compound import reference, evaluate, BASE_SOURCE, operators, capabilities_for
from research.rsi_discovery_20260916.v2.compound_campaign import operational_gate


def test_polynomial_map_accounts_for_both_cross_terms_and_zero_operators():
    a, b, c = ((0, 0),), ((0, 1),), ((0, 2),)
    value = reference([{a: 1, b: 2}, {a: -3, c: 1}, {}])
    assert value[((1, 0), (0, 0))][(0, 0)] == 1
    assert value[((1, 1), (0, 1))][(0, 0)] == 4
    assert value[((1, 0), (0, 0))][(0, 1)] == -6
    assert value[((1, 0), (0, 2))][(0, 1)] == 1
    assert value[((1, 2), (0, 1))][(0, 1)] == 2
    assert all(2 not in coordinate for row in value.values() for coordinate in row)


def test_polynomial_map_keeps_integers_above_float_precision_exact():
    a = ((0, 0),)
    coefficient = 2**71 + 37
    value = reference([{a: coefficient}])
    assert value[((1, 0), (0, 0))][(0, 0)] == coefficient * coefficient


def test_baseline_runs_through_independent_polynomial_gate():
    spec = {"seed": 893, "base": {"family": "collective", "spatial": 3, "weight": -1},
            "unique_words": 9, "operators": 4, "terms": 3, "relations": True, "zero_operator": True}
    value = evaluate({"source": BASE_SOURCE, "tasks": [spec], "repeats": 2})
    assert value["status"] == "verified_encoded_claim"
    assert value["tasks"][0]["operator_count"] == 4


def test_wrong_polynomial_program_returns_a_scoped_counterexample():
    spec = {"seed": 131, "base": {"family": "collective", "spatial": 3, "weight": -1},
            "unique_words": 6, "operators": 3, "terms": 2}
    value = evaluate({"source": "def propose(payload):\n    return {}\n", "tasks": [spec], "repeats": 1})
    assert value["status"] == "counterexample_to_encoded_claim"
    assert "expected" in value and "row" in value


def test_frozen_library_includes_existing_exact_monomial_map():
    words = [((0, 0),), ((0, 2),), ((0, 0),)]
    capability = capabilities_for()["monomial_map"]
    actual = capability({"words": words, "cache": {}})
    assert actual == reference([{word: 1} for word in words])


def test_inherited_execution_speed_alone_cannot_pass_discovery_gate():
    favorable = {"comparable": True, "D_over_A_wins": 6, "paired_sign_p_one_sided": 1 / 64,
                 "median_speedups": {"A/D": 3, "A/B": 2, "C/D": 2}}
    common = {**favorable, "median_speedups": {"A/D": 1, "A/B": 1, "C/D": 1}}
    assert not operational_gate(favorable, common, stable=True, generated_and_consumed=True)
    assert not operational_gate(favorable, favorable, stable=True, generated_and_consumed=False)
    assert operational_gate(favorable, favorable, stable=True, generated_and_consumed=True)
