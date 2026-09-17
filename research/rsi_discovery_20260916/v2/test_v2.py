from dataclasses import replace
from fractions import Fraction as F

import pytest

from experiments.marginal_hunt_car import add, mono
from experiments.marginal_symbolic import number_shift, product
from research.rsi_discovery_20260916.ledger import Ledger
from research.rsi_discovery_20260916.v2.algebra import (BASE_SOURCE, oracle, gram_reference,
    task_words, ladder_identity, product_cases, validate_map)
from research.rsi_discovery_20260916.v2.ideal import IdealSpan, number_null_cross_terms, verify_witness
from research.rsi_discovery_20260916.v2.language import compile_program
from research.rsi_discovery_20260916.v2.records import Research, Method, State, Action, choose_action


def test_ladder_oracle_independent_of_normal_ordering():
    for l, r in product_cases(194, 120):
        assert ladder_identity(l, r, oracle(l, r))
    assert not ladder_identity(((0, 0),), ((1, 0),), {(): -1})


def test_generated_grammar_supports_real_algorithms():
    fn = compile_program('def helper(x):\n    return x + 1\ndef propose(payload):\n    x = 0\n    while x < payload:\n        x = helper(x)\n    return x\n', {})
    assert fn(8) == 8


@pytest.mark.parametrize("source", [
    "import os\ndef propose(payload):\n    return 1",
    "def propose(payload):\n    return payload.__class__",
    "def propose(payload):\n    return (lambda: 1)()",
    "def propose(payload):\n    return 2 ** payload",
])
def test_grammar_refuses_escape_constructs(source):
    with pytest.raises(ValueError):
        compile_program(source, {})


def test_map_has_both_gram_cross_terms():
    words = [((0, 0),), ((0, 1),)]
    result = gram_reference(words)
    assert result[((1, 0), (0, 1))][(0, 1)] == 1
    assert result[((1, 1), (0, 0))][(0, 1)] == 1
    validate_map(result, words)
    fn = compile_program(BASE_SOURCE, {"oracle": oracle})
    assert fn({"words": words}) == result


def test_support_graphs_are_direct_and_transfer_changes_inputs():
    a = task_words({"family": "collective", "spatial": 4})
    b = task_words({"family": "star", "spatial": 5, "mode_offset": 31, "mode_stride": 3, "adjoint": True})
    assert a and b and len(b) > len(a)
    assert min(i for w in b for _, i in w) >= 31


def test_energy_scalar_never_enters_vanishing_ideal():
    for label in ("energy", "free_scalar", "identity", "dishonest_alias"):
        with pytest.raises(ValueError):
            IdealSpan([(label, mono(()))])
    assert not verify_witness(mono(()), [("energy", mono(()))], {0: 1})
    assert not verify_witness(mono(()), [("dishonest_alias", mono(()))], {0: 1})
    n = number_shift(4, 2)
    with pytest.raises(ValueError, match="span the energy scalar"):
        IdealSpan([("constraint_one", n), ("constraint_two", add(n, mono(())))])


def test_exact_span_does_not_silently_add_columns():
    n = number_shift(4, 2)
    span = IdealSpan([("N-n", n)])
    good = span.witness({w: 7 * c for w, c in n.items()})
    assert good["equivalent_in_declared_ideal"]
    assert verify_witness({w: 7 * c for w, c in n.items()}, span.columns, good["coordinates"])
    assert not span.witness(mono(((1, 0), (0, 0))))["equivalent_in_declared_ideal"]


def test_physical_null_does_not_imply_truncated_sos_equivalence():
    n = number_shift(4, 2)
    span = IdealSpan([("N-n", n)])
    result = number_null_cross_terms(4, 2, mono(((0, 1),)), mono(((0, 0),)), span)
    assert result["physical_fixed_N_equality"]
    assert result["includes_both_cross_terms_and_null_square"]
    assert not result["equivalent_in_declared_ideal"]
    # Explicitly admitting this independently justified null constraint changes the relaxation.
    expanded = IdealSpan([("N-n", n), ("explicit_new_constraint", result["difference"])])
    accepted = expanded.witness(result["difference"])
    assert accepted["equivalent_in_declared_ideal"]
    assert verify_witness(result["difference"], expanded.columns, accepted["coordinates"])


def research(tmp_path):
    return Research(Ledger(tmp_path), {"H": "fixed", "N": 2, "tolerance": "1/625"}, {"fingerprint": "env1"})


def test_dependency_dag_includes_all_context_reads(tmp_path):
    r = research(tmp_path)
    a = r.observe(action=Action.OPEN, state=State.ACCEPTED, claim="lemma")
    b = r.observe(action=Action.COUNTEREXAMPLE, state=State.REFUTED, claim="separate conjecture")
    c = r.observe(action=Action.COMBINE, state=State.ACCEPTED, claim="method", dependencies=(a,), reads=(b,))
    with pytest.raises(ValueError, match="before acquisition"):
        r.replay([a, c], "env1")
    assert all(not x["needs_live_rerun"] for x in r.replay([a, b, c], "env1"))


def test_unvisited_and_changed_environment_are_unknown(tmp_path):
    r = research(tmp_path)
    a = r.observe(action=Action.OPEN, state=State.ACCEPTED, claim="one observed branch")
    assert r.replay([a], "different")[0]["state"] == State.UNKNOWN.value
    assert r.replay(["missing"], "env1")[0]["needs_live_rerun"]


def test_method_is_automatically_callable_only_after_acquisition(tmp_path):
    r = research(tmp_path)
    source = r.ledger.blob("def propose(payload):\n    return payload + 1\n")
    evidence = r.observe(action=Action.PROVE, state=State.ACCEPTED, claim="bounded test")
    method = Method("increment", source, {"payload": "integer"}, {"result": "integer"},
                    ("scoped evidence",), ("algebra",), (evidence,), {"seconds": 1}, (),
                    ("unproved universal contract",), "propose(integer)", "enabling")
    key = r.retain(method)
    assert key in r.applicable("algebra")
    assert r.activate(key, {})(3) == 4
    with pytest.raises(ValueError):
        r.activate(key, {}, available=[])
    assert not r.applicable("unrelated")
    assert r.ledger.audit()["chain_valid"]


def test_policy_uses_observed_improvement_and_keeps_unvisited_branches_unknown():
    archive = [{"before": State.ACCEPTED.value, "action": Action.SEARCH.value,
                "measured": True, "accepted": True, "improves_parent": False, "cost_seconds": 5},
               {"before": State.ACCEPTED.value, "action": Action.COMBINE.value,
                "measured": True, "accepted": True, "improves_parent": True, "cost_seconds": 5}]
    assert choose_action(State.ACCEPTED.value, archive) == Action.COMBINE
    assert choose_action(State.ACCEPTED.value, archive, fixed=True) == Action.SEARCH
    assert choose_action(State.REPAIRABLE.value, archive) == Action.DEBUG
