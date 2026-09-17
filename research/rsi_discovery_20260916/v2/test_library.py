import pytest
from dataclasses import replace

from research.rsi_discovery_20260916.ledger import Ledger
from research.rsi_discovery_20260916.v2.library import MethodLibrary
from research.rsi_discovery_20260916.v2.records import Research, Method, State, Action


def method(ledger, name, code, evidence, deps=()):
    return Method(name, ledger.blob(code), {"payload": "integer"}, {"result": "integer"},
        ("scoped test",), ("algebra",), (evidence,), {}, deps, ("new inputs need a gate",), "propose(payload)", "enabling")


def test_later_process_reopens_and_automatically_installs_acquired_operation(tmp_path):
    ledger = Ledger(tmp_path)
    r = Research(ledger, {"test": "fixed"}, {"fingerprint": "original"})
    evidence = r.observe(action=Action.PROVE, state=State.ACCEPTED, claim="increment on scoped inputs")
    first = r.retain(method(ledger, "increment", "def propose(payload):\n    return payload + 1\n", evidence))
    early = ledger.events()[-1]["event_hash"]
    second_evidence = r.observe(action=Action.COMBINE, state=State.ACCEPTED, claim="composition", reads=(evidence,))
    second = r.retain(method(ledger, "twice", "def propose(payload):\n    return increment(increment(payload))\n", second_evidence, (evidence,)))
    reopened = MethodLibrary(Ledger(tmp_path))
    assert reopened.activate(second, {})(4) == 6
    assert set(reopened.applicable("algebra")) == {first, second}
    assert set(MethodLibrary(Ledger(tmp_path), through=early).applicable("algebra")) == {first}


def test_later_method_cannot_access_unread_capability(tmp_path):
    ledger = Ledger(tmp_path)
    r = Research(ledger, {"test": "fixed"}, {"fingerprint": "original"})
    first_evidence = r.observe(action=Action.PROVE, state=State.ACCEPTED, claim="first")
    r.retain(method(ledger, "increment", "def propose(payload):\n    return payload + 1\n", first_evidence))
    independent = r.observe(action=Action.PROVE, state=State.ACCEPTED, claim="independent bounded test")
    second = r.retain(method(ledger, "bad_dependency", "def propose(payload):\n    return increment(payload)\n", independent))
    fn = MethodLibrary(Ledger(tmp_path)).activate(second, {})
    with pytest.raises(NameError):
        fn(1)


def test_unknown_horizon_does_not_expose_later_library(tmp_path):
    with pytest.raises(ValueError, match="horizon"):
        MethodLibrary(Ledger(tmp_path), through="unvisited")


def test_conditional_representation_survives_reopen_without_becoming_callable(tmp_path):
    ledger = Ledger(tmp_path)
    r = Research(ledger, {"test": "fixed physical target"}, {"fingerprint": "original"})
    evidence = r.observe(action=Action.PROVE, state=State.CONDITIONAL,
                         claim="Exact defect enclosure is too wide for the target")
    branch = replace(method(ledger, "unresolved", "def propose(payload):\n    return payload\n", evidence),
                     archive="high_risk")
    key = r.retain(branch)
    reopened = MethodLibrary(Ledger(tmp_path))
    assert key in reopened.methods
    assert key not in reopened.applicable("algebra")
    with pytest.raises(ValueError, match="unresolved branch"):
        reopened.activate(key, {})
