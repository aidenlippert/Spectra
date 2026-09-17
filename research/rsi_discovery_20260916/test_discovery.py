from copy import deepcopy
from fractions import Fraction as F
import json
from pathlib import Path
import sqlite3
import sys

import pytest

from research.rsi_discovery_20260916 import candidates
from research.rsi_discovery_20260916.campaign import attempt, check_sources, compare_arms
from research.rsi_discovery_20260916.dialect import compile_candidate
from research.rsi_discovery_20260916.evaluator import verify_certificate, verify_grams, reference_gram
from research.rsi_discovery_20260916.ledger import Ledger, canonical, digest
from research.rsi_discovery_20260916.lemmas import cross_term_counterexample, car_lemma
from research.rsi_discovery_20260916.policy import Policy, replay
from research.rsi_discovery_20260916.worker import integer_gram, exact_margin


def fixture():
    return {"modes": 2, "particles": 1, "hamiltonian": [
        {"word": [[1, 0], [0, 0]], "coefficient": "1"},
        {"word": [[1, 1], [0, 1]], "coefficient": "2"}]}


def certificate():
    return {"kind": "discovered_full_N_certificate_v1", "fixture_sha256": digest(fixture()),
            "blocks": [
                {"alpha_count": 0, "lower_Ha": "1999/1000", "factor": [[0]], "factor_denominator": 1, "upper_vector": [1]},
                {"alpha_count": 1, "lower_Ha": "999/1000", "factor": [[0]], "factor_denominator": 1, "upper_vector": [1]}]}


def test_exact_gate_and_refusals():
    cert = certificate()
    receipt = verify_certificate(fixture(), cert)
    assert receipt["width_mHa"] == 1
    mutations = [
        lambda c: c.update(fixture_sha256="wrong"),
        lambda c: c["blocks"].pop(),
        lambda c: c["blocks"][1].update(lower_Ha="1001/1000"),
        lambda c: c["blocks"][1].update(lower_Ha="0"),
        lambda c: c["blocks"][1].update(factor_denominator=0),
        lambda c: c["blocks"][1].update(factor_denominator=True),
        lambda c: c["blocks"][1].update(factor=[[False]]),
        lambda c: c["blocks"][1].update(upper_vector=[0]),
    ]
    for mutate in mutations:
        bad = deepcopy(cert)
        mutate(bad)
        with pytest.raises(ValueError):
            verify_certificate(fixture(), bad)


@pytest.mark.parametrize("source", [
    "import os\ndef propose(payload): return 1",
    "def propose(payload): return payload.__class__",
    "def propose(payload): return __import__('os')",
    "def propose(payload): return open('/tmp/x')",
    "def propose(payload): return (lambda: 1)()",
    "def propose(payload): return 2 ** 999999",
    "@print\ndef propose(payload): return 1",
])
def test_restricted_candidate_cannot_access_host(source):
    with pytest.raises((ValueError, NameError)):
        compile_candidate(source, {})({})


def test_exact_component_large_signed_and_shape_failures():
    factors = [[[0]], [[2**120], [-2**100, 3]], [[-2], [0, 4], [7, 0, -8]]]
    payloads = [{"factor": f} for f in factors]
    outputs = [integer_gram(f) for f in factors]
    assert verify_grams(payloads, outputs)["instances"] == 3
    assert all(o == reference_gram(f) for f, o in zip(factors, outputs))
    outputs[1][0][0] += 1
    with pytest.raises(ValueError, match="equivalence"):
        verify_grams(payloads, outputs)
    with pytest.raises(ValueError):
        reference_gram([[True]])
    with pytest.raises(ValueError):
        verify_grams(payloads, [])


def test_acquired_remainder_is_equivalent_to_existing_checker():
    from research.compact_response_20260913.closure import check_factor
    h = [[10, 2], [2, 8]]
    factor = [[2], [1, 1]]
    for use_fast in (False, True):
        assert exact_margin(h, 1, F(0), factor, 1, use_fast) == check_factor(h, 1, F(0), factor, 1)


def test_ledger_crash_visibility_integrity_and_immutable_events(tmp_path):
    ledger = Ledger(tmp_path)
    key = ledger.blob("exact source")
    assert ledger.read_blob(key) == b"exact source"
    operation = ledger.append("operation_started", label="interrupted")
    ledger.append("attempt_started", id="interrupted_candidate")
    with ledger.measure("failed"):
        pass
    assert ledger.audit()["unfinished_operations"] == [operation]
    assert ledger.audit()["unfinished_attempts"] == ["interrupted_candidate"]
    with pytest.raises(sqlite3.IntegrityError):
        ledger.db.execute("DELETE FROM events")
    (tmp_path / "objects" / key).write_text("corruption")
    with pytest.raises(ValueError, match="Corrupted"):
        ledger.audit()


def history():
    return [
        {"id": "0", "parent": None, "family": "root", "depth": 0, "utility": 1., "cost_seconds": 2., "task_id": "t"},
        {"id": "1", "parent": "0", "family": "bad", "depth": 1, "utility": 0., "cost_seconds": 3., "task_id": "t"},
        {"id": "2", "parent": "1", "family": "good", "depth": 2, "utility": 5., "cost_seconds": 7., "task_id": "t"},
    ]


def test_replay_respects_ancestry_failures_stopping_and_no_future_visibility():
    class Spy(Policy):
        def choose(self, frontier, observed):
            assert all(set(n) == {"id", "parent", "family", "depth"} for n in frontier)
            assert all(n["id"] not in {v["id"] for v in frontier} for n in observed)
            return super().choose(frontier, observed)
    receipt = replay(history(), Spy(), 100)
    assert receipt["visited"] == ["0", "1", "2"]
    assert receipt["charged_seconds"] == 12
    assert replay(history(), Policy(patience=1), 100)["visited"] == ["0", "1"]
    assert replay(history(), Policy(), 3)["overshoot_seconds"] == 2
    with pytest.raises(ValueError, match="ancestors"):
        replay(history()[1:], Policy(), 100)
    cyclic = history()
    cyclic[0]["parent"] = "2"
    with pytest.raises(ValueError, match="Cyclic"):
        replay(cyclic, Policy(), 100)


def test_replay_refuses_unrecorded_actions():
    class Unrecorded(Policy):
        def choose(self, frontier, observed):
            return [{"id": "invented"}]
    with pytest.raises(ValueError, match="invented"):
        replay(history(), Unrecorded(), 100)


def test_counterexample_and_finite_car_proof():
    result = cross_term_counterexample(3)
    assert result["status"] == "refuted"
    x, y = map(F, result["witness"].values())
    assert 3*x*y > x*x+y*y
    assert cross_term_counterexample(2)["status"] == "conjecture"
    assert car_lemma()["status"] == "verified_instance"


def test_real_subprocess_certificate_path_and_failure_costs(tmp_path):
    ledger = Ledger(tmp_path)
    task = {"id": "tiny-fixed-task", "name": "tiny", "fixture": fixture()}
    valid = attempt(ledger, task, candidates.constructor(), "reference")
    assert valid["status"] == "verified_instance"
    assert valid["receipt"]["numerical_imports"] == []
    invalid = attempt(ledger, task, 'def propose(payload):\n    return {"blocks": []}\n', "invalid")
    assert invalid["status"] == "rejected"
    assert invalid["cost_seconds"] > 0
    assert ledger.audit()["unfinished_operations"] == []
    assert ledger.costs()["failed_operations"] == 1


def test_frozen_source_mutation_is_refused(tmp_path):
    import hashlib
    path = tmp_path / "gate.py"
    path.write_text("immutable")
    sources = {str(path): hashlib.sha256(path.read_bytes()).hexdigest()}
    check_sources(sources)
    path.write_text("changed acceptance")
    with pytest.raises(ValueError, match="Frozen"):
        check_sources(sources)


def test_incomplete_libraries_are_not_reused(tmp_path):
    from research.rsi_discovery_20260916.reuse import load_library
    ledger = Ledger(tmp_path)
    ledger.append("retained_method", category="constructor", name="development_selected", source_sha256=ledger.blob("unverified"))
    with pytest.raises(ValueError, match="completed"):
        load_library(tmp_path)


def test_malformed_program_output_is_retained_as_rejection(tmp_path):
    ledger = Ledger(tmp_path)
    task = {"id": "tiny-fixed-task", "name": "tiny", "fixture": fixture()}
    result = attempt(ledger, task, 'def propose(payload):\n    return 0\n', "malformed")
    assert result["status"] == "rejected"
    assert result["cost_seconds"] > 0
    assert len(ledger.events("attempt_finished")) == 1


@pytest.mark.parametrize("acquired", [False, True])
def test_adaptive_constructor_reaches_the_independent_gate(tmp_path, acquired):
    ledger = Ledger(tmp_path)
    task = {"id": "tiny", "name": "tiny", "fixture": fixture()}
    row = attempt(ledger, task, candidates.adaptive_constructor(acquired), "adaptive", acquired=acquired)
    assert row["status"] == "verified_instance"


def test_failed_arm_cannot_claim_a_cost_advantage():
    rows = [{"family": "baseline", "task_id": "t", "status": "verified_instance", "cost_seconds": 10},
            {"family": "candidate", "task_id": "t", "status": "inconclusive", "cost_seconds": 1}]
    assert compare_arms(rows, "baseline", "candidate", 0)["paid_cost_advantage_seconds"] is None
    rows[1]["status"] = "verified_instance"
    assert compare_arms(rows, "baseline", "candidate", 20)["paid_cost_advantage_seconds"] == -11


def test_reuse_charges_proof_reactivation_and_controller_overhead():
    from research.rsi_discovery_20260916.reuse import reuse_accounting
    rows = [{"family": "without_acquired_kernel", "task_id": "t", "status": "verified_instance", "cost_seconds": 10},
            {"family": "with_acquired_kernel", "task_id": "t", "status": "verified_instance", "cost_seconds": 5}]
    result = reuse_accounting(rows, prior_seconds=20, elapsed=24, generation_seconds=4)
    assert result["library_reactivation_and_overhead_seconds"] == 5
    assert result["paid_cost_advantage_seconds"] == -20
