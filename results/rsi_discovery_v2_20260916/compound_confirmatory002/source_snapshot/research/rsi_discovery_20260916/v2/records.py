"""Typed scientific state, operational methods and causal dependency replay.

The history stores observations, not counterfactual outcomes. A changed environment
invalidates outcome reuse even if the source of a candidate is identical.
"""
from dataclasses import dataclass, asdict
from enum import Enum
import hashlib
import importlib.metadata
from pathlib import Path
import platform
import sys

from research.rsi_discovery_20260916.ledger import digest
from research.rsi_discovery_20260916.v2.language import compile_program


class State(str, Enum):
    OPEN = "open_mechanism"
    CONJECTURE = "unproved_conjecture"
    ACCEPTED = "verified_encoded_claim"
    REFUTED = "counterexample_to_encoded_claim"
    REPAIRABLE = "implementation_failure"
    NUMERICAL = "numerical_failure"
    REPRESENTATION = "representation_obstruction"
    COST = "resource_obstruction"
    CONDITIONAL = "conditional_proof_obligation"
    SUSPENDED = "suspended"
    UNKNOWN = "unobserved"


class Action(str, Enum):
    OPEN = "open_mechanism"
    SEARCH = "continue_search"
    DEBUG = "isolate_bug"
    COUNTEREXAMPLE = "seek_counterexample"
    PROVE = "prove_encoded_statement"
    COMBINE = "combine_methods"
    TRANSFER = "test_transfer"
    SUSPEND = "suspend_branch"


@dataclass(frozen=True)
class Method:
    name: str
    source_sha256: str
    inputs: dict
    outputs: dict
    assumptions: tuple
    domain: tuple
    evidence: tuple
    cost_regime: dict
    dependencies: tuple
    failure_cases: tuple
    usage: str
    archive: str  # reliable, enabling, high_risk

    def __post_init__(self):
        if self.archive not in ("reliable", "enabling", "high_risk") or not self.evidence:
            raise ValueError("Methods require typed archives and evidence")


def environment(inputs, evaluator_files, library, model="gpt-6-astra", effort="high"):
    packages = {}
    for name in ("numpy", "scipy", "python-flint", "cvxpy", "quimb", "pyscf", "pytest"):
        try:
            packages[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            packages[name] = None
    value = {"inputs": inputs, "evaluator": {str(p): hashlib.sha256(Path(p).read_bytes()).hexdigest()
             for p in evaluator_files}, "library": library, "python": sys.version,
             "packages": packages, "hardware": {"platform": platform.platform(),
             "machine": platform.machine(), "processor": platform.processor()},
             "model": model, "provider": "native_codex_subscription", "reasoning_effort": effort,
             "cache_protocol": "fresh candidate cache per timed cold run; explicit warm run separately",
             "threads": 1}
    binary = Path("/Applications/ChatGPT.app/Contents/Resources/codex")
    value["native_client_sha256"] = hashlib.sha256(binary.read_bytes()).hexdigest() if binary.exists() else None
    return {"fingerprint": digest(value), **value}


class Research:
    def __init__(self, ledger, target, env):
        self.ledger, self.target, self.env = ledger, digest(target), env["fingerprint"]
        self.nodes = {}
        self.methods = {}
        ledger.append("research_environment", target=target, target_sha256=self.target, environment=env)

    def observe(self, *, action, state, claim, dependencies=(), reads=(), evidence=None, cost=None):
        dependencies = tuple(dict.fromkeys((*dependencies, *reads)))
        if any(d not in self.nodes for d in dependencies):
            raise ValueError("Cannot read or depend on an unavailable research event")
        key = self.ledger.append("research_node", target_sha256=self.target, environment_sha256=self.env,
            action=Action(action).value, state=State(state).value, claim=claim,
            dependencies=dependencies, reads=tuple(reads), evidence=evidence or {}, cost=cost or {})
        self.nodes[key] = self.ledger.events("research_node")[-1]
        return key

    def retain(self, method):
        if any(d not in self.nodes for d in (*method.evidence, *method.dependencies)):
            raise ValueError("Unacquired method dependency")
        if method.archive != "high_risk" and not all(self.nodes[d]["state"] == State.ACCEPTED.value for d in method.evidence):
            raise ValueError("An executable method needs accepted scoped evidence")
        self.ledger.read_blob(method.source_sha256)
        key = digest(asdict(method))
        self.ledger.append("retained_method", method_sha256=key, method=asdict(method))
        self.methods[key] = method
        return key

    def applicable(self, domain, available=None):
        available = set(self.nodes) if available is None else set(available)
        return {key: method for key, method in self.methods.items()
                if domain in method.domain and method.archive != "high_risk"
                and set((*method.dependencies, *method.evidence)) <= available}

    def activate(self, key, capabilities, available=None):
        method = self.methods[key]
        if key not in self.applicable(method.domain[0], available):
            raise ValueError("Method is not acquired in this replay")
        return compile_program(self.ledger.read_blob(method.source_sha256).decode(), capabilities)

    def replay(self, order, environment_sha256):
        acquired, result = set(), []
        for key in order:
            if key not in self.nodes:
                result.append({"id": key, "state": State.UNKNOWN.value, "needs_live_rerun": True})
                continue
            node = self.nodes[key]
            if not set(node["dependencies"]) <= acquired:
                raise ValueError("Replay would expose a dependency before acquisition")
            matching = environment_sha256 == node["environment_sha256"]
            result.append({"id": key, "state": node["state"] if matching else State.UNKNOWN.value,
                           "needs_live_rerun": not matching})
            # Only a live observation under the requested environment can unlock descendants.
            if matching:
                acquired.add(key)
        return result


def choose_action(state, archive, fixed=False):
    """Typed scheduling learned from completed traces; it never fabricates branch outcomes."""
    if fixed:
        return Action.SEARCH
    choices = {State.REPAIRABLE.value: Action.DEBUG, State.CONJECTURE.value: Action.PROVE,
               State.REFUTED.value: Action.COUNTEREXAMPLE, State.ACCEPTED.value: Action.COMBINE,
               State.CONDITIONAL.value: Action.PROVE, State.COST.value: Action.SEARCH,
               State.REPRESENTATION.value: Action.OPEN}
    action = choices.get(state, Action.SEARCH)
    eligible = [r for r in archive if r.get("before") == state and r.get("measured")]
    if eligible:
        # This changes a future live action, not the recorded result of an unvisited branch.
        by_action = {}
        for r in eligible:
            stats = by_action.setdefault(r["action"], [0, 0.0])
            stats[0] += int(r["accepted"] and r.get("improves_parent", True))
            stats[1] += r["cost_seconds"]
        action = Action(max(by_action, key=lambda a: (by_action[a][0] + 0.5) / (by_action[a][1] + 1.0)))
    return action
