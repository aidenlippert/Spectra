"""Causal replay of recorded branches. Unrecorded actions have no outcome.

Policies see only visible action metadata and outcomes they already paid for.
Replay includes proposal, failed evaluation and checking costs stored per node.
Historical selection is explicitly separate from later, fresh evaluations.
"""
from dataclasses import asdict, dataclass
from itertools import product
import math
import time


@dataclass(frozen=True)
class Policy:
    name: str = "fixed"
    branch_width: int = 1
    depth_first: bool = False
    patience: int = 999
    max_attempts: int = 20
    adaptive: bool = False
    family_order: tuple = ()

    def choose(self, frontier, observed):
        rewards, counts = {}, {}
        for node in observed:
            family = node["family"]
            rewards[family] = rewards.get(family, 0) + node["utility"]
            counts[family] = counts.get(family, 0) + 1

        def key(action):
            family = action["family"]
            order = self.family_order.index(family) if family in self.family_order else len(self.family_order)
            value = 0
            if self.adaptive:
                n = counts.get(family, 0)
                value = (rewards.get(family, 0) / n + math.sqrt(2 * math.log(2 + len(observed)) / n)) if n else float("inf")
            return (order, -value, -action["depth"] if self.depth_first else action["depth"], action["id"])

        return sorted(frontier, key=key)[:self.branch_width]


def replay(nodes, policy, budget_seconds):
    if budget_seconds <= 0 or policy.branch_width < 1 or policy.max_attempts < 1 or policy.patience < 1:
        raise ValueError("Positive bounded replay allocations required")
    by_id = {n["id"]: n for n in nodes}
    if len(by_id) != len(nodes):
        raise ValueError("Duplicate recorded action identity")
    for n in nodes:
        if n["parent"] and n["parent"] not in by_id:
            raise ValueError("Replay cannot synthesize missing ancestors")
        if not math.isfinite(n["cost_seconds"]) or n["cost_seconds"] < 0 or not math.isfinite(n["utility"]):
            raise ValueError("Invalid recorded cost or outcome")
        seen = {n["id"]}
        parent = n["parent"]
        while parent:
            if parent in seen:
                raise ValueError("Cyclic discovery history")
            seen.add(parent)
            parent = by_id[parent]["parent"]
    observed, completed = [], set()
    cost, stale, best = 0.0, 0, 0.0
    while len(observed) < policy.max_attempts and cost < budget_seconds and stale < policy.patience:
        frontier = [{k: n[k] for k in ("id", "parent", "family", "depth")}
                    for n in nodes if n["id"] not in completed and (n["parent"] is None or n["parent"] in completed)]
        if not frontier:
            break
        chosen = policy.choose(frontier, list(observed))
        allowed = {n["id"] for n in frontier}
        if not chosen or len({n["id"] for n in chosen}) != len(chosen) or any(n["id"] not in allowed for n in chosen):
            raise ValueError("Unsupported replay action; no invented outcome")
        # Branch groups are serial on the declared one-process local resource.
        for action in chosen:
            if len(observed) >= policy.max_attempts or cost >= budget_seconds or stale >= policy.patience:
                break
            node = by_id[action["id"]]
            cost += node["cost_seconds"]
            observed.append(node)
            completed.add(node["id"])
            if node["utility"] > best:
                best, stale = node["utility"], 0
            else:
                stale += 1
    return {"best_utility": best, "charged_seconds": cost, "visited": [n["id"] for n in observed],
            "overshoot_seconds": max(0, cost - budget_seconds),
            "scope": "Recorded outcomes only; no prediction for untried discoveries"}


def train(ledger, histories, budget_seconds):
    with ledger.measure("policy_development", replay_histories=len(histories)) as metered:
        families = sorted({n["family"] for nodes in histories for n in nodes})
        candidates = [Policy(), Policy(name="adaptive", adaptive=True)]
        for width, depth, patience, preferred in product((1, 2), (False, True), (2, 4), ((), *[(f,) for f in families])):
            candidates.append(Policy(name="replay_selected", branch_width=width, depth_first=depth,
                                     patience=patience, family_order=preferred))
        records = []
        for candidate in candidates:
            results = [replay(nodes, candidate, budget_seconds) for nodes in histories]
            score = sum(r["best_utility"] for r in results)
            spent = sum(r["charged_seconds"] for r in results)
            records.append({"policy": asdict(candidate), "utility": score, "charged_seconds": spent, "replays": results})
        best = min(records, key=lambda r: (-r["utility"], r["charged_seconds"]))
        metered.update(policies_evaluated=len(records), recorded_actions=sum(map(len, histories)))
    key = ledger.blob(__import__("json").dumps(records, sort_keys=True))
    ledger.append("policy_selection", policy=best["policy"], evidence_sha256=key,
                  development_task_ids=sorted({n["task_id"] for nodes in histories for n in nodes}),
                  scope="Replay optimization on development histories; fresh validation required")
    policy = dict(best["policy"])
    policy["family_order"] = tuple(policy["family_order"])
    return Policy(**policy), records
