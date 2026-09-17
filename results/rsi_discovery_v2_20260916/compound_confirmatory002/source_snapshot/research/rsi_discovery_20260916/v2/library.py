"""Reopen retained methods in a later process without hand-editing its prompt.

Activation returns a proposer, never an acceptance decision. A new environment
still needs its own live semantic check. Historical outcomes are evidence about
their original inputs, not counterfactual outcomes on the new task.
"""
from research.rsi_discovery_20260916.ledger import digest
from research.rsi_discovery_20260916.v2.language import compile_program


class MethodLibrary:
    def __init__(self, ledger, through=None):
        self.ledger = ledger
        ledger.audit()
        self.nodes, self.methods = {}, {}
        seen = set()
        found_cursor = through is None
        for event in ledger.events():
            if event["kind"] == "research_node":
                if not set(event["dependencies"]) <= self.nodes.keys():
                    raise ValueError("Retained history exposes a dependency before acquisition")
                self.nodes[event["event_hash"]] = event
            elif event["kind"] == "retained_method" and "method_sha256" in event:
                method = event["method"]
                if digest(method) != event["method_sha256"]:
                    raise ValueError("Retained manifest changed")
                if not set(method["dependencies"] + method["evidence"]) <= self.nodes.keys():
                    raise ValueError("Retained method has an unavailable dependency")
                ledger.read_blob(method["source_sha256"])
                if method["archive"] != "high_risk" and any(self.nodes[e]["state"] != "verified_encoded_claim" for e in method["evidence"]):
                    raise ValueError("Callable method lacks accepted scoped evidence")
                self.methods[event["method_sha256"]] = method
            seen.add(event["event_hash"])
            if event["event_hash"] == through:
                found_cursor = True
                break
        if not found_cursor:
            raise ValueError("Unknown historical horizon")

    def applicable(self, domain):
        return {key: value for key, value in self.methods.items()
                if domain in value["domain"] and value["archive"] != "high_risk"}

    def ancestors(self, ids):
        result, pending = set(), list(ids)
        while pending:
            key = pending.pop()
            if key in result:
                continue
            if key not in self.nodes:
                raise ValueError("Unknown dependency")
            result.add(key)
            pending.extend(self.nodes[key]["dependencies"])
        return result

    def activate(self, key, base_capabilities):
        method = self.methods[key]
        if method["archive"] == "high_risk":
            raise ValueError("An unresolved branch is not an accepted executable capability")
        allowed = self.ancestors(method["dependencies"] + method["evidence"])
        scope = dict(base_capabilities)
        # Dependencies are activated in their original acquisition order, never by name alone.
        for other_key, other in self.methods.items():
            if other_key == key:
                break
            if other["archive"] != "high_risk" and set(other["dependencies"] + other["evidence"]) <= allowed:
                if other["name"] in scope:
                    raise ValueError("Retained capability cannot shadow a base operation")
                scope[other["name"]] = compile_program(self.ledger.read_blob(other["source_sha256"]).decode(), scope)
        return compile_program(self.ledger.read_blob(method["source_sha256"]).decode(), scope)
