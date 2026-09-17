"""Conjecture, exact counterexample, finite CAR proof and Lean proof records."""
from fractions import Fraction as F
import json
from pathlib import Path
import re
import shutil

from experiments.marginal_symbolic import add, canonical, mono, product, scale, encode
from research.rsi_discovery_20260916.candidates import LEAN_PAIR, LEAN_COMPOSITION
from research.rsi_discovery_20260916.process import run


def cross_term_counterexample(coefficient, timeout_ms=3000):
    import z3
    x, y = z3.Reals("x y")
    s = z3.Solver()
    s.set(timeout=timeout_ms)
    s.add(coefficient * x * y > x * x + y * y)
    answer = s.check()
    if answer == z3.sat:
        model = s.model()
        a, b = (F(str(model.eval(v))) for v in (x, y))
        if not coefficient * a * b > a * a + b * b:
            raise ValueError("Counterexample did not replay exactly")
        return {"status": "refuted", "witness": {"x": str(a), "y": str(b)},
                "domain": "rational scalars", "coefficient": coefficient}
    return {"status": "conjecture", "solver_result": str(answer),
            "domain": "real scalars", "coefficient": coefficient,
            "note": "A solver result alone is not promoted to a checked general lemma"}


def car_lemma():
    n0, n1 = (mono(((1, i), (0, i))) for i in (0, 1))
    pair = product(n0, n1)
    witness = add(product(pair, pair), scale(pair, -1))
    if canonical(witness):
        raise ValueError("Pair occupation projector identity failed")
    # The proof artifact is an exact CAR polynomial identity on these labels.
    return {"status": "verified_instance", "statement": "(n_0 n_1)^2 = n_0 n_1",
            "assumptions": ["two specified fermionic modes obey CAR"],
            "witness": {"lhs": encode(product(pair, pair)), "rhs": encode(pair)},
            "general_renaming_rule_added": False}


def discover(ledger):
    receipts = []
    for coefficient in (3, 2):
        with ledger.measure("lemma_counterexample_search", coefficient=coefficient) as cost:
            receipt = cross_term_counterexample(coefficient)
            cost.update(receipt)
        receipts.append(ledger.append("lemma", **receipt))
    with ledger.measure("finite_car_proof") as cost:
        receipt = car_lemma()
        cost.update(receipt)
    receipts.append(ledger.append("lemma", **receipt))
    lean = shutil.which("lean")
    if lean is None:
        fallback = Path.home() / ".elan/bin/lean"
        lean = str(fallback) if fallback.exists() else None
    for name, source, statement in (
        ("cross_term_bounds", LEAN_PAIR, "For all integers x,y, |2xy| <= x^2+y^2"),
        ("compose_lower_errors", LEAN_COMPOSITION, "For integers, L <= E+a+b implies L-a-b <= E"),
    ):
        source_key = ledger.blob(source)
        path = ledger.root / f"{name}.lean"
        path.write_text(source)
        if lean is None:
            receipt = {"status": "conjecture", "reason": "Lean unavailable"}
        else:
            with ledger.measure("lean_proof_check", theorem=name) as cost:
                result = run([lean, str(path)], ledger.root / f"lean_{name}", 45)
                out = Path(result["stdout"]).read_text()
                err = Path(result["stderr"]).read_text()
                cost.update(result, stdout_sha256=ledger.blob(out), stderr_sha256=ledger.blob(err))
                # No sorryAx, custom axioms, unsafe declarations or extra imports
                # are permitted for these frozen arithmetic templates.
                matches = re.findall(r"depends on axioms: \[([^\]]*)\]", out)
                axioms = {a.strip() for group in matches for a in group.split(",") if a.strip()}
                audited = bool(matches) or "does not depend on any axioms" in out
                proved = result["exit_code"] == 0 and audited and axioms <= {"propext", "Quot.sound", "Classical.choice"}
                receipt = {"status": "general_lemma" if proved else "conjecture",
                           "axiom_report": out, "allowed_foundational_axioms": sorted(axioms), "checker": lean, "process": result}
                if not proved:
                    cost["status"] = "failed"
        record = ledger.append("lemma", name=name, statement=statement, source_sha256=source_key,
                               assumptions=["Lean Int semantics"], **receipt)
        receipts.append(record)
        if receipt["status"] == "general_lemma":
            ledger.append("retained_method", category="lemma", name=name, evidence=record,
                          source_sha256=source_key, scope=statement, acquisition_costs="entire preceding lemma campaign")
    return receipts
