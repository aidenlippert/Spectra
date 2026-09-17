"""Summarize complete receipts; missing results and token counts remain unknown."""
import json
from pathlib import Path
import statistics
import sys

from research.rsi_discovery_20260916.ledger import Ledger, digest


def read(path):
    return json.loads(path.read_text())


def is_offloaded(path):
    # Darwin sys/stat.h: SF_DATALESS. Reading such a file can block on FileProvider.
    return bool(getattr(path.stat(), "st_flags", 0) & 0x40000000)


def audit_available_evidence(ledger):
    """Audit the chain and resident objects without claiming skipped bytes passed."""
    previous = "0" * 64
    for kind, record, parent, key in ledger.db.execute(
            "SELECT kind,record,previous,hash FROM events ORDER BY seq"):
        if parent != previous or key != digest([kind, json.loads(record), parent]):
            raise ValueError("Broken history chain")
        previous = key
    verified = 0
    unavailable = []
    for path in sorted((ledger.root / "objects").iterdir()):
        if is_offloaded(path):
            unavailable.append({"object": path.name, "reason": "macOS SF_DATALESS; bytes not re-read"})
            continue
        ledger.read_blob(path.name)
        verified += 1
    starts = {e["event_hash"] for e in ledger.events("operation_started")}
    ends = {e["operation"] for e in ledger.events("operation_finished")}
    attempts = {e["id"] for e in ledger.events("attempt_started")}
    completed = {e["id"] for e in ledger.events("attempt_finished")}
    return {"chain_valid": True, "all_blob_bytes_verified": not unavailable,
            "blobs_verified": verified, "unavailable_blobs": unavailable,
            "unfinished_operations": sorted(starts - ends),
            "unfinished_attempts": sorted(attempts - completed),
            "events": ledger.db.execute("SELECT count(*) FROM events").fetchone()[0]}


def main():
    root = Path(sys.argv[1]).resolve()
    result_file = root / "confirmatory001/result.json"
    experiment = read(result_file) if result_file.exists() else None
    followup_path = root / "compound_confirmatory002/result.json"
    followup = read(followup_path) if followup_path.exists() else None
    frontier_path = root / "frontier_search003/result.json"
    frontier = read(frontier_path) if frontier_path.exists() else None
    correlation_path = root / "correlation_search003/result.json"
    correlation = read(correlation_path) if correlation_path.exists() else None
    operations, usages, audits = [], [], []
    for path in sorted(root.rglob("history.sqlite")):
        ledger = Ledger(path.parent)
        try:
            audit = audit_available_evidence(ledger)
            audits.append({"directory": str(path.parent.relative_to(root)), **audit})
            rows = ledger.events("operation_finished")
        finally:
            ledger.db.close()
        for row in rows:
            row["ledger_directory"] = str(path.parent.relative_to(root))
        operations.extend(rows)
        if len(audits) % 25 == 0:
            print(f"Audited {len(audits)} ledger chains and their available objects", flush=True)
    model_rows = [r for r in operations if r["label"] == "model_proposal"]
    for row in model_rows:
        usages.append(row["result"].get("usage"))
    accounting = {"metered_operation_wall_seconds_sum": sum(r["wall_seconds"] for r in operations),
        "wall_sum_is_not_elapsed_time": "Native inference overlaps; this is a sum of charged operation walls.",
        "self_cpu_seconds": sum(r["self_cpu_seconds"] for r in operations),
        "child_cpu_seconds": sum(r["child_cpu_seconds"] for r in operations),
        "largest_recorded_process_family_rss_bytes": max((r["result"].get("process_family_peak_rss_bytes", 0) for r in operations), default=0),
        "memory_scope": "Local child-family high-water marks, not additive. macOS does not provide the Linux address-space limit; provider accelerator memory/work is unavailable.",
        "operations": len(operations), "failed_or_timed_out_operations": sum(r["result"]["status"] in ("failed", "timeout") for r in operations),
        "native_model_calls": len(model_rows), "model_calls_with_unavailable_token_usage": sum(u is None for u in usages),
        "known_input_tokens": sum(u.get("input_tokens", 0) for u in usages if u),
        "known_output_tokens": sum(u.get("output_tokens", 0) for u in usages if u),
        "reported_reasoning_tokens_not_added_again_to_output": sum(u.get("reasoning_output_tokens", 0) for u in usages if u),
        "money_purchased_by_runner": 0, "subscription_dollar_allocation": None,
        "complete_lifetime_cost_known": False,
        "exclusions": ["parent implementation/reasoning tokens unavailable", "initial dependency installation and unmetered development commands",
            "earlier v1 development", "prototype validations outside a ledger", "report maintenance", "unknown token usage on interrupted model calls"],
        "audit_scope": "All event hash chains and all resident content-addressed objects. Offloaded objects are enumerated as unavailable, never counted as byte-verified.",
        "freshly_verified_blobs": sum(a["blobs_verified"] for a in audits),
        "unavailable_blobs": sum(len(a["unavailable_blobs"]) for a in audits),
        "unfinished_operations": sum(len(a["unfinished_operations"]) for a in audits),
        "unfinished_attempts": sum(len(a["unfinished_attempts"]) for a in audits),
        "audits": audits}
    accounting["study_costs"] = {name: {
        "metered_operation_wall_seconds_sum": sum(r["wall_seconds"] for r in operations if r["ledger_directory"].startswith(name + "/") or r["ledger_directory"] == name),
        "native_model_calls": sum(r["ledger_directory"].startswith(name + "/") or r["ledger_directory"] == name for r in model_rows)}
        for name in ("confirmatory001", "compound_confirmatory001", "compound_confirmatory002", "later_process_conformance001",
                     "frontier_search001", "frontier_search002", "frontier_search003", "correlation_search001",
                     "correlation_search002", "correlation_search003", "runtime_recovery001", "runtime_recovery002")}
    (root / "accounting.json").write_text(json.dumps(accounting, indent=2) + "\n")
    physical = []
    for name in ("nonenumerating006", "nonenumerating_reference001", "physical_transfer_water001", "physical_transfer_h6001"):
        path = root / name / "result.json"
        if path.exists():
            record = read(path)
            physical.append({"run": name, "status": record["status"],
                "width_mHa": record["acceptance"].get("width_mHa"),
                "target_met": record["acceptance"].get("target_1p6mHa_met", False),
                "metered_seconds": record["costs"]["metered_wall_seconds"]})
    incremental = []
    if experiment:
        for p in sorted((root / "confirmatory001").glob("campaign_*/result.json")):
            c = read(p); directory = p.parent.relative_to(root)
            evaluation = [r for r in operations if r["ledger_directory"] == str(directory)]
            def eval_cost(labels):
                return sum(r["wall_seconds"] for r in evaluation if r.get("metadata", {}).get("label_id") in labels)
            acquisition = c["acquisition_cost"]["metered_wall_seconds"] + eval_cost(["acquisition"])
            dev = {a: sum(t["cost"]["metered_wall_seconds"] for t in item["attempts"]) + eval_cost([a + "_0", a + "_1"])
                   for a, item in c["arms"].items()}
            future = {a: eval_cost(["transfer_" + a]) for a in c["arms"]}
            incremental.append({"campaign": c["index"] + 1, "acquisition_seconds": acquisition,
                "discovery_seconds": dev, "complete_recorded_transfer_seconds": future,
                "D_minus_A_upfront_seconds_excluding_shared_setup": acquisition + dev["D"] - dev["A"],
                "D_net_saved_seconds_after_recorded_transfer_excluding_shared_setup":
                    future["A"] - future["D"] - acquisition - dev["D"] + dev["A"],
                "scope": "Conditional accounting comparison only. Shared policy/pilot development and complete lifetime costs are additional."})
    (root / "incremental_costs.json").write_text(json.dumps(incremental, indent=2) + "\n")
    followup_incremental = None
    if followup:
        acquisition = sum(accounting["study_costs"][name]["metered_operation_wall_seconds_sum"]
                          for name in ("confirmatory001", "later_process_conformance001", "compound_confirmatory001"))
        reuse = followup["assigned_backend"]["comparisons"]
        comparable = all(row["comparable"] for row in reuse)
        saved = sum(row["construction_seconds"]["A"] - row["construction_seconds"]["D"] for row in reuse) if comparable else None
        followup_incremental = {"known_preceding_acquisition_cost_seconds": acquisition,
            "current_followup_recorded_research_cost_seconds": accounting["study_costs"]["compound_confirmatory002"]["metered_operation_wall_seconds_sum"],
            "single_cold_construction_per_heldout_task_saved_seconds": saved,
            "construction_savings_minus_preceding_acquisition_seconds": saved - acquisition if saved is not None else None,
            "complete_lifetime_payback_proved": False,
            "scope": "Charged operation-wall convention. One median cold construction per held-out task; this does not multiply the savings by the three timing repetitions. New research, activation, verification and maintenance are additional costs, not savings. Parent/setup/provider accelerator work and interrupted-call tokens remain unavailable."}
        (root / "followup_incremental_costs.json").write_text(json.dumps(followup_incremental, indent=2) + "\n")
    passed = bool(experiment and experiment["predeclared_operational_test_passed"])
    followup_passed = bool(followup and followup["predeclared_operational_test_passed"])
    lines = ["# Spectra RSI discovery v2 — evidence receipt", "",
        ("A scoped discovery-improvement test passed; the exact scope and negative results are below." if passed or followup_passed else
         "Neither completed scoped discovery-improvement study passed." if experiment and followup else
         "The first scoped discovery test failed. The corrected follow-up is still running." if experiment else "The confirmatory experiment is still running."),
        "This is not a proof of universal or indefinite recursive self-improvement. Complete lifetime payback is not established.", "",
        "The implementation now supports executable algebra programs, operational method packages, typed research states, full recorded dependencies, live reruns under changed environments, and separate reliable/enabling/high-risk archives. It connects generated methods to the existing non-enumerating MPS/SOS pipeline and preserves its original exact acceptance path.", "",
        "## The tested chain", "",
        "1. Generation one independently proposes and checks a general CAR compiler in each campaign.",
        "2. Generation two discovers executable symmetric Gram-map algorithms using a frozen or acquired library and a fixed or replay-informed policy. These are programs, not precision parameters.",
        "3. Every selected program is frozen before transfer to unseen support graphs, larger dictionaries, adjoints and injective mode relabelings. Every transfer map is compared with the original exact semantics.", "",
        "A = frozen library/fixed policy; B = acquired library/fixed policy; C = frozen library/replay-informed policy; D = acquired library/replay-informed policy. All use the same native GPT-6-Astra/high backend, task inputs and per-call resource limits. Five independent acquisition/search campaigns are the units of comparison; candidates and repeated timing samples are not independent campaigns.", ""]
    if experiment:
        lines += ["## Frozen four-arm experiment", "", "| Campaign | A/D speedup | A/B | C/D | Comparable |", "|---|---:|---:|---:|---|"]
        for r in experiment["comparisons"]:
            s = r["speedups"]
            lines.append(f"| {r['campaign']} | {s.get('A/D', float('nan')):.3f} | {s.get('A/B', float('nan')):.3f} | {s.get('C/D', float('nan')):.3f} | {r['comparable']} |")
        lines += ["", "A ratio above one favors the denominator. These ratios measure cold construction of exactly equal maps, including fresh method caches, not whole molecular solution time.", "",
            f"One-sided paired sign-test p: `{experiment['paired_sign_test_one_sided_p']}`. Predeclared operational gate passed: `{passed}`. Source snapshot unchanged during the run: `{experiment['source_snapshot_unchanged']}`.", "",
            "The primary gate required five D/A wins, median speedup at least 1.2, exact D transfer, a new generated D method in every campaign, and median A/B and C/D construction-time ratios at least 1.2. Secondary policy contrasts are descriptive. Shared model weights, related algebra tasks and one machine limit generalization; this is not a population-level claim about all scientific discovery.", ""]
    conformance_path = root / "later_process_conformance001/result.json"
    if conformance_path.exists():
        conformance = read(conformance_path)
        accepted = [r for r in conformance["results"] if r["status"] == "verified_encoded_claim"]
        lines += ["## Later-process reuse", "",
            f"{len(accepted)} retained methods reopened from their original ledgers and passed all 146 actual clustered-frame operator groups at 8 and 12 modes. Three arms had no retained generated method and remained unobserved. This is additional semantic conformance, not 146 independent research campaigns. The protocol was frozen before this conformance run, after the original study had entered transfer.", ""]
    lines += ["## New representation problem with a fair baseline", "",
        "The follow-up retains the fastest exact monomial-map program from the first study (campaign 2, arm A) and asks six fresh discovery campaigns to invent maps for sparse polynomial operators. Thus the acquired primitive itself came from a frozen-library control. This is conditional on one shared acquired library, not six new independent acquisitions.", "",
        "Every arm receives the existing ordinary monomial-map interface. B/D receive the acquired backend through that same interface. After selection, every program is frozen and tested on new support families, zero/dependent operators and integers beyond floating-point precision. A second transfer executes every program with the same ordinary backend; inherited kernel execution speed alone cannot pass the discovery gate.", "",
        "The initial follow-up was invalidated before any development evaluation or transfer because its controls omitted the existing monomial-map capability. Its 12 started model calls remain charged. The corrected protocol requires six D/A wins, raw median A/D≥1.25, common-backend median A/D≥1.2, A/B and C/D≥1.2 under both backends, exact outputs and actual acquisition consumption. The two completed RSI studies form a Bonferroni family with threshold 0.025 for the follow-up; the discarded design supplies no hypothesis-test evidence.", ""]
    if followup:
        lines += ["| Campaign | Assigned-backend A/D | Common-backend A/D | Assigned A/B | Common A/B |",
                  "|---|---:|---:|---:|---:|"]
        for raw, common in zip(followup["assigned_backend"]["comparisons"], followup["common_ordinary_backend"]["comparisons"]):
            a, b = raw["speedups"], common["speedups"]
            lines.append(f"| {raw['campaign']} | {a.get('A/D', float('nan')):.3f} | {b.get('A/D', float('nan')):.3f} | {a.get('A/B', float('nan')):.3f} | {b.get('A/B', float('nan')):.3f} |")
        lines += ["", f"Corrected operational gate passed: `{followup_passed}`. Assigned-backend D/A wins: `{followup['assigned_backend']['D_over_A_wins']}/6`; common-backend wins: `{followup['common_ordinary_backend']['D_over_A_wins']}/6`. Intersection-union one-sided p: `{followup['intersection_union_p_one_sided']}`; two-study adjusted p: `{followup['two_study_Bonferroni_p']}`. Protocol source unchanged during execution: `{followup['source_snapshot_unchanged']}`.", ""]
        lines += [f"The raw median A/D ratio was {followup['assigned_backend']['median_speedups']['A/D']:.3f}; the common-backend ratio was {followup['common_ordinary_backend']['median_speedups']['A/D']:.3f}. One D program independently implemented its own compiler and made zero calls to the acquired primitive, failing the additional consumption requirement. These are useful execution results but do not pass the discovery-improvement criterion.", "",
            "The replay-informed policy selected the same continue-search action as the fixed policy in every second-round follow-up attempt. Thus this study does not demonstrate a policy contribution. The first study did exercise different repair actions. The separate policy-action audit preserves those distinctions.", ""]
    else:
        lines += ["Corrected discovery/transfer is pending; no positive finding is inferred from development timings.", ""]
    lines += ["## Actual molecular workflow", "", "| Run | Exact width (mHa) | ≤1.6 mHa | Metered seconds |", "|---|---:|---|---:|"]
    for r in physical:
        width = "unknown" if r["width_mHa"] is None else f"{r['width_mHa']:.6f}"
        lines.append(f"| {r['run']} | {width} | {r['target_met']} | {r['metered_seconds']:.3f} |")
    lines += ["", "The full fixed-N gate combines separately checked singlet and all-nonsinglet lower certificates and charges the original-H spin defect. Upper bounds come from exact charge-MPS contractions. These runs construct no full molecular occupation sector and read no ground-state or trajectory teacher. The H4 reference and generated-map runs produce identical certified endpoints. Their similar total times do not establish whole-pipeline acceleration.", "",
        "The water/H6 outcomes are separate tests of physical certificate quality. An exact coefficient map can transfer successfully while a bounded MPS/SOS search fails its chemical-accuracy target. Such a failure is not a proof of a representational impossibility.", "",
        "## Mathematical and frontier gates", "",
        "- Three explicit operator dictionaries reduce from three operators to two with all six Gram cross-term witnesses reconstructed in the declared number ideal. An intentionally narrower ideal refuses optimization equivalence despite a valid physical null. Scalar energy columns and column sets spanning identity are rejected.",
        "- Five Lean transport/cross-term/error lemmas pass with foundational axioms only. The full Python programs are not formally verified.",
        "- A new checked-template operation verifies each novel local algebra template once, then transports its immutable witness by increasing mode relabeling. Exact comparisons on real coefficient maps include every cache-construction/checking cost.",
        "- On the original amplitude-0.5, four-phase H8 control task, the four-observable baseline computes exact defining dynamical defects and rigorous propagation/rounding bounds without a full-state teacher. Its enclosure is only [-2,2], insufficient for D≤-0.6. The target and uncertainty budgets are unchanged; the subsequent executable representation search is reported below.", "",
        "See [DERIVATIONS.md](" + str(Path(__file__).with_name("DERIVATIONS.md")) + ") for the encoded statements and mathematical assumptions.", "",
        "## Teacher-free representation search", "",
        "Candidate programs can now submit up to fourteen additional sparse Hermitian observables. Identity and the target D are prepended by the checker. Exact matrices and residuals are derived from the unchanged H8 inputs, in a separate process that never executes proposal code. A further fresh exact replay checks the certificate. Every unsuccessful representation remains a typed high-risk method with its unresolved obligation; it cannot be activated as an accepted capability.", ""]
    if frontier:
        lines += ["The first two native control proposals exceeded their 150-second client deadline. One completed response contained executable code in the mechanism field and prose in the source field. Its code was recovered verbatim, with the original response and transformation preserved; it was not retroactively counted as within budget. The source schema was clarified for subsequent model calls.", "",
            "The recovered algorithm constructed a larger observable basis, but independent exact acceptance exceeded its 180-second deadline. A second live run reduced the allowed extra observables to four without changing the algorithm, Hamiltonian, controls, target or uncertainty. That exact check also exceeded the budget. Neither unfinished check supplies a certified interval or proves impossibility.", ""]
        for i, row in enumerate(frontier["proposals"]):
            lines.append(f"- Candidate {i + 1}: `{row['status']}`; interval `{row.get('D_interval')}`; dimensions `{row.get('dimensions')}`; error `{row.get('error')}`.")
        lines += ["", f"Original robust target proved by this exploratory representation search: `{frontier['target_proved']}`. This is not a controlled RSI comparison.", ""]
    else:
        lines += ["The candidate interface has passed separate-process integration on a five-observable program. Native representation proposals are pending; the integration example is not counted as a learned discovery.", ""]
    lines += ["## Executable correlation selection", "",
        "A native proposal chooses local and coupled collective orbital supports directly from Hamiltonian coefficients before constructing cubic maps. The support count matches the width-three window baseline. The retained exact map compiler is consumed by the live MPS/SOS workflow; the fixed original-H full-N checker determines each interval. One source is frozen before both H4 development and H6 transfer, with no transfer repair. This is a small physical exploration, not another RSI significance test.", ""]
    if correlation:
        if correlation.get("outcomes"):
            lines += ["| Case | Actual width (mHa) | State | Workflow seconds |", "|---|---:|---|---:|"]
            for row in correlation["outcomes"]:
                width = row["acceptance"].get("width_mHa")
                cost = row.get("costs") or {}
                display = "unknown" if width is None else f"{width:.6f}"
                lines.append(f"| {row['case']} | {display} | {row['status']} | {cost.get('metered_wall_seconds', float('nan')):.3f} |")
            lines += ["", "The unchanged program selected local two-orbital supports while keeping the existing coupled collective-pair channels. H6's accepted interval narrowed from 287.010 to 10.685 mHa at the same declared solver budget; it still misses the 1.6 mHa target. Recorded workflow times were 158.379 seconds for the earlier window baseline and 162.215 seconds for the new selection. H4 remained chemically accurate but did not improve the earlier interval. These are single exploratory runs on a shared machine, not a demonstrated speed advantage or causal RSI result.", "",
                "An initial historical-path lookup failed before model dispatch. The next run produced the frozen program, but both physical executions hit environment-loading deadlines while numerical-library files in Documents were offloaded. Those failures stay charged. The final physical results rerun the same source after restoring all 34 pinned package versions in a local runtime outside Documents; no transfer feedback was used to edit the program.", ""]
        else:
            lines += [f"The native correlation proposal did not complete: `{correlation.get('status')}`. No physical capability is inferred.", ""]
    else:
        lines += ["Native proposal and physical execution are pending.", ""]
    lines += ["## Costs and limitations", "",
        f"Recorded operation walls sum to {accounting['metered_operation_wall_seconds_sum']:.3f} seconds, including {accounting['native_model_calls']} native model calls and {accounting['failed_or_timed_out_operations']} failed/timed-out operations. Concurrent inference means this sum is not elapsed wall time. Known completed-call usage is {accounting['known_input_tokens']:,} input tokens and {accounting['known_output_tokens']:,} output tokens; {accounting['model_calls_with_unavailable_token_usage']} calls have unavailable usage.", "",
        "The full raw costs, failed attempts, exact outputs and hash-chain audits remain in the result directories. `incremental_costs.json` separates per-campaign discovery/acquisition from later transfer and reports a conditional net cost. Shared policy training and prototype costs are additional. Parent implementation tokens, some initial setup, report maintenance and interrupted-call token counts are unavailable, so a complete lifetime efficiency or dollar-payback claim would be unsupported.", "",
        f"The final audit verified all {len(audits)} ledger hash chains and re-read {accounting['freshly_verified_blobs']} content-addressed objects. It could not re-read {accounting['unavailable_blobs']} objects marked offloaded by macOS; those contents are not counted as freshly verified. Unfinished recorded operations: {accounting['unfinished_operations']}; unfinished attempts: {accounting['unfinished_attempts']}. Archived study-source snapshots and all 98 copied Spectra dependency modules were separately re-read and matched their recorded hashes.", "",
        "Processor work, model tokens and local process-family memory are reported separately in accounting.json. Local memory readings are high-water marks, not additive allocation totals or provider accelerator usage. The second study charges the entire first study and later conformance as acquisition, rather than charging only its selected successful proposal.", "",
        "Validated domain coverage is limited to the executed algebra, representation, molecular-bound and control-defect tasks. No inverse-design system, laboratory adapter, or experimental-science capability is claimed. Negative finite studies do not prove that recursive improvement is impossible.", "",
        "Final regression validation passed 60 tests; four additional reporting-audit tests passed after handling offloaded evidence. The restored runtime imports the required numerical libraries and Z3 solves the checked integer example. The installer reports a Z3 wheel-platform-tag warning (macosx_13_3_arm64 on macOS 26.4.1 arm64); its native runtime check and the regression suite passed. This packaging warning is preserved in verification.json rather than silently reported as a clean installer check.", "",
        "No purchased cloud compute, API-key extraction, laboratory action, or GitHub push was used. The original dirty Spectra checkout was preserved; work is in the isolated `codex/rsi-discovery` checkout.", "",
        "## Receipts", "",
        f"- [Protocol]({root / 'confirmatory001/protocol.json'})",
        f"- [Four-arm result]({root / 'confirmatory001/result.json'})",
        f"- [Corrected follow-up protocol]({root / 'compound_confirmatory002/protocol.json'})",
        f"- [Corrected follow-up result]({followup_path})",
        f"- [Invalidated design receipt]({root / 'compound_confirmatory001/design_invalidation.json'})",
        f"- [Later-process conformance]({conformance_path})",
        f"- [Accounting and audits]({root / 'accounting.json'})",
        f"- [Incremental cost comparison]({root / 'incremental_costs.json'})",
        f"- [Follow-up acquisition and reuse costs]({root / 'followup_incremental_costs.json'})",
        f"- [Full-N molecular acceptance]({root / 'nonenumerating006/acceptance.json'})",
        f"- [Exact dictionary witnesses]({root / 'quotient004/result.json'})",
        f"- [Checked operation and Lean receipts]({root / 'proof_operations001/result.json'})",
        f"- [Original control-task result]({root / 'frontier003/result.json'})",
        f"- [Executable representation search]({frontier_path})",
        f"- [Direct correlation-selection run]({correlation_path})",
        f"- [Initial control deadlines]({root / 'frontier_search001/result.json'})",
        f"- [Recovered control-program evaluation]({root / 'frontier_search002/result.json'})",
        f"- [Suspended control obligation]({root / 'frontier_search003/suspension.json'})",
        f"- [Policy-action audit]({root / 'policy_action_audit.json'})",
        f"- [Final verification]({root / 'verification.json'})",
        f"- [Final audit recovery and limitation]({root / 'final_audit_recovery.json'})",
        f"- [Source and dependency integrity]({root / 'source_integrity.json'})", ""]
    report = Path(__file__).with_name("RESULTS.md")
    report.write_text("\n".join(lines))
    print(json.dumps({"report": str(report), "confirmatory_complete": experiment is not None,
                      "first_test_passed": passed, "followup_test_passed": followup_passed,
                      "physical": physical, "accounting": {k: v for k, v in accounting.items() if k != "audits"}}, indent=2))


if __name__ == "__main__":
    main()
