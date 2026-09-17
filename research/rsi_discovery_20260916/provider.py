"""Use the supported Codex CLI login, never extract authentication secrets."""
import json
from pathlib import Path
import shutil

from research.rsi_discovery_20260916.dialect import compile_candidate
from research.rsi_discovery_20260916.process import run, ROOT

SCHEMA = {"type": "object", "additionalProperties": False,
          "properties": {"hypothesis": {"type": "string"}, "source": {"type": "string"}},
          "required": ["hypothesis", "source"]}


def propose(ledger, parent_source, feedback, number=0, timeout=120):
    directory = ledger.root / f"model_{number}"
    directory.mkdir(exist_ok=False)
    schema = directory / "schema.json"
    schema.write_text(json.dumps(SCHEMA))
    answer = directory / "answer.json"
    prompt = ("Propose one cheaper molecular certificate constructor in the restricted Python dialect below. "
              "Return only the specified JSON. Do not use tools, access files or run commands. "
              "This is a scientific proposal; never alter the independent checker or claim an improvement without evidence. "
              "The program is exactly one function propose(payload), no imports, no decorators, no private names. "
              "Capabilities: rational(n,den=1) or rational(string), eigenpair(integer_matrix,den,subset_bool) returns "
              "(lowest_float_eigenvalue,float_vector), cholesky(integer_matrix,den,rational_shift,digits) returns "
              "rounded triangular integer factor. Use only basic arithmetic, lists, dictionaries, loops and standard "
              "len/range/zip/sum/int/round/str/min/max. Powers require literal integer exponents 0..16; "
              "do not use variable exponents. Each block must have alpha_count, lower_Ha as rational string, "
              "factor_denominator integer, triangular integer factor and nonzero integer upper_vector. "
              "The full fixed-N independent checker requires a positive semidefinite residual after factor subtraction "
              "and width <= 1/625 Hartree. Lowering factor precision may fail its exact diagonal-dominance check. "
              "Try a concrete algorithmic change to save arithmetic; numerical success alone is insufficient. "
              "Current executable parent:\n" + parent_source + "\nDevelopment observations:\n" + json.dumps(feedback))
    prompt_hash = ledger.blob(prompt)
    executable = shutil.which("codex")
    if executable is None:
        ledger.append("proposal_failure", reason="Codex CLI unavailable", prompt_sha256=prompt_hash)
        return None
    with ledger.measure("model_proposal", provider="codex_subscription_cli", prompt_sha256=prompt_hash) as metered:
        result = run([executable, "exec", "--ignore-user-config", "--ephemeral", "--sandbox", "read-only",
                      "--skip-git-repo-check", "--json", "--color", "never", "-C", str(directory),
                      "--output-schema", str(schema), "--output-last-message", str(answer), "-"],
                     directory / "process", timeout, stdin=prompt, limit_files=False)
        metered.update(result)
        usage = None
        for line in Path(result["stdout"]).read_text().splitlines():
            try:
                event = json.loads(line)
            except ValueError:
                continue
            if event.get("type") == "turn.completed":
                usage = event.get("usage")
        metered["usage"] = usage
        metered["stdout_sha256"] = ledger.blob(Path(result["stdout"]).read_bytes())
        metered["stderr_sha256"] = ledger.blob(Path(result["stderr"]).read_bytes())
        if result["exit_code"] or not answer.exists():
            metered["status"] = "failed"
            return None
        try:
            value = json.loads(answer.read_text())
            compile_candidate(value["source"], {})
        except (ValueError, SyntaxError, KeyError) as error:
            metered.update(status="failed", error=str(error))
            return None
        metered["answer_sha256"] = ledger.blob(answer.read_bytes())
        value["source_sha256"] = ledger.blob(value["source"])
        value["cost_seconds"] = result["wall_seconds"]
        ledger.append("model_proposal", prompt_sha256=prompt_hash, **value)
        return value
