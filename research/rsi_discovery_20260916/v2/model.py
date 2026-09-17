"""Fixed native subscription backend, with complete prompt/output receipts."""
import json
from pathlib import Path
import shutil

from research.rsi_discovery_20260916.process import run

SCHEMA = {"type": "object", "additionalProperties": False, "properties": {
    "hypothesis": {"type": "string"}, "mechanism": {"type": "string", "description": "Plain-language algorithm summary, not executable code."},
    "source": {"type": "string", "description": "Complete executable Python source code, including def propose(payload) and every helper. This field is code, not source attribution."}, "assumptions": {"type": "string"},
    "failure_cases": {"type": "string"}},
    "required": ["hypothesis", "mechanism", "source", "assumptions", "failure_cases"]}

LANGUAGE = """Return the JSON schema. Do not use tools, access files, or execute commands.
Put complete executable pure Python in the source field, with entry point propose(payload).
Put the plain-language explanation in mechanism. Helper functions are permitted.
No imports, decorators, annotations, default arguments, private names, lambda, exceptions or I/O.
Allowed: for/while/if, comprehensions, lists, tuples, dictionaries, sets, ordinary arithmetic,
literal powers 0..16; builtins len/range/enumerate/zip/sum/min/max/abs/round/int/float/str/bool/
list/tuple/dict/set/sorted/reversed/any/all. Methods append/extend/get/items/keys/values/pop/
copy/setdefault/update/clear/add/discard/sort/index/count. Rational attributes numerator/denominator.
Your program runs under a hard time limit. Preserve exact semantics; do not specialize to supplied examples.
Do not claim measured improvement. Return one coherent new algorithm, not an evaluation harness.
"""


def propose(ledger, prompt, label, timeout=120, model="gpt-6-astra", effort="high"):
    directory = ledger.root / "models" / label
    directory.mkdir(parents=True, exist_ok=False)
    (directory / "schema.json").write_text(json.dumps(SCHEMA))
    answer = directory / "answer.json"
    full = LANGUAGE + "\n" + prompt
    prompt_hash = ledger.blob(full)
    bundled = Path("/Applications/ChatGPT.app/Contents/Resources/codex")
    executable = str(bundled) if bundled.is_file() else shutil.which("codex")
    if executable is None:
        raise RuntimeError("Native Codex backend unavailable")
    with ledger.measure("model_proposal", label_id=label, model=model, effort=effort, executable=executable,
                        prompt_sha256=prompt_hash) as receipt:
        result = run([executable, "exec", "--ignore-user-config", "--ephemeral", "--sandbox", "read-only",
                      "--skip-git-repo-check", "--json", "--color", "never", "--model", model,
                      "-c", 'model_reasoning_effort="' + effort + '"', "-C", str(directory),
                      "--output-schema", str(directory / "schema.json"),
                      "--output-last-message", str(answer), "-"], directory / "process", timeout,
                     stdin=full, limit_files=False)
        receipt.update(result)
        events = []
        for line in Path(result["stdout"]).read_text().splitlines():
            try:
                events.append(json.loads(line))
            except ValueError:
                pass
        receipt["usage"] = next((e.get("usage") for e in reversed(events) if e.get("type") == "turn.completed"), None)
        receipt["stdout_sha256"] = ledger.blob(Path(result["stdout"]).read_bytes())
        receipt["stderr_sha256"] = ledger.blob(Path(result["stderr"]).read_bytes())
        receipt["tool_actions"] = [e for e in events if e.get("item", {}).get("type") in ("command_execution", "mcp_tool_call", "web_search")]
        if result["exit_code"] or not answer.exists() or receipt["tool_actions"]:
            receipt["status"] = "timeout" if result["timeout"] else "failed"
            return None
        value = json.loads(answer.read_text())
        receipt["answer_sha256"] = ledger.blob(answer.read_bytes())
        value["source_sha256"] = ledger.blob(value["source"])
        value["model_wall_seconds"] = result["wall_seconds"]
        value["prompt_sha256"] = prompt_hash
        value["usage"] = receipt["usage"]
        ledger.append("scientific_proposal", label_id=label, **value)
        return value
