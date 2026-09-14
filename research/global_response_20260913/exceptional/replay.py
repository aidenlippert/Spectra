"""Replay the exceptional high-occupation diagnostic without numerical packages."""
from fractions import Fraction
import json
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from research.composable_response_20260913.joint import load_case, check_gap

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "results/global_response_20260913/exceptional"

def replay(name):
    data, tail, _, _ = load_case(name)
    cert = json.loads((OUT / f"{name}_high3.json").read_text())["certificate"]
    # check_gap uses only exact rational arithmetic after loading the frozen inputs.
    return check_gap(data, tail, cert)

def replay_low_exception(name):
    data, tail, _, _ = load_case(name)
    cert = json.loads((OUT / f"{name}_low_exception_full_complement.json").read_text())["certificate"]
    return check_gap(data, tail, cert)

if __name__ == "__main__":
    for name in ("h6", "h8"):
        print(json.dumps({"case": name, "high3": replay(name), "low_exception": replay_low_exception(name)}, sort_keys=True))
