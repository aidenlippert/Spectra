"""Small, explicit Lambda Cloud pool controller.

This module does not launch or terminate at import time.  Termination accepts
only instance IDs previously recorded by this tool in the local registry.
"""
import argparse, json, os, re, sys, time
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

API = "https://cloud.lambda.ai/api/v1"
KEY = Path.home() / ".config/spectra/lambda-api-key"
ROOT = Path(__file__).resolve().parents[2]
REGISTRY = ROOT / "results/lambda_runs/direct_discovery_pool/pool.json"

def _key():
    value = KEY.read_text().strip()
    if not value: raise RuntimeError("Lambda API key file is empty")
    return value

def _json(method, path, body=None):
    payload = None if body is None else json.dumps(body).encode()
    req = Request(API + path, data=payload, method=method,
                  headers={"Authorization": "Bearer " + _key(),
                           "User-Agent": "spectra-research/1.0",
                           "Content-Type": "application/json"})
    with urlopen(req, timeout=30) as response:
        return json.loads(response.read())

def _safe(value):
    text = json.dumps(value)
    text = re.sub(r"https?://[^\"\s]*jupyter[^\"\s]*", "[redacted-jupyter-link]", text, flags=re.I)
    text = re.sub(r"(Bearer\s+)[^\"\s]+", r"\1[redacted]", text, flags=re.I)
    return json.loads(text)

def _read_registry():
    if not REGISTRY.exists(): return {"instances": [], "events": []}
    return json.loads(REGISTRY.read_text())

def _write_registry(data):
    REGISTRY.parent.mkdir(parents=True, exist_ok=True)
    REGISTRY.write_text(json.dumps(data, indent=2) + "\n")

def status():
    raw = _json("GET", "/instances"); rows = raw.get("data", raw) if isinstance(raw, dict) else raw
    allow = ("id", "name", "status", "ip", "region", "instance_type")
    print(json.dumps([{k: row[k] for k in allow if k in row} for row in rows], indent=2))

def types():
    raw = _json("GET", "/instance-types"); data = raw.get("data", raw)
    allowed = {k: v for k, v in data.items() if k == "gpu_1x_a10"} if isinstance(data, dict) else [x for x in data if x.get("name") == "gpu_1x_a10"]
    print(json.dumps(_safe(allowed), indent=2))

def launch(args):
    if args.instance_type.upper() != "A10": raise SystemExit("Only A10 is allowed")
    raw_types = _json("GET", "/instance-types"); types_data = raw_types.get("data", raw_types)
    spec = types_data.get("gpu_1x_a10") if isinstance(types_data, dict) else next((x for x in types_data if x.get("name")=="gpu_1x_a10"), None)
    price = spec.get('instance_type', {}).get('price_cents_per_hour', 999999) if spec else 999999
    regions = {r['name'] for r in spec.get('regions_with_capacity_available', [])} if spec else set()
    if price > 129 or args.region not in regions:
        raise SystemExit("A10 price/capacity gate refused launch")
    body = {"name": args.name, "region_name": args.region, "instance_type_name": "gpu_1x_a10", "ssh_key_names": ["aiden-mac"], "quantity": 1}
    result = _json("POST", "/instance-operations/launch", body); data = result.get("data", result)
    reg = _read_registry(); reg.setdefault("instances", []).extend(
        [{"id": i, "name": args.name, "region": args.region, "registered_at": time.time()}
         for i in data.get("instance_ids", [])])
    reg.setdefault("events", []).append({"action": "launch", "name": args.name, "region": args.region, "timestamp": time.time()})
    _write_registry(reg); print(json.dumps(_safe(result), indent=2))

def terminate(args):
    reg = _read_registry(); known = {str(x.get("id")) for x in reg.get("instances", [])}
    if args.id not in known: raise SystemExit("Refusing unknown or unregistered instance ID")
    result = _json("POST", "/instance-operations/terminate", {"instance_ids": [args.id]})
    reg.setdefault("events", []).append({"action": "terminate", "id": args.id, "timestamp": time.time()})
    _write_registry(reg); print(json.dumps(_safe(result), indent=2))

def register(args):
    raw = _json("GET", "/instances"); rows = raw.get("data", raw)
    row = next((x for x in rows if x.get("id") == args.id), None)
    if not row or row.get("instance_type", {}).get('name') != "gpu_1x_a10": raise SystemExit("Refusing absent or non-A10 instance")
    reg = _read_registry()
    if args.id not in {r['id'] for r in reg.setdefault('instances', [])}:
        reg['instances'].append({"id": args.id, "name": row.get("name"), "region": row.get("region"), "registered_at": time.time()})
        reg.setdefault("events", []).append({"action":"register","id":args.id,"timestamp":time.time()})
        _write_registry(reg)
    print(json.dumps({"registered":args.id}))

def main(argv=None):
    p = argparse.ArgumentParser(); sub = p.add_subparsers(dest="command", required=True)
    sub.add_parser("status").set_defaults(fn=lambda args: status()); sub.add_parser("types").set_defaults(fn=lambda args: types())
    l = sub.add_parser("launch"); l.add_argument("--name", required=True); l.add_argument("--region", default="us-west-1"); l.add_argument("--instance-type", default="A10"); l.set_defaults(fn=launch)
    t = sub.add_parser("terminate"); t.add_argument("id"); t.set_defaults(fn=terminate)
    r = sub.add_parser("register"); r.add_argument("id"); r.set_defaults(fn=register)
    args = p.parse_args(argv); args.fn(args)
if __name__ == "__main__": main()
