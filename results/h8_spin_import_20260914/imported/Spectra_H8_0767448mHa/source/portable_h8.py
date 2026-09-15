"""Portable adapter for the pinned H8 inputs; original files stay unchanged."""
import argparse
import hashlib
import json
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parent
ORIGINAL_ROOT = "/Users/aidenlippert/Documents/Spectra/"
PREPARED = ROOT / "results/sector_quotient_20260914/prepared_linear_closure"


def verify_manifest(path):
    manifest = json.loads(path.read_text())
    for relative, record in manifest["files"].items():
        target = ROOT / relative
        digest = hashlib.sha256()
        with target.open("rb") as stream:
            for chunk in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(chunk)
        if target.stat().st_size != record["size"] or digest.hexdigest() != record["sha256"]:
            raise ValueError("Changed or incomplete file: " + relative)
    return len(manifest["files"])


def portable_metadata(meta):
    meta = dict(meta)
    for key in ("fixture", "seed", "compact_frame", "original_seed"):
        value = meta.get(key)
        if isinstance(value, str) and value.startswith(ORIGINAL_ROOT):
            meta[key] = str(ROOT / value[len(ORIGINAL_ROOT):])
    return meta


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("action", choices=("check", "solve"))
    parser.add_argument("--seconds", type=int, default=30)
    parser.add_argument("--tag", default="pro_trial")
    args = parser.parse_args()
    if not 1 <= args.seconds <= 600:
        parser.error("Use a positive solve budget of at most 600 seconds.")
    if not args.tag or any(c not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-" for c in args.tag):
        parser.error("Tag must contain only letters, digits, underscores, or hyphens.")
    os.chdir(ROOT)
    count = verify_manifest(ROOT / "HANDOFF_MANIFEST.json")
    cache_manifest = ROOT / "CACHE_MANIFEST.json"
    if cache_manifest.exists():
        count += verify_manifest(cache_manifest)
    import numpy as np
    from scipy import sparse
    from research.sector_quotient_20260914.search import Operator
    metadata = portable_metadata(json.loads((PREPARED / "frame.json").read_text()))
    with np.load(PREPARED / "seed.npz", allow_pickle=False) as seed:
        for key in seed.files:
            if seed[key].dtype.hasobject:
                raise ValueError("Unexpected object array")
        arrays = len(seed.files)
    for path in PREPARED.glob("*.npz"):
        if path.name != "seed.npz":
            sparse.load_npz(path)
    for path in PREPARED.glob("*.npy"):
        np.load(path, mmap_mode="r", allow_pickle=False)
    operator = Operator(PREPARED, metadata)
    print(json.dumps({"verified_files": count, "seed_arrays": arrays,
                      "coefficient_rows": len(operator.rhs),
                      "gram_entries": sum(q.size for q in operator.Q),
                      "cache_present": cache_manifest.exists(),
                      "energy_replayed": False}, indent=2), flush=True)
    if args.action == "check":
        return
    if not cache_manifest.exists():
        raise ValueError("Extract Spectra-H8-solver-cache.zip before the conditioned solve.")
    from research.sector_quotient_20260914 import eliminated
    class PortableOperator(Operator):
        def __init__(self, prepared, meta):
            super().__init__(prepared, portable_metadata(meta))
    eliminated.Operator = PortableOperator
    eliminated.OUT = ROOT / "pro_runs"
    eliminated.run(
        args.tag, seconds=args.seconds, max_cg=20, mu=.005, fixed_mu=True,
        prepared_path=str(PREPARED),
        normal_cache=str(ROOT / "results/sector_quotient_20260914/normal_linear_closure"),
        restart=str(ROOT / "results/sector_quotient_20260914/candidates/maximize_fixed_mu/checkpoint.npz"),
    )
    print("New outputs are untrusted proposals and require independent exact replay.")


if __name__ == "__main__":
    main()
