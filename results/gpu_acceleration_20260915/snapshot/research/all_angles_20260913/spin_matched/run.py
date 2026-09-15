"""Same-host matched X ablation on the original active-space Hamiltonian."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from experiments.marginal_symbolic import decode
import research.certificate_scaling.spin_invariant_discovery as backend


def run(fixture_path, out, variant, seconds):
    if variant not in ('restricted_Sz_X', 'full_number_X'):
        raise ValueError('Unknown multiplier family')
    raw = Path(fixture_path).read_bytes()
    fixture = json.loads(raw)
    original_weight = backend.spin_weight2
    try:
        if variant == 'full_number_X':
            backend.spin_weight2 = lambda word: 0
        result = backend.run(decode(fixture['hamiltonian'], fixture['modes'], 4),
                             fixture['modes'], fixture['particles'], out,
                             solver='SCS', seconds=seconds)
    finally:
        backend.spin_weight2 = original_weight
    receipt = {'variant': variant, 'fixture_sha256': hashlib.sha256(raw).hexdigest(),
               'fixture_path': str(fixture_path), 'original_H': True,
               'seconds': seconds,
               'comparison': 'Same SOS word family and scaling RULE; added X columns change numerical row scaling.',
               'source_receipt': 'receipt.json', 'exact_lower': result['exact']['lower']}
    (Path(out)/'comparison.json').write_text(json.dumps(receipt, indent=2)+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fixture', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--variant', choices=['restricted_Sz_X', 'full_number_X'], required=True)
    parser.add_argument('--seconds', type=float, default=120)
    args = parser.parse_args()
    run(args.fixture, args.out, args.variant, args.seconds)
