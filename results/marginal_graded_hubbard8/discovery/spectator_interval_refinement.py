"""Audit nested exact bounds on the same fixed spectator certificate family."""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import json

ROOT = Path(__file__).resolve().parents[3]
FIXED = ('kind', 'half_vector', 'charged_vector', 'ratio', 'theta_half',
         'theta_joint', 'diagonal_shapes', 'W', 'range_two_density_profile')


def compare(old, new, old_certificate, new_certificate):
    for receipt in (old, new):
        if not receipt['accepted'] or not receipt.get('spectator_hopping'):
            raise ValueError('Accepted spectator family required')
    for key in FIXED:
        if old_certificate[key] != new_certificate[key]:
            raise ValueError('Fixed family differs: ' + key)
    lo0, lo1 = (F(r['accepted_periodic_lower']) for r in (old, new))
    hi0, hi1 = (F(r['periodic_family_upper']) for r in (old, new))
    for r, lo, hi in ((old, lo0, hi0), (new, lo1, hi1)):
        if hi < lo or F(r['family_gap']) != hi - lo:
            raise ValueError('Invalid exact interval')
    if lo1 < lo0 or hi1 > hi0 or hi1 - lo1 >= hi0 - lo0:
        raise ValueError('Strict nested interval refinement required')
    return {
        'accepted': True,
        'old_periodic_lower': str(lo0), 'new_periodic_lower': str(lo1),
        'old_ceiling': str(hi0), 'new_ceiling': str(hi1),
        'lower_improvement': str(lo1 - lo0),
        'ceiling_improvement': str(hi0 - hi1),
        'old_gap': str(hi0 - lo0), 'new_gap': str(hi1 - lo1),
        'new_gap_float': float(hi1 - lo1),
        'gap_reduction_factor': str((hi0 - lo0) / (hi1 - lo1)) if hi1 > lo1 else None,
        'scope': 'Strictly nested exact bounds on the optimal periodic lower certificate in the identical fixed family. Family ceilings are not physical ground-energy upper bounds. No exact optimum attainment or general representability claim.',
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('case_directory', type=Path)
    folder = parser.parse_args().case_directory.resolve()
    old_dir, new_dir = folder / 'completed_basis', folder / 'thermal/final'
    files = {Path(__file__).resolve()}
    values = []
    for directory in (old_dir, new_dir):
        rp = directory / 'range_two_family_limit_replay.json'
        cp = directory / 'range_two_family_limit_certificate.json'
        receipt = json.loads(rp.read_text())
        for name, expected in receipt['source_sha256'].items():
            if hashlib.sha256((ROOT / name).read_bytes()).hexdigest() != expected:
                raise ValueError('Stale family source: ' + name)
        values.append((receipt, json.loads(cp.read_text())))
        files.update((rp, cp))
    result = compare(values[0][0], values[1][0], values[0][1], values[1][1])
    result['source_sha256'] = {
        str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sorted(files)
    }
    (new_dir / 'interval_refinement.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'accepted': True, 'gap': result['new_gap_float'],
                      'reduction_factor': float(F(result['gap_reduction_factor']))}))


if __name__ == '__main__':
    main()
