"""Bounded numerical five-window overlap probe for the six-site witness."""
from pathlib import Path
import hashlib, json, math, time
import numpy as np
from scipy.linalg import eigvalsh

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / 'results/marginal_graded_hubbard8'
SOURCE = BASE / 'six_site_projector/certificate.json'


def columns(vector, windows=5, block_sites=6):
    amps = {int(s): int(a) for s, a in vector.items() if int(a)}
    norm = sum(a*a for a in amps.values())
    groups = {}
    for offset in range(windows):
        for env in range(4 ** (windows - 1)):
            left = env & ((1 << (2 * offset)) - 1)
            right = env >> (2 * offset)
            col = {left | (s << (2 * offset)) |
                   (right << (2 * (offset + block_sites))): a
                   for s, a in amps.items()}
            sample = next(iter(col))
            support_sites = windows + block_sites - 1
            sector = (sum((sample >> (2*i)) & 1 for i in range(support_sites)),
                      sum((sample >> (2*i+1)) & 1 for i in range(support_sites)))
            groups.setdefault(sector, []).append(col)
    return norm, groups


def main():
    start = time.monotonic()
    cert = json.loads(SOURCE.read_text())
    norm, groups = columns(cert['vector'])
    maxima = []
    dimensions = []
    for sector, cols in sorted(groups.items()):
        gram = np.array([[sum(a * right.get(s, 0) for s, a in left.items())
                          for right in cols] for left in cols], dtype=np.float64)
        value = float(eigvalsh(gram / norm, subset_by_index=[len(cols)-1, len(cols)-1])[0])
        maxima.append(value)
        dimensions.append({'sector': list(sector), 'dimension': len(cols), 'maximum': value})
    maximum = max(maxima)
    ceiling = math.ceil((maximum + 1e-8) * 10**9) / 10**9
    out = {'kind': 'five_window_overlap_probe_v1', 'windows': 5,
           'block_sites': 6, 'column_count': sum(x['dimension'] for x in dimensions),
           'maximum_sector_dimension': max(x['dimension'] for x in dimensions),
           'numerical_sum_maximum': maximum, 'rational_ceiling': f'{ceiling:.9f}',
           'average_fidelity_ceiling': ceiling / 5,
           'accepted_m4_average_ceiling': 0.5424350565,
           'improvement_over_m4': 0.5424350565 - ceiling / 5,
           'seconds': time.monotonic() - start,
           'source_sha256': hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
           'sectors': dimensions,
           'scope': 'Numerical-only low-rank Gram probe; no exact PSD certificate or physical conclusion.'}
    target = BASE / 'six_site_projector/five_window_probe.json'
    target.write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps(out))


if __name__ == '__main__':
    main()
