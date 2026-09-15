"""Fresh original-H10 spin-adapted SOS discovery and new residual proof."""
from pathlib import Path
import hashlib
import json
import resource
import sys
import time

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from experiments.marginal_symbolic import decode
from research.certificate_scaling.spin_invariant_discovery import run as discover
from research.certificate_scaling.wedge_spectral_bound import run as spectral


def main():
    started = time.monotonic()
    fixture = ROOT/'results/certificate_scaling/active_space_ladder/h10/fixture.json'
    raw = fixture.read_bytes()
    expected = 'ec437b52d210f08233a8cc853d0f905dae52bd1f2a923285a25c76ac99385f92'
    if hashlib.sha256(raw).hexdigest() != expected:
        raise ValueError('Original-H10 fixture hash mismatch')
    f = json.loads(raw)
    out = ROOT/'results/wave2_20260913/control/h10_spin_scs600'
    result = discover(decode(f['hamiltonian'], f['modes'], 4), f['modes'], f['particles'], out,
                      solver='SCS', seconds=600)
    spectral_out = out.parent/'h10_spin_scs600_spectral'
    proof = spectral(out/'certificate.json', spectral_out)
    record = {'fixture': str(fixture.relative_to(ROOT)), 'fixture_sha256': expected,
              'lower_certificate': str((out/'certificate.json').relative_to(ROOT)),
              'spectral_witness': str((spectral_out/'witness.json').relative_to(ROOT)),
              'source_factors_used': False, 'source_upper_used': False,
              'lower': proof['lower'], 'lower_float': proof['lower_float'],
              'discovery_seconds': result['wall_seconds'],
              'new_residual_proof_seconds': proof['wall_seconds'],
              'whole_process_seconds': time.monotonic()-started,
              'peak_rss_kib_linux': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
              'scope': 'Fresh original-H10 baseline using existing spin-invariant mixed-cubic SOS. No old certificate/proof/state input; no sparse or scalable discovery claim.'}
    (out.parent/'fresh_h10_lower.json').write_text(json.dumps(record, indent=2)+'\n')
    print(json.dumps(record), flush=True)


if __name__ == '__main__':
    main()
