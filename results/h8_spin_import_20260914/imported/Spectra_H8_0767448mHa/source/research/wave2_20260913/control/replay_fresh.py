"""Independent standard-library replay; may prepare while a new proof exports."""
import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from experiments.marginal_symbolic import decode
from research.certificate_scaling.wedge_spectral_bound import extract, replay_residual
from research.certificate_scaling.streaming_reference_upper import upper


def main():
    p = argparse.ArgumentParser()
    for name in ('certificate', 'fixture', 'proof', 'reference', 'out'):
        p.add_argument('--'+name, type=Path, required=True)
    args = p.parse_args()
    start = time.monotonic()
    raw = args.certificate.read_bytes(); cert = json.loads(raw)
    fraw = args.fixture.read_bytes(); fixture = json.loads(fraw)
    if sha256(fraw).hexdigest() != 'ec437b52d210f08233a8cc853d0f905dae52bd1f2a923285a25c76ac99385f92':
        raise ValueError('Original-H10 fixture mismatch')
    if any(cert[k] != fixture[k] for k in ('modes', 'particles')):
        raise ValueError('Sector mismatch')
    if decode(cert['hamiltonian'], cert['modes'], 4) != decode(fixture['hamiltonian'], fixture['modes'], 4):
        raise ValueError('Hamiltonian mismatch')
    residual, original = extract(cert)
    extracted = time.monotonic()
    print(json.dumps({'stage':'exact_SOS_and_residual_replayed','seconds':extracted-start}), flush=True)
    while not args.proof.exists():
        if time.monotonic()-extracted > 1200:
            raise TimeoutError('New residual proof did not arrive')
        time.sleep(5)
    praw = args.proof.read_bytes(); proof = json.loads(praw)
    if proof['certificate_sha256'] != sha256(raw).hexdigest():
        raise ValueError('Residual proof belongs to another certificate')
    waited = time.monotonic()-extracted
    lower_record = replay_residual(cert, residual, proof)
    rraw = args.reference.read_bytes(); reference = json.loads(rraw)
    hi, upper_record = upper(fixture, reference['independent_upper'])
    if hi != Fraction(reference['upper']):
        raise ValueError('Upper receipt mismatch')
    lo = Fraction(lower_record['lower']); width = hi-lo
    if width < 0:
        raise ValueError('Inconsistent certified endpoints')
    result = {'method':'fresh_spin_SOS_spectral_plus_fresh_PT2_upper',
              'lower':str(lo),'upper':str(hi),'width':str(width),
              'lower_float':float(lo),'upper_float':float(hi),'width_float':float(width),
              'passes_0_0016_Ha':width <= Fraction(16,10000),
              'fixture_sha256':sha256(fraw).hexdigest(),
              'certificate_sha256':sha256(raw).hexdigest(),
              'proof_sha256':sha256(praw).hexdigest(),
              'reference_sha256':sha256(rraw).hexdigest(),
              'lower_replay':lower_record,'upper_replay':upper_record,
              'SOS_extract_seconds':extracted-start,'proof_wait_seconds':waited,
              'active_replay_seconds':time.monotonic()-start-waited,
              'wall_seconds':time.monotonic()-start,
              'scope':'Independent python -S exact original-H10 two-sided replay. Both lower and upper newly discovered this wave. No FCI coefficients or old lower/proof reused. No scalability or physical-model error claim.'}
    args.out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('lower_replay','upper_replay')}),flush=True)


if __name__ == '__main__':
    main()
