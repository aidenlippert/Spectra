"""Independent replay of the small molecular gap certificate."""
import json
import sys
import time

from research.positive_response_20260913.molecular_diagnostic import ROOT, OUT, sector_gap_check


def run():
    start = time.monotonic(); src = ROOT/'results/molecular_collective_20260913/campaign/h6'
    fixture_bytes = (src/'fixture.json').read_bytes(); tail_bytes = (src/'rank_10/tail.json').read_bytes()
    cert_bytes = (OUT/'sector_gap_certificate.json').read_bytes()
    result = sector_gap_check(json.loads(fixture_bytes), json.loads(tail_bytes), json.loads(cert_bytes))
    forbidden = [name for name in ('numpy', 'scipy', 'cvxpy', 'pyscf') if name in sys.modules]
    if forbidden:
        raise AssertionError('Numerical imports during exact sector-gap replay')
    result.update(numerical_packages_loaded=forbidden, certificate_bytes=len(cert_bytes),
        inherited_fixture_and_tail_bytes=len(fixture_bytes)+len(tail_bytes),
        complete_data_load_and_check_seconds=time.monotonic()-start)
    (OUT/'sector_gap_replay.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'gap_Ha': result['accepted_delta_float_Ha'], 'certificate_bytes': len(cert_bytes),
        'data_load_and_check_seconds': result['complete_data_load_and_check_seconds']}), flush=True)


if __name__ == '__main__':
    run()
