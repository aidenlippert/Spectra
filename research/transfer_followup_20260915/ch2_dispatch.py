"""Use lowercase job identifiers; retain the original pre-execution rejections."""
import hashlib
import subprocess
import time
from research.transfer_solver_20260915.budget import ROOT, OUT, dump
from research.transfer_solver_20260915.campaign import STD, CHEM, PREFIX


def run():
    folder = OUT/'ch2_model_study'
    source = ROOT/'research/ch2_model_study_20260915/study.py'
    dump(folder/'dispatch_repair.json', {'original_rejected_dispatches': 4,
        'reason': 'Uppercase S in job names violated the existing runner name rule before the chemistry processes started',
        'change': 'Lowercase job identifiers only', 'physical_method_source_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
        'original_dispatch_overhead': 'Not separately timed; no chemistry worker started in those attempts'})
    start = time.monotonic()
    results = []
    for basis in ('cc-pvdz', 'cc-pvtz'):
        for spin in (0, 1):
            code = subprocess.run([STD, '-B', '-S', '-m', PREFIX+'budget', '--name', f'ch2_{basis}_s{spin}',
                '--seconds', '360', '--', CHEM, '-B', '-m', 'research.ch2_model_study_20260915.study',
                'case', '--basis', basis, '--spin', str(spin)], cwd=ROOT).returncode
            results.append({'basis': basis, 'spin': spin, 'exit_code': code})
    dump(folder/'corrected_dispatch_execution.json', {'seconds': time.monotonic()-start, 'cases': results})


if __name__ == '__main__':
    run()
