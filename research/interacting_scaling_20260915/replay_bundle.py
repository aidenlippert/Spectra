"""Replay an accepted local-orbital bundle without writing into the sealed case."""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import shutil
from research.interacting_scaling_20260915.complete import run as complete
from research.interacting_scaling_20260915.budget import dump


def run(original, source, output):
    output.mkdir(parents=True, exist_ok=False)
    working = output/'input'
    working.mkdir()
    for name in ('fixture.json', 'upper.json', 'nonsinglet.json', 'rotation.json'):
        shutil.copyfile(source/name, working/name)
    for name in ('mps/state.json', 'solve/export/certificate.json'):
        destination = working/name
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source/name, destination)
    complete(original, working, 'solve', output/'exact')
    recomputed = json.loads((working/'original_interval.json').read_text())
    claimed = json.loads((source/'original_interval.json').read_text())
    if any(F(claimed[k]) != F(recomputed[k]) for k in ('lower_Ha', 'upper_Ha')):
        raise ValueError('Accepted endpoints did not reproduce from the actual bundle')
    dump(output/'comparison.json', {'exact_endpoints_reproduced': True,
        'width_mHa': recomputed['width_mHa'], 'target_met': recomputed['target_met'],
        'source_case_unchanged': str(source), 'original_model': str(original)})


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('original', type=Path)
    p.add_argument('case', type=Path); p.add_argument('new_output', type=Path)
    a = p.parse_args(); run(a.original.resolve(), a.case.resolve(), a.new_output.resolve())
