"""Explicit conditional development cases; later cold cases use fresh inputs."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
from research.interacting_scaling_20260915.budget import ROOT, OUT, dump
from research.interacting_scaling_20260915.dictionary import windows


def initialize(name, source, widths, collective=True, complete=False, collective_widths=(), nonsinglet_source=None):
    case = OUT/'cases'/name
    case.mkdir(parents=True, exist_ok=False)
    inputs = {}
    for relative in ('fixture.json', 'mps/state.json', 'upper.json', 'nonsinglet.json'):
        path = nonsinglet_source if relative == 'nonsinglet.json' and nonsinglet_source else source/relative
        if not path.exists():
            raise ValueError(('Required conditional input missing', str(path)))
        target = case/relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(path, target)
        inputs[relative] = {'source': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}
    if (source/'rotation.json').exists():
        shutil.copyfile(source/'rotation.json', case/'rotation.json')
    data = json.loads((case/'fixture.json').read_text())
    clusters = [c for width in widths for c in windows(data['modes']//2, width)]
    design = {'clusters': clusters, 'collective_pairs': collective, 'complete': complete,
        'collective_clusters': [c for width in collective_widths for c in windows(data['modes']//2, width)],
        'max_Gram_entries': 2000000, 'max_coefficient_rows': 180000,
        'global_particle_number_only': True, 'all_declared_Gram_cross_terms_retained': True}
    dump(case/'design.json', design)
    dump(case/'input_dependencies.json', {'conditional_development_case': True,
        'inputs': inputs, 'source_discovery_and_construction_cost_additional': True,
        'old_Gram_or_coefficient_map_input': False, 'old_winning_factor_input': False})
    return case


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('name')
    p.add_argument('source', type=Path)
    p.add_argument('--widths', type=int, nargs='*', default=[2])
    p.add_argument('--no-collective', action='store_true')
    p.add_argument('--complete', action='store_true')
    p.add_argument('--collective-widths', type=int, nargs='*', default=[])
    p.add_argument('--nonsinglet-source', type=Path)
    a = p.parse_args()
    print(initialize(a.name, a.source.resolve(), a.widths, not a.no_collective, a.complete, a.collective_widths, a.nonsinglet_source))
