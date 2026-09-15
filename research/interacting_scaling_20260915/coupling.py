"""Exact support partition of a molecular Hamiltonian; only global N is fixed."""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import shutil
from experiments.marginal_symbolic import add, decode, encode, hermitian, scale
from research.molecular_collective_20260913.core import digest


def coupled_fixture(original, partition, coupling):
    from research.nvidia_followup_20260915.strict_replay import require_supported_sector
    require_supported_sector(original)
    m, n = original['modes'], original['particles']
    flat = [i for fragment in partition for i in fragment]
    if not partition or any(not fragment for fragment in partition) or sorted(flat) != list(range(m//2)):
        raise ValueError('Fragments must form a disjoint complete spatial-orbital partition')
    if any(type(i) is not int for i in flat) or not isinstance(coupling, F) or not 0 <= coupling <= 1:
        raise ValueError('Exact coupling in [0,1] and integer orbital labels required')
    h = decode(original['hamiltonian'], m, 4)
    if not hermitian(h) or any(sum(2*c-1 for c, i in w) for w in h):
        raise ValueError('Hermitian number-conserving molecular Hamiltonian required')
    fragment_of = {i: k for k, fragment in enumerate(partition) for i in fragment}
    inside, between = {}, {}
    for word, value in h.items():
        touched = {fragment_of[i//2] for c, i in word}
        (inside if len(touched) <= 1 else between)[word] = value
    total = add(inside, scale(between, coupling))
    if not hermitian(inside) or not hermitian(between) or add(inside, between) != h:
        raise ValueError('Support split failed exact reconstruction')
    data = {'kind': 'exact_molecular_coupling_diagnostic_v1', 'modes': m, 'particles': n,
        'hamiltonian': encode(total), 'original_fixture_sha256': digest(original),
        'coupling': str(coupling), 'spatial_partition': partition,
        'global_particle_number_only': True,
        'physical_original_endpoint': coupling == 1}
    record = {'within_fragment_terms': len(inside), 'between_fragment_terms': len(between),
        'between_coefficient_l1_Ha': str(sum(abs(c) for c in between.values())),
        'between_terms_changing_fragment_charge': sum(any(sum((2*c-1) for c, i in w if i//2 in fragment) for fragment in partition) for w in between),
        'lambda_one_reconstructs_original_exactly': add(inside, between) == h,
        'local_particle_number_constraints_added': False}
    return data, record


def initialize(source, name, coupling, fragment_width=2, widths=(2, 3), previous=None):
    from research.interacting_scaling_20260915.budget import OUT, dump
    from research.interacting_scaling_20260915.dictionary import windows
    data = json.loads((source/'fixture.json').read_text())
    s = data['modes']//2
    if type(fragment_width) is not int or not 1 <= fragment_width <= s:
        raise ValueError('Invalid fragment width')
    partition = [list(range(i, min(s, i+fragment_width))) for i in range(0, s, fragment_width)]
    changed, receipt = coupled_fixture(data, partition, coupling)
    # Retain the literal original fixture at the physical endpoint. Diagnostic
    # metadata lives in its own receipt and cannot change the bound's identity.
    fixture = data if coupling == 1 else changed
    case = OUT/'cases'/name
    case.mkdir(parents=True, exist_ok=False)
    dump(case/'fixture.json', fixture)
    dump(case/'coupling.json', {**receipt, 'coupling': str(coupling), 'partition': partition,
        'source': str(source), 'original_fixture_sha256': digest(data),
        'physical_endpoint': coupling == 1, 'continuation_from': str(previous) if previous else None})
    dump(case/'design.json', {'clusters': [c for width in widths for c in windows(s, width)],
        'collective_pairs': True, 'complete': False, 'max_Gram_entries': 2000000,
        'max_coefficient_rows': 180000, 'global_particle_number_only': True})
    if coupling == 1 and (source/'rotation.json').exists():
        shutil.copyfile(source/'rotation.json', case/'rotation.json')
    if previous:
        old_data = json.loads((previous/'fixture.json').read_text())
        state = json.loads((previous/'mps/state.json').read_text())
        if state['fixture_sha256'] != digest(old_data) or any(old_data[k] != fixture[k] for k in ('modes', 'particles')):
            raise ValueError('Continuation requires the same verified orbital and particle domain')
        prior = json.loads((previous/'coupling.json').read_text())
        if prior['original_fixture_sha256'] != digest(data) or prior['partition'] != partition:
            raise ValueError('Continuation cannot cross an undeclared orbital partition or original input')
        state['fixture_sha256'] = digest(fixture)
        dump(case/'continuation_seed.json', state)
    return case


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('source', type=Path); p.add_argument('name'); p.add_argument('coupling', type=F)
    p.add_argument('--fragment-width', type=int, default=2)
    p.add_argument('--widths', type=int, nargs='+', default=[2, 3])
    p.add_argument('--previous', type=Path)
    a = p.parse_args()
    print(initialize(a.source.resolve(), a.name, a.coupling, a.fragment_width, a.widths, a.previous))
