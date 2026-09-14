"""Standard-library lower and physical upper replay for a range-two target."""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import json
import sys
import time

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from experiments.marginal_projector_extendibility import replay
from experiments.marginal_range_two_transfer import compile_block, contract, enclose


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('certificate', type=Path)
    parser.add_argument('--psd-witnesses', type=Path, required=True)
    args = parser.parse_args()
    started = time.monotonic()
    c = json.loads(args.certificate.read_text())
    if c.get('kind') not in ('hubbard_projector_extension_v6','hubbard_projector_extension_v7','hubbard_projector_extension_v8','hubbard_projector_extension_v9','hubbard_projector_extension_v10','hubbard_projector_extension_v11','hubbard_projector_extension_v12','hubbard_projector_extension_v13','hubbard_projector_extension_v14'):
        raise ValueError('Explicit range-two certificate required')
    target = c['target']
    witness_data = json.loads(args.psd_witnesses.read_text())
    lower = replay(c, psd_witnesses=witness_data['witnesses'])
    lower_seconds = time.monotonic()-started
    if lower['local_sum_dimensions'] != 4096 or len(lower['local_sectors']) != 94:
        raise ValueError('Incomplete all-Fock local coverage')
    print('exact range-two lower accepted', flush=True)
    source_path = ROOT/'results/marginal_graded_hubbard8/density_transfer/certificate.json'
    source = json.loads(source_path.read_text())
    # This fixed physical trial state is transferred without optimizing its filter.
    # Its source target selects a recipe, not the Hamiltonian being evaluated.
    choices = [case for case in source['targets'] if case['target'] == {'U':'4','t':'1','V':'1/2'}]
    if len(choices) != 1 or c['chain_sites'] != source['sites'] or c['chain_sites'] % 8:
        raise ValueError('Fixed physical recipe and multiple-of-eight chain required')
    case = choices[0]
    compiled = compile_block(source['upper'], source['hamiltonian'], **{key: F(value) for key, value in target.items()})
    exact = contract(compiled, case['a'], case['b'], 3)
    small = enclose(compiled, case['a'], case['b'], 3, 160)
    if not F(small['energy_lower']) <= F(exact['energy']) <= F(small['energy_upper']):
        raise ValueError('Exact24-site physical energy is outside its enclosure')
    upper = enclose(compiled, case['a'], case['b'], c['chain_sites']//8, 160)
    if upper['target'] != target or upper['sites'] != c['chain_sites']:
        raise ValueError('Physical upper target or size mismatch')
    lo, hi = F(lower['open_lower_density']), F(upper['upper_per_site'])
    if lo > hi:
        raise ValueError('Lower exceeds the physical upper')
    files = {Path(__file__).resolve(), args.certificate.resolve(), args.psd_witnesses.resolve(), source_path}
    for module in tuple(sys.modules.values()):
        path = getattr(module, '__file__', None)
        if path and str(Path(path).resolve()).startswith(str(ROOT/'experiments')+'/'):
            files.add(Path(path).resolve())
    result = {'accepted': True, 'target': target, 'sites': c['chain_sites'],
              'lower_per_site': str(lo), 'upper_per_site': str(hi), 'width_per_site': str(hi-lo),
              'lower_replay': lower, 'upper_replay': upper,
              'physical_source_replay': compiled['source_upper_replay'],
              'exact_24_site': exact, 'enclosed_24_site': small,
              'transferred_filter': {'source_target':case['target'],'a':case['a'],'b':case['b']},
              'lower_seconds': lower_seconds, 'seconds': time.monotonic()-started,
              'source_sha256': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(files)},
              'scope': 'Exact half-filled open U,t,V,W chain interval with next-nearest density coupling. All-Fock lower and a freshly evaluated physical filtered-state upper. Projector sources and the filter recipe are transferred from the W0 experiment; physical profiles were rescaled and the scalar lower threshold recomputed for the new U,V,W target. Correction coefficients are prescribed by the candidate; the separate frozen-transfer comparison checks which coefficients were preserved or removed. No coefficient optimization was performed. The fixed48-dimensional upper recurrence includes six-site dressed patches. Historical norm-perturbation comparison uses N-2 added pairs of norm one. No generic molecular, general representability or requested-accuracy scaling conclusion.'}
    args.certificate.with_name('range_two_replay.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({key: value for key, value in result.items() if key not in ('lower_replay', 'upper_replay', 'physical_source_replay', 'exact_24_site', 'enclosed_24_site', 'source_sha256')}), flush=True)


if __name__ == '__main__':
    main()
