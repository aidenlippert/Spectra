"""Fresh exact overlapping-projector lower and physical transfer upper proof."""
from pathlib import Path
from fractions import Fraction as F
import sys, json, time, hashlib, argparse

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from experiments.marginal_projector_extendibility import replay, projector_bound
from experiments.marginal_boundary_transfer import compile_block, contract, enclose
import density_window_ceiling

BASE = ROOT / 'results/marginal_graded_hubbard8'
OUT = BASE / 'six_site_projector'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--certificate', default='certificate.json')
    args = parser.parse_args()
    start = time.monotonic()
    inputs = [OUT/args.certificate, BASE/'density_transfer/certificate.json',
              BASE/'density_transfer/window_ceiling_certificate.json']
    source_paths = {Path(__file__).resolve(), Path(density_window_ceiling.__file__).resolve()}
    source_paths.update(Path(m.__file__).resolve() for name, m in tuple(sys.modules.items())
                        if name.startswith('experiments.') and getattr(m, '__file__', None))
    source_paths.update(inputs)
    source_paths.add(BASE/'discovery/singlet_moment_energy.py')
    hashes = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
              for p in sorted(source_paths)}
    c = json.loads(inputs[0].read_text())
    upper_input = json.loads(inputs[1].read_text())
    case = upper_input['targets'][0]
    if c['target'] != case['target'] or c['chain_sites'] != upper_input['sites']:
        raise ValueError('Lower and upper targets must match')
    density_window_ceiling.main()
    ceiling = json.loads((BASE/'density_transfer/window_ceiling_replay.json').read_text())
    if c['vector'] != json.loads(inputs[2].read_text())['physical_state']:
        raise ValueError('Penalty must use the actual obstruction witness')
    lower = replay(c)
    overlaps = [projector_bound(c['vector'], m, B, 6) for m, B in
                [(2, '1375515741/1000000000'), (3, '933300153/500000000'),
                 (4, '1084870113/500000000')] if m != c['windows']]
    overlaps.append(lower['overlap'])
    target = {k: F(v) for k, v in c['target'].items()}
    compiled = compile_block(upper_input['upper'], upper_input['hamiltonian'], **target)
    exact = contract(compiled, case['a'], case['b'], 3)
    small = enclose(compiled, case['a'], case['b'], 3, 160)
    if not F(small['energy_lower']) <= F(exact['energy']) <= F(small['energy_upper']):
        raise ValueError('Exact small-chain energy outside rounded enclosure')
    upper = enclose(compiled, case['a'], case['b'], c['chain_sites']//8, 160)
    lo, hi = F(lower['open_lower_density']), F(upper['upper_per_site'])
    family = F(ceiling['profile_density_ceiling'])
    if not family < lo <= hi:
        raise ValueError('New interval must be nonempty and cross the profile ceiling')
    if F(upper['width_per_site']) >= F(1, 10**30):
        raise ValueError('Trial energy enclosure too wide')
    old = F(case['lower_certificate']['lower'])/5 - F(5, 2*c['chain_sites'])
    for p, digest in hashes.items():
        if hashlib.sha256((ROOT/p).read_bytes()).hexdigest() != digest:
            raise ValueError('Source or input changed during replay: '+p)
    result = {'accepted': True, 'target': c['target'], 'sites': c['chain_sites'],
              'lower_per_site': str(lo), 'upper_per_site': str(hi),
              'width_per_site': str(hi-lo), 'previous_lower_per_site': str(old),
              'lower_improvement_per_site': str(lo-old),
              'fraction_previous_interval_closed': str((lo-old)/(hi-old)),
              'profile_density_ceiling': str(family),
              'improvement_above_profile_ceiling': str(lo-family),
              'lower_replay': lower, 'overlap_replays': overlaps,
              'upper_replay': upper, 'exact_24_site': exact, 'enclosed_24_site': small,
              'source_sha256': hashes, 'seconds': time.monotonic()-start,
              'scope': 'Specified half-filled open U4,t1,V1/2 chain. Full-Fock six-site positivity and all-state overlapping-projector inequality give the lower bound. Independently reconstructed physical filtered-block state gives the upper. General representability, accuracy-cost convergence, and arbitrary Hamiltonian transfer remain open.'}
    output_name = inputs[0].stem.replace('certificate', 'independent_replay')+'.json'
    (OUT/output_name).write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: result[k] for k in ('accepted', 'lower_per_site', 'upper_per_site',
                                           'width_per_site', 'seconds')}), flush=True)


if __name__ == '__main__':
    main()
