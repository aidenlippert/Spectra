"""Fresh exact cap on the frozen old recipe, compared with joint reoptimization."""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import json
import sys
from spin_coherence_fixed_limit_replay import replay as fixed_replay, VECTORS
from experiments.marginal_pure_coherence import LABELS
from pair_transfer_separation import read_verified

ROOT = Path(__file__).resolve().parents[3]


def source_hashes(files):
    for module in tuple(sys.modules.values()):
        name = getattr(module, '__file__', None)
        if name:
            path = Path(name).resolve()
            if path.is_relative_to(ROOT / 'experiments') or path.is_relative_to(ROOT / 'results/marginal_graded_hubbard8/discovery'):
                files.add(path)
    return {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(files)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('directory', type=Path)
    parser.add_argument('--previous-directory', type=Path, required=True)
    parser.add_argument('--case', choices=('W_zero', 'W_plus_1'), required=True)
    args = parser.parse_args()
    folder, previous = args.directory.resolve(), args.previous_directory.resolve()
    ep, op = folder / 'range_two_replay.json', previous / 'range_two_replay.json'
    new_energy, old_energy = read_verified(ep), read_verified(op)
    match_path = folder / 'previous_family_comparison.json'
    match = read_verified(match_path)
    if (match['source_sha256'].get(str(ep.relative_to(ROOT))) != hashlib.sha256(ep.read_bytes()).hexdigest()
            or match['source_sha256'].get(str(op.relative_to(ROOT))) != hashlib.sha256(op.read_bytes()).hexdigest()):
        raise ValueError('Joint and frozen energies are not bound to the same matched-family comparison')
    cp = previous / 'profile_joint_r1_2_certificate.json'
    c = json.loads(cp.read_text())
    pp = ROOT / 'results/marginal_graded_hubbard8/spin_coherence' / args.case / 'fixed_limit/family_proposal.json'
    proposal = json.loads(pp.read_text())
    original = ROOT / proposal['seed_certificate']
    if cp.read_bytes() != original.read_bytes() or list(LABELS.values()) != list(VECTORS):
        raise ValueError('Frozen recipe or new operators changed')
    result = fixed_replay(c, proposal)
    if F(result['accepted_seed_lower']) != F(old_energy['lower_replay']['periodic_lower_density']):
        raise ValueError('Frozen lower differs from accepted old energy')
    files = {Path(__file__).resolve(), cp, pp, original, op, match_path}
    result['source_sha256'] = source_hashes(files)
    fp = previous / 'fixed_new_coefficient_limit_replay.json'
    if fp.exists():
        raise ValueError('Refusing to overwrite fixed-family proof')
    fp.write_text(json.dumps(result, indent=2) + '\n')
    lower, cap = F(new_energy['lower_replay']['periodic_lower_density']), F(result['fixed_recipe_family_ceiling'])
    output = {'accepted': True, 'joint_periodic_lower': str(lower), 'frozen_old_recipe_ceiling': str(cap),
              'joint_excess_above_frozen_ceiling': str(lower-cap), 'joint_excess_float': float(lower-cap),
              'strict_joint_improvement_beyond_frozen_recipe': lower > cap,
              'source_sha256': source_hashes({Path(__file__).resolve(), ep, op, fp, cp, pp, original, match_path}),
              'scope': 'Comparison of a freshly reconstructed exact ceiling for arbitrary two new coefficients with every older field frozen, against the accepted jointly optimized energy. This is not separation from the reoptimized previous full family.'}
    (folder / 'frozen_limit_comparison.json').write_text(json.dumps(output, indent=2) + '\n')
    print(json.dumps({'accepted': True, 'joint_excess': float(lower-cap)}))


if __name__ == '__main__':
    main()
