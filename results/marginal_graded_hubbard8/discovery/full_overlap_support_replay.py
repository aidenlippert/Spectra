"""Audit new coherence support against every prior nondiagonal telescope."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import importlib
import json
import sys
from full_overlap_telescope import ROOT, build
from experiments.marginal_local_hubbard_block import _sector
from experiments.marginal_charged_projectors import charged_vectors

BASE = ROOT/'results/marginal_graded_hubbard8/full_overlap'


def main():
    files = {Path(__file__).resolve()}; cases = {}; entries = []; matrices = []
    for case in ('W_zero', 'W_plus_1'):
        path = BASE/case/'telescope_replay.json'; receipt = json.loads(path.read_text()); files.add(path)
        if receipt.get('accepted') is not True:
            raise ValueError('Accepted telescope required')
        for name, expected in receipt['source_sha256'].items():
            p = ROOT/name; files.add(p)
            if hashlib.sha256(p.read_bytes()).hexdigest() != expected:
                raise ValueError('Stale telescope source')
        vector = {int(s): a for s, a in receipt['five_site_vector'].items()}
        denominator, _, matrix = build(vector)
        support = receipt['new_support_entry']; row, column = support['row'], support['column']
        if F(matrix.get((row, column), 0), denominator) != F(support['matrix_element']):
            raise ValueError('New support entry changed')
        matrices.append({key: F(value, denominator) for key, value in matrix.items()})
        entries.append((row, column))
        cp = ROOT/f'results/marginal_graded_hubbard8/three_spectator/{case}/final/range_two_family_limit_certificate.json'
        family = json.loads(cp.read_text()); files.add(cp)
        if {_sector(int(s), 6) for s, a in family['half_vector'].items() if a} != {(3, 3)}:
            raise ValueError('Old half projector support changed')
        charged, _ = charged_vectors(family['charged_vector'])
        if any(any(s.bit_count() not in (5, 7) for s in v) for v in charged):
            raise ValueError('Old charged projector support changed')
        cases[case] = {'new_support_entry': support, 'old_projectors_vanish': True}
    checked = []; counts = {}
    for name, allowed in (('hopping_telescope', {2}), ('spin_telescope', {0, 4}),
                          ('spectator_hopping', {2}), ('pair_transfer', {4}),
                          ('two_spectator_hopping', {2}), ('three_spectator_hopping', {2})):
        module = importlib.import_module('experiments.marginal_'+name)
        labels = module.LABELS if hasattr(module, 'LABELS') else ('single',)
        counts[name] = len(labels)
        for label in labels:
            action = module.actions({label: 1}) if hasattr(module, 'LABELS') else module.actions()
            if len(action) != 4096 or any((r ^ c).bit_count() not in allowed
                                         for c, image in enumerate(action) for r in image):
                raise ValueError('Old telescope has unexpected Fock support')
            if any(action[c].get(r, 0) for r, c in entries):
                raise ValueError('New entry belongs to an old telescope')
            checked.append(name+':'+label)
    minor = [[matrix.get(entry, F(0)) for matrix in matrices] for entry in entries]
    determinant = minor[0][0]*minor[1][1]-minor[0][1]*minor[1][0]
    if not determinant:
        raise ValueError('Two witnesses do not add independent directions modulo the prior support')
    for module in tuple(sys.modules.values()):
        name = getattr(module, '__file__', None)
        if name and str(Path(name).resolve()).startswith(str(ROOT/'experiments')+'/'):
            files.add(Path(name).resolve())
    files.update(Path(sys.modules[name].__file__).resolve()
                 for name in ('full_overlap_telescope', 'full_overlap_density'))
    result = {'accepted': True, 'cases': cases, 'prior_telescope_counts': counts,
              'prior_telescope_count': len(checked), 'all_4096_state_supports_checked': True,
              'new_direction_minor': [[str(v) for v in row] for row in minor],
              'new_direction_minor_determinant': str(determinant),
              'independent_new_directions_modulo_prior_family': 2,
              'source_sha256': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(files)},
              'scope': 'Two independent matrix directions beyond ENERGYv16/FAMILYv12. Exact exhaustive support checks for all71 prior nondiagonal telescopes, plus actual fixed-projector sector checks. Diagonal operators and physical one-body hopping cannot contribute to these eight-bit offdiagonal entries. Independence does not imply an improved energy bound.'}
    output = BASE/'support_replay.json'
    if output.exists():
        raise ValueError('Refusing to overwrite support receipt')
    output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: result[k] for k in ('accepted', 'prior_telescope_count',
                                           'new_direction_minor_determinant', 'independent_new_directions_modulo_prior_family')}))


if __name__ == '__main__':
    main()
