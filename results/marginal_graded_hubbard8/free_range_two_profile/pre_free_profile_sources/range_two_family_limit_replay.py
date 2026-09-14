"""Independent exact family-limit replay matched to a range-two energy proof."""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import json
import sys
import time

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from experiments.marginal_range_two_family_limit import replay


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('directory', type=Path)
    folder = parser.parse_args().directory.resolve()
    started = time.monotonic()
    proposal_path = folder/'diagonal_family_limit_proposal.json'
    source = json.loads(proposal_path.read_text())
    if source['kind'] != 'joint_diagonal_range2_family_limit_proposal_v1':
        raise ValueError('Explicit range-two proposal required')
    source['kind'] = 'joint_diagonal_range2_family_limit_v1'
    source['scope'] = 'Exact candidate input for fixed range-two profile. Acceptance requires independent physical replay.'
    certificate_path = folder/'range_two_family_limit_certificate.json'
    certificate_path.write_text(json.dumps(source, indent=2)+'\n')
    result = replay(source)
    energy_path = folder/'profile_joint_r1_2_certificate.json'
    energy_receipt_path = folder/'range_two_replay.json'
    energy = json.loads(energy_path.read_text())
    accepted = json.loads(energy_receipt_path.read_text())
    if not accepted['accepted'] or accepted['source_sha256'][str(energy_path.relative_to(ROOT))] != hashlib.sha256(energy_path.read_bytes()).hexdigest():
        raise ValueError('Energy proof does not match the current certificate')
    if set(energy['target']) != {'U', 't', 'V', 'W'} or any(F(energy['target'][k]) != v for k, v in [('U', F(4)), ('t', F(1)), ('V', F(1, 2)), ('W', F(source['W']))]):
        raise ValueError('Fixed physical target differs')
    if energy['vector'] != source['half_vector'] or energy['joint']['vector'] != source['charged_vector'] or F(energy['joint']['ratio']) != F(source['ratio']):
        raise ValueError('Projector sources or ratio differ')
    if F(energy['projector_sum_ceiling'])/energy['windows'] != F(source['theta_half']) or F(energy['joint']['projector_sum_ceiling'])/energy['joint']['windows'] != F(source['theta_joint']):
        raise ValueError('Fixed projector ceilings differ')
    if list(map(F, energy['local_window']['range_two_density_profile'])) != list(map(F, result['range_two_density_profile'])):
        raise ValueError('Range-two profile differs from the fixed dual profile')
    diagonal = {int(s): F(v) for s, v in energy['telescoping_diagonal'].items()}
    rebuilt = {}
    for shape in source['diagonal_shapes']:
        shape = {int(s): F(v) for s, v in shape.items()}
        if set(shape) & set(rebuilt):
            raise ValueError('Comparison requires disjoint supplied shapes')
        ratios = {diagonal.get(s, F(0))/v for s, v in shape.items() if v}
        if len(ratios) != 1:
            raise ValueError('Energy diagonal is outside the fixed span')
        coefficient = ratios.pop()
        rebuilt.update({s: coefficient*v for s, v in shape.items()})
    if {s: v for s, v in rebuilt.items() if v} != diagonal:
        raise ValueError('Energy has uncovered diagonal entries')
    lower = F(accepted['lower_replay']['periodic_lower_density'])
    upper = F(result['periodic_family_upper'])
    if upper < lower:
        raise ValueError('Fixed-family cap contradicts the accepted energy lower')
    files = {Path(__file__).resolve(), proposal_path, certificate_path, energy_path, energy_receipt_path}
    for module in tuple(sys.modules.values()):
        path = getattr(module, '__file__', None)
        if path and str(Path(path).resolve()).startswith(str(ROOT/'experiments')+'/'):
            files.add(Path(path).resolve())
    receipt = {'accepted': True, 'family_replay': result, 'accepted_periodic_lower': str(lower),
               'periodic_family_upper': str(upper), 'family_gap': str(upper-lower),
               'family_gap_float': float(upper-lower), 'seconds': time.monotonic()-started,
               'source_sha256': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(files)},
               'scope': 'Exact fixed-profile range-two family cap, matched to accepted energy sources/ceilings/ratio/shape span. Not a physical energy upper or a cap over other range-two profiles.'}
    (folder/'range_two_family_limit_replay.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps({key: value for key, value in receipt.items() if key not in ('family_replay', 'source_sha256')}), flush=True)


if __name__ == '__main__':
    main()
