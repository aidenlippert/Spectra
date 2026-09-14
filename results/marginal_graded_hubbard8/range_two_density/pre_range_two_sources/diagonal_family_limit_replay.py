"""Independent standard-library replay of the proposed nine-shape dual limit."""
from pathlib import Path
from fractions import Fraction as F
import hashlib
import json
import argparse
import sys
import time

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from experiments.marginal_diagonal_family_limit import replay


def main():
    started = time.monotonic()
    folder = ROOT/'results/marginal_graded_hubbard8/joint_projector/signed_density/symmetric_diagonals_9/extra_31'
    parser=argparse.ArgumentParser();parser.add_argument('--directory',type=Path,default=folder)
    folder=parser.parse_args().directory.resolve()
    proposal_path = folder/'diagonal_family_limit_proposal.json'
    source = json.loads(proposal_path.read_text())
    if source['kind'] != 'joint_diagonal_family_limit_proposal_v1':
        raise ValueError('Expected an explicitly untrusted numerical proposal')
    source['kind'] = 'joint_diagonal_family_limit_v1'
    source['scope'] = 'Candidate exact input for the bounded nine-shape family verifier. Acceptance is recorded separately by fresh physical-expectation replay.'
    certificate_path = folder/'diagonal_family_limit_certificate.json'
    certificate_path.write_text(json.dumps(source, indent=2)+'\n')
    result = replay(source)
    energy_path = folder/'profile_joint_r1_2_certificate.json'
    energy_receipt_path = folder/'profile_joint_r1_2_certificate_accelerated_replay.json'
    if not energy_receipt_path.exists():energy_receipt_path=folder/'congruence_witnesses_replay.json'
    energy = json.loads(energy_path.read_text())
    receipt = json.loads(energy_receipt_path.read_text())
    if not receipt['accepted'] or receipt['source_sha256'][str(energy_path.relative_to(ROOT))] != hashlib.sha256(energy_path.read_bytes()).hexdigest():
        raise ValueError('Energy receipt is not bound to the current certificate')
    if energy['vector'] != source['half_vector'] or energy['joint']['vector'] != source['charged_vector']:
        raise ValueError('Fixed projector sources differ')
    if F(energy['joint']['ratio']) != F(source['ratio']):
        raise ValueError('Fixed charged ratio differs')
    if F(energy['projector_sum_ceiling'])/energy['windows'] != F(source['theta_half']) or F(energy['joint']['projector_sum_ceiling'])/energy['joint']['windows'] != F(source['theta_joint']):
        raise ValueError('Fixed ceilings differ')
    # Current shapes have disjoint support. Verify the energy correction is
    # in their span exactly, rather than assuming that from folder names.
    correction = {int(s): F(v) for s, v in energy['telescoping_diagonal'].items()}
    rebuilt = {}
    coefficients = []
    for shape in source['diagonal_shapes']:
        shape = {int(s): F(v) for s, v in shape.items()}
        if set(rebuilt) & set(shape):
            raise ValueError('This comparison driver requires disjoint shape supports')
        values = {correction.get(s, F(0))/v for s, v in shape.items() if v}
        if len(values) != 1:
            raise ValueError('Energy correction lies outside the certified span')
        coefficient = values.pop()
        coefficients.append(str(coefficient))
        rebuilt.update({s: coefficient*v for s, v in shape.items()})
    if {s: v for s, v in rebuilt.items() if v} != correction:
        raise ValueError('Energy correction has uncovered entries')
    lower = F(receipt['lower_replay']['periodic_lower_density'])
    upper = F(result['periodic_family_upper'])
    if upper < lower:
        raise ValueError('Family upper contradicts the accepted lower certificate')
    files = {Path(__file__).resolve(), certificate_path, proposal_path, energy_path, energy_receipt_path}
    for module in tuple(sys.modules.values()):
        path = getattr(module, '__file__', None)
        if path and str(Path(path).resolve()).startswith(str(ROOT/'experiments')+'/'):
            files.add(Path(path).resolve())
    output = {'accepted': True, 'family_replay': result, 'accepted_periodic_lower': str(lower),
              'periodic_family_upper': str(upper), 'remaining_family_gap': str(upper-lower),
              'remaining_family_gap_float': float(upper-lower), 'energy_shape_coefficients': coefficients,
              'seconds': time.monotonic()-started,
              'source_sha256': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(files)},
              'scope': 'Exact family-limit replay from physical determinant actions and all nine diagonal expectations, with explicit source/ratio/ceiling/span match to the accepted energy certificate. This caps attainable lower certificates and is not a physical energy upper bound.'}
    (folder/'diagonal_family_limit_replay.json').write_text(json.dumps(output, indent=2)+'\n')
    print(json.dumps({k: v for k, v in output.items() if k not in ('family_replay', 'source_sha256')}), flush=True)


if __name__ == '__main__':
    main()
