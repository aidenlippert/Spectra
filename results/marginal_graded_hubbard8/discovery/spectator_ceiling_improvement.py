"""Compare two independently accepted ceilings for the identical fixed family."""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import json

ROOT = Path(__file__).resolve().parents[3]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('case_directory', type=Path)
    folder = parser.parse_args().case_directory.resolve()
    old_dir, new_dir = folder / 'limit_refined', folder / 'completed_basis'
    files = {Path(__file__).resolve()}
    values = []
    for directory in (old_dir, new_dir):
        receipt_path = directory / 'range_two_family_limit_replay.json'
        certificate_path = directory / 'range_two_family_limit_certificate.json'
        receipt = json.loads(receipt_path.read_text())
        if not receipt['accepted'] or not receipt.get('spectator_hopping'):
            raise ValueError('Accepted spectator-family comparison required')
        for name, expected in receipt['source_sha256'].items():
            if hashlib.sha256((ROOT / name).read_bytes()).hexdigest() != expected:
                raise ValueError('Stale family source: ' + name)
        values.append((receipt, json.loads(certificate_path.read_text())))
        files.update((receipt_path, certificate_path))
    old, new = [item[0] for item in values]
    for key in ('kind','half_vector','charged_vector','ratio','theta_half','theta_joint','diagonal_shapes','W','range_two_density_profile'):
        if values[0][1][key] != values[1][1][key]:
            raise ValueError('Fixed family differs: ' + key)
    if F(old['accepted_periodic_lower']) != F(new['accepted_periodic_lower']):
        raise ValueError('Comparison must keep the accepted lower fixed')
    improvement = F(old['periodic_family_upper']) - F(new['periodic_family_upper'])
    if improvement <= 0:
        raise ValueError('No strict ceiling improvement')
    result = {
        'accepted': True,
        'old_ceiling': old['periodic_family_upper'], 'new_ceiling': new['periodic_family_upper'],
        'exact_improvement': str(improvement), 'improvement_float': float(improvement),
        'remaining_gap': new['family_gap'], 'remaining_gap_float': new['family_gap_float'],
        'gap_reduction_fraction': str(improvement / F(old['family_gap'])),
        'source_sha256': {str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(files)},
        'scope': 'Strictly tighter upper limit on lower certificates in the identical fixed spectator family, at the same accepted lower. This is not a physical ground-energy upper or proof of unrestricted numerical optimality.',
    }
    (new_dir / 'ceiling_improvement.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'accepted':True,'improvement':float(improvement),'remaining_gap':new['family_gap_float']}))


if __name__ == '__main__':
    main()
