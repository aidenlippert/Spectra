"""Exact lower improvement and signed comparison with the preceding family ceiling."""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import json

ROOT = Path(__file__).resolve().parents[3]


def read_verified(path):
    result = json.loads(path.read_text())
    if result.get('accepted') is not True: raise ValueError('Accepted receipt required')
    for source,expected in result['source_sha256'].items():
        if hashlib.sha256((ROOT/source).read_bytes()).hexdigest() != expected:
            raise ValueError('Stale proof source: '+source)
    return result


def match(new, old, family):
    if new['kind'] != 'hubbard_projector_extension_v17' or old['kind'] != 'hubbard_projector_extension_v16':
        raise ValueError('Explicit coherent-projector extension and preceding three-spectator certificate required')
    if new['target'] != old['target'] or new['chain_sites'] != old['chain_sites']:
        raise ValueError('Target or size changed')
    for key in ('vector','windows','projector_sum_ceiling'):
        if new[key] != old[key]: raise ValueError('Fixed half projector changed')
    for key in ('vector','windows','ratio','projector_sum_ceiling'):
        if new['joint'][key] != old['joint'][key]: raise ValueError('Fixed joint projector changed')
    shapes = family['diagonal_shapes']
    for c in (new,old):
        diagonal = {int(s):F(v) for s,v in c['telescoping_diagonal'].items()}
        rebuilt = {}
        for shape in shapes:
            shape = {int(s):F(v) for s,v in shape.items()}
            if set(shape)&set(rebuilt): raise ValueError('Overlapping sparse shapes')
            ratios = {diagonal.get(s,F(0))/v for s,v in shape.items() if v}
            if len(ratios) != 1: raise ValueError('Changed sparse correction span')
            ratio = ratios.pop(); rebuilt.update({s:ratio*v for s,v in shape.items()})
        if {s:v for s,v in rebuilt.items() if v} != diagonal:
            raise ValueError('Uncovered sparse correction')


def main():
    parser = argparse.ArgumentParser();parser.add_argument('directory',type=Path)
    parser.add_argument('--previous-directory',type=Path)
    args=parser.parse_args();folder=args.directory.resolve()
    old_folder=args.previous_directory.resolve() if args.previous_directory else folder/'previous_family'
    ep = folder/'range_two_replay.json'; op = old_folder/'range_two_replay.json'
    fp = old_folder/'range_two_family_limit_replay.json'
    energy,old_energy,family_receipt = [read_verified(p) for p in (ep,op,fp)]
    cp = folder/'profile_joint_r1_2_certificate.json'; ocp = old_folder/cp.name
    fcp = old_folder/'range_two_family_limit_certificate.json'
    new,old,family = [json.loads(p.read_text()) for p in (cp,ocp,fcp)]
    match(new,old,family)
    lower = F(energy['lower_replay']['periodic_lower_density'])
    cap = F(family_receipt['periodic_family_upper'])
    if not family_receipt.get('three_spectator_hopping'):
        raise ValueError('Preceding three-spectator family cap required')
    if lower <= F(old_energy['lower_replay']['periodic_lower_density']):
        raise ValueError('No strict improvement over the preceding lower certificate')
    files = [Path(__file__).resolve(),ep,op,fp,cp,ocp,fcp]
    result = {'accepted':True,'strict_family_separation':lower>cap,'periodic_coherent_projector_lower':str(lower),
              'preceding_family_ceiling':str(cap),'exact_separation':str(lower-cap),
              'separation_float':float(lower-cap),
              'lower_improvement':str(lower-F(old_energy['lower_replay']['periodic_lower_density'])),
              'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
              'scope':'Exact lower improvement at fixed target, size, projector sources, ratio, ceilings and sparse span. A positive signed separation proves the new lower exceeds the preceding full three-spectator family ceiling; a nonpositive value does not establish such separation. No enlarged-family optimality or physical interpretation of family ceilings.'}
    (folder/'previous_family_comparison.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'accepted':True,'separation':float(lower-cap)}))


if __name__ == '__main__': main()
