"""Exact separation from the entire preceding fixed spectator family."""
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
    if new['kind'] != 'hubbard_projector_extension_v14' or old['kind'] != 'hubbard_projector_extension_v13':
        raise ValueError('Explicit pair extension and preceding spectator certificate required')
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
    folder = parser.parse_args().directory.resolve(); old_folder = folder/'previous_family'
    ep = folder/'range_two_replay.json'; op = old_folder/'range_two_replay.json'
    fp = old_folder/'range_two_family_limit_replay.json'
    energy,old_energy,family_receipt = [read_verified(p) for p in (ep,op,fp)]
    cp = folder/'profile_joint_r1_2_certificate.json'; ocp = old_folder/cp.name
    fcp = old_folder/'range_two_family_limit_certificate.json'
    new,old,family = [json.loads(p.read_text()) for p in (cp,ocp,fcp)]
    match(new,old,family)
    lower = F(energy['lower_replay']['periodic_lower_density'])
    cap = F(family_receipt['periodic_family_upper'])
    if not family_receipt.get('spectator_hopping') or lower <= cap:
        raise ValueError('No strict separation from the full preceding family')
    files = [Path(__file__).resolve(),ep,op,fp,cp,ocp,fcp]
    result = {'accepted':True,'periodic_pair_lower':str(lower),
              'preceding_family_ceiling':str(cap),'exact_separation':str(lower-cap),
              'separation_float':float(lower-cap),
              'lower_improvement':str(lower-F(old_energy['lower_replay']['periodic_lower_density'])),
              'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
              'scope':'New accepted pair-transfer lower strictly exceeds a freshly replayed ceiling for every certificate in the preceding fixed spectator family. Same physical target, size, projector sources, ratio, ceilings and sparse span. No optimality claim for the enlarged family or physical interpretation of the old family ceiling.'}
    (folder/'strict_family_separation.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'accepted':True,'separation':float(lower-cap)}))


if __name__ == '__main__': main()
