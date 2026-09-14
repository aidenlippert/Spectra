"""Verify exact frozen transfer and compare accepted pair/no-pair thresholds."""
from pathlib import Path
from fractions import Fraction as F
import copy
import hashlib
import json
from pair_transfer_separation import read_verified

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT/'results/marginal_graded_hubbard8/pair_transfer'


def main():
    source = BASE/'W_zero/profile_joint_r1_2_certificate.json'
    cp = BASE/'held_out/profile_joint_r1_2_certificate.json'
    ap = BASE/'held_out/without_pair/profile_joint_r1_2_certificate.json'
    ep = cp.with_name('range_two_replay.json'); aep = ap.with_name('range_two_replay.json')
    source_receipt = source.with_name('range_two_replay.json')
    original,new,ablation = [json.loads(p.read_text()) for p in (source,cp,ap)]
    read_verified(source_receipt); energy,no_pair = [read_verified(p) for p in (ep,aep)]
    expected = copy.deepcopy(original)
    target = new['target']; old_target = original['target']
    if target != {'U':'5','t':'1','V':'1/4','W':'-1/5'}:
        raise ValueError('Explicit held-out target required')
    expected['target'] = target
    local = expected['local_window']
    local['U'] = str(F(5,6)*F(target['U']))
    local['onsite_profile'] = [str(F(v)*F(target['U'])/F(old_target['U'])) for v in local['onsite_profile']]
    local['V'] = target['V']
    local['density_profile'] = [str(F(v)*F(target['V'])/F(old_target['V'])) for v in local['density_profile']]
    local['range_two_density_profile'] = [str(F(v)+F(5,4)*(F(target['W'])-F(old_target['W']))) for v in local['range_two_density_profile']]
    expected['penalized_lower'] = new['penalized_lower']
    if expected != new: raise ValueError('Frozen source, correction or profile recipe changed')
    expected_ablation = copy.deepcopy(new)
    expected_ablation.pop('pair_transfer'); expected_ablation['kind']='hubbard_projector_extension_v13'
    expected_ablation['penalized_lower'] = ablation['penalized_lower']
    if expected_ablation != ablation: raise ValueError('Ablation changed other coefficients')
    lower = F(energy['lower_per_site']); old_lower = F(no_pair['lower_per_site'])
    if energy['upper_per_site'] != no_pair['upper_per_site']:
        raise ValueError('Physical trial-state upper changed')
    files = [Path(__file__).resolve(),Path(__file__).with_name('pair_transfer_separation.py'),
             source,cp,ap,ep,aep,source_receipt]
    result = {'accepted':True,'target':target,'with_pair_open_lower':str(lower),
              'without_pair_open_lower':str(old_lower),'exact_improvement':str(lower-old_lower),
              'improvement_float':float(lower-old_lower),'physical_upper':energy['upper_per_site'],
              'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
              'scope':'Frozen projector sources, penalties and correction coefficients transferred to a new U,V,W target on the same chain geometry. Only physical profiles and scalar local threshold changed. Comparison is between two accepted fixed-coefficient certificates; it is not a ceiling on a reoptimized no-pair family. No generic molecular, geometry or scalability claim.'}
    (BASE/'held_out/frozen_transfer_comparison.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'accepted':True,'improvement':float(lower-old_lower)}))


if __name__ == '__main__': main()
