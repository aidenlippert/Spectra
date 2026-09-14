"""Exact matched-target fixed-recipe contribution of both new operators."""
from pathlib import Path
from fractions import Fraction as F
import argparse,copy,hashlib,json
from pair_transfer_separation import read_verified

ROOT=Path(__file__).resolve().parents[3]

def main():
    parser=argparse.ArgumentParser();parser.add_argument('directory',type=Path)
    directory=parser.parse_args().directory.resolve()
    cp=directory/'profile_joint_r1_2_certificate.json';ap=directory/'without_coherent/profile_joint_r1_2_certificate.json'
    rp=directory/'range_two_replay.json';arp=ap.with_name('range_two_replay.json')
    energy,ablation=[read_verified(p) for p in (rp,arp)]
    c,a=[json.loads(p.read_text()) for p in (cp,ap)]
    if c['kind']!='hubbard_projector_extension_v17':raise ValueError('New energy version required')
    expected=copy.deepcopy(c);expected.pop('coherent_projector');expected['kind']='hubbard_projector_extension_v16';expected['penalized_lower']=a['penalized_lower']
    if expected!=a or energy['upper_per_site']!=ablation['upper_per_site']:raise ValueError('Ablation changed other coefficients or target')
    difference=F(energy['lower_per_site'])-F(ablation['lower_per_site'])
    files=[Path(__file__).resolve(),Path(__file__).with_name('pair_transfer_separation.py'),cp,ap,rp,arp]
    result={'accepted':True,'with_new_open_lower':energy['lower_per_site'],'without_new_open_lower':ablation['lower_per_site'],
            'exact_fixed_recipe_contribution':str(difference),'fixed_recipe_contribution_float':float(difference),
            'new_terms_help_this_fixed_recipe':difference>0,
            'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
            'scope':'Signed difference of accepted matched-target energy certificates, with only both coherent-projector coefficients removed and scalar threshold recomputed. A positive value proves a benefit for this fixed recipe, not separation from a reoptimized older family.'}
    (directory/'ablation_comparison.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'accepted':True,'fixed_recipe_contribution':float(difference)}))

if __name__=='__main__':main()
