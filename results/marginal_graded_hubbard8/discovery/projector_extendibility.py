"""Standard-library replay of compact quantitative extension certificates."""
from pathlib import Path
from fractions import Fraction as F
import hashlib
import json
import sys
import time

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT))
from experiments.marginal_projector_extendibility import replay
from experiments.marginal_window_family_bound import replay as replay_family

OUT = ROOT/'results/marginal_graded_hubbard8/projector_extendibility'


def main():
    OUT.mkdir(exist_ok=True)
    family_path = ROOT/'results/marginal_graded_hubbard8/weighted_window_family/certificate.json'
    original = json.loads(family_path.read_text())
    family = replay_family(original)
    old_ceiling = F(family['maximum_local_minimum_upper'])/3
    certificates = []
    for m,ceiling in ((2,'1.380638468'),(3,'1.901418483'),(4,'2.223786409')):
        certificates.append({'kind':'hubbard_projector_extension_v1',
            'a':original['a'],'b':original['b'],'vector':original['upper_vector'],
            'windows':m,'projector_sum_ceiling':ceiling,'penalty':'0.3157',
            'penalized_lower':'-1.724725','chain_sites':1000000})
    certificates.append(dict(certificates[-1],a='0.5613',b='0.8171',
                             penalty='0.4386',penalized_lower='-1.62983'))
    (OUT/'certificates.json').write_text(json.dumps(certificates,indent=2)+'\n')
    start = time.monotonic()
    receipts = [replay(c) for c in certificates]
    for r in receipts:
        improvement = F(r['periodic_lower_density'])-old_ceiling
        if improvement <= 0:
            raise ValueError('Proposed extension does not break the old family ceiling')
        r['improvement_over_four_site_boundary_family_ceiling'] = str(improvement)
    sources = ['experiments/marginal_projector_extendibility.py',
        'experiments/marginal_window_family_bound.py',
        'experiments/marginal_local_hubbard_block.py',
        'experiments/marginal_transfer_verify.py',
        'experiments/marginal_polynomial_sos.py',
        'results/marginal_graded_hubbard8/discovery/projector_extendibility.py',
        'results/marginal_graded_hubbard8/weighted_window_family/certificate.json',
        'results/marginal_graded_hubbard8/projector_extendibility/certificates.json']
    result = {'accepted':True,'seconds':time.monotonic()-start,
        'four_site_boundary_family_ceiling_density':str(old_ceiling),
        'certificates':receipts,
        'source_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sources},
        'scope':'Fresh exact local and projector Gram positivity. Improves the four-site boundary-correction ceiling; remains weaker than the separately certified six-site energy lower bound. No optimizer convergence, family optimality, or general representability claim.'}
    (OUT/'independent_replay.json').write_text(json.dumps(result,indent=2)+'\n')
    for r in receipts:
        print(json.dumps({'windows':r['overlap']['windows'],
            'fidelity_ceiling':r['overlap']['average_fidelity_ceiling'],
            'periodic_lower_density':float(F(r['periodic_lower_density'])),
            'open_lower_density':float(F(r['open_lower_density'])),
            'maximum_gram_psd':r['overlap']['maximum_psd_dimension']}),flush=True)


if __name__ == '__main__':
    main()
