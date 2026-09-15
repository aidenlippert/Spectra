"""Export executable retained models and check the supporting-family obstruction."""
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json
import time

from experiments.marginal_symbolic import verify,decode,encode,add,scale,product
from research.molecular_collective_20260913.core import (
    extract,digest,tail_replay,support_replay,factor_operators,retained_polynomial,one_matrix,matmul)

ROOT=Path(__file__).resolve().parents[2]


def run(campaign,out):
    start=time.monotonic();campaign=campaign.resolve();summary=json.loads((campaign/'summary.json').read_text());models=[]
    for system in summary['systems']:
        directory=campaign/system['name'];data=json.loads((directory/'fixture.json').read_text())
        tail=json.loads((directory/f"rank_{system['best_rank']}/tail.json").read_text());p=extract(data)
        receipt=tail_replay(data,tail)
        model={'kind':'fixed_N_density_square_model_v1','modes':p['modes'],'particles':p['particles'],
               'hamiltonian_definition':'one_body + (1/2) sum_k weight_k (sum_a vector_ka Q_a)^2',
               'feature_definition':'For p=q, Q_a=sum_spin a_p,spin^dagger a_p,spin. For p<q add its Hermitian transpose. Spin orbitals are interleaved (2p,2p+1).',
               'features':[list(pair) for pair in p['pairs']],'one_body':encode(p['one']),
               'factors':tail['factors'],'source_fixture_sha256':digest(data),
               'lower_operator_shift_Ha':receipt['lower_operator_shift_Ha'],
               'upper_operator_shift_Ha':receipt['upper_operator_shift_Ha'],
               'scope':'Retained model and operator sandwich valid in this fixed-N sector; full Hamiltonian is required to replay its source certificate.'}
        restored=decode(model['one_body'],p['modes'],2)
        for weight,q in factor_operators(p,tail):restored=add(restored,scale(product(q,q),weight/2))
        if restored!=retained_polynomial(p,tail):raise ValueError('Retained model export changed its polynomial')
        path=directory/'retained_model.json';raw=json.dumps(model,indent=2)+'\n';path.write_text(raw)
        models.append({'system':system['name'],'path':str(path.relative_to(ROOT)),
                       'bytes':len(raw.encode()),'compact_bytes':len(json.dumps(model,separators=(',',':')).encode()),
                       'sha256':hashlib.sha256(raw.encode()).hexdigest()})
    # A previously independent exact lower bound makes the dual obstruction
    # a gap to the true ground energy, not just a gap to an approximate upper.
    directory=campaign/'h6';data=json.loads((directory/'fixture.json').read_text())
    prior_path=ROOT/'results/molecular_identity_20260913/h6_coupled/round_1/certificate.json'
    prior=json.loads(prior_path.read_text())
    if (prior['modes'],prior['particles'])!=(data['modes'],data['particles']) or decode(prior['hamiltonian'],12,4)!=decode(data['hamiltonian'],12,4):
        raise ValueError('Prior lower proof has another Hamiltonian or sector')
    t=time.monotonic();lower=verify(prior);prior_seconds=time.monotonic()-t;obstructions=[]
    for rank in (10,21):
        folder=directory/f'rank_{rank}';tail=json.loads((folder/'tail.json').read_text());support=json.loads((folder/'support.json').read_text())
        result=support_replay(data,tail,support);gap=F(lower['lower'])-F(result['support_family_ceiling_Ha'])
        if gap<=0:raise ValueError('No certified supporting-family obstruction')
        obstructions.append({'requested_rank':rank,'actual_factors':result['factors'],
                             'certified_true_ground_lower_Ha':lower['lower'],
                             'exact_support_family_ceiling_Ha':result['support_family_ceiling_Ha'],
                             'true_ground_minus_any_support_lower_at_least_Ha':str(gap),
                             'support_path':str((folder/'support.json').relative_to(ROOT)),
                             'support_sha256':hashlib.sha256((folder/'support.json').read_bytes()).hexdigest(),
                             'tail_path':str((folder/'tail.json').relative_to(ROOT)),
                             'tail_sha256':hashlib.sha256((folder/'tail.json').read_bytes()).hexdigest()})
    p=extract(data);tail=json.loads((directory/'rank_10/tail.json').read_text())
    ops=[one_matrix(q,12) for _,q in factor_operators(p,tail)];largest=(F(0),None)
    for i,a in enumerate(ops):
        for j in range(i+1,len(ops)):
            ab=matmul(a,ops[j]);ba=matmul(ops[j],a)
            for row in range(12):
                for col in range(12):
                    c=ab[row][col]-ba[row][col]
                    if abs(c)>largest[0]:largest=(abs(c),{'patterns':[i,j],'one_body_entry':[row,col],'exact_commutator_coefficient':str(c)})
    result={'models':models,'prior_lower_source':str(prior_path.relative_to(ROOT)),
            'prior_lower_sha256':hashlib.sha256(prior_path.read_bytes()).hexdigest(),
            'prior_lower_replay_seconds':prior_seconds,'obstructions':obstructions,
            'noncommuting_pattern_witness':largest[1],'wall_seconds':time.monotonic()-start,
            'scope':'Exact obstruction for linear supports of the supplied square factors plus the stated scalar tail bound. Does not exclude joint fermionic constraints, alternative factors, or general compression.'}
    out.write_text(json.dumps(result,indent=2)+'\n');print([(r['requested_rank'],float(F(r['true_ground_minus_any_support_lower_at_least_Ha']))) for r in obstructions]);return result


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--campaign',type=Path,required=True);parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args();run(args.campaign,args.out)
