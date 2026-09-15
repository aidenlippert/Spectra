"""Exact inherited-SOS terminal diagnosis and the H8 target counterexample.

Imported lower proofs are explicitly identified. This is not new small-factor
discovery and does not use an inherited expectation as a lower bound.
"""
from fractions import Fraction as F
import hashlib
import json
import sys
import time

from experiments.marginal_symbolic import decode
from research.certificate_scaling import wedge_spectral_bound as spectral
from research.certificate_scaling.streaming_reference_upper import upper
from research.compact_response_20260913 import program
from research.composable_response_20260913 import joint
from research.global_response_20260913 import global_program as gp


def lower(data,name):
    start=time.monotonic();meta=json.loads((gp.ROOT/f'results/certificate_scaling/cubic_precision/intervals/final_{name}.json').read_text())
    path=gp.ROOT/meta['certificate'];raw=path.read_bytes();cert=json.loads(raw)
    proofpath=gp.ROOT/meta['proof'];proofraw=proofpath.read_bytes();proof=json.loads(proofraw)
    sha=hashlib.sha256(raw).hexdigest()
    if sha!=meta['certificate_sha256'] or hashlib.sha256(proofraw).hexdigest()!=meta['proof_sha256'] or proof['certificate_sha256']!=sha:
        raise ValueError('Inherited lower binding failed')
    if (data['modes'],data['particles'])!=(cert['modes'],cert['particles']) or decode(data['hamiltonian'],data['modes'],4)!=decode(cert['hamiltonian'],cert['modes'],4):
        raise ValueError('Inherited proof Hamiltonian or physical sector mismatch')
    residual,base=spectral.extract(cert);receipt=spectral.replay_residual(cert,residual,proof)
    return F(receipt['lower']),{'certificate':meta['certificate'],'certificate_sha256':sha,'certificate_bytes':len(raw),
        'residual_proof':meta['proof'],'residual_proof_bytes':len(proofraw),'receipt':receipt,
        'lower_discovery_uses_Hamiltonian_and_sector_only':True,'lower_discovery_enumerates_fixed_N_determinants':False,
        'inherited_full_cubic_Gram_discovery':True,'new_small_factor_discovery':False,'replay_seconds':time.monotonic()-start}


def scalar_terminal_margin(global_lower,b,penalties):
    if any(v<0 for v in penalties):raise ValueError('Negative response allowance')
    gamma=global_lower-b
    if gamma<=0:return {'positive':False,'reason':'Imported global lower does not exceed this target','gamma_Ha':str(gamma)}
    margin=gamma-sum(penalties,F(0))
    return {'positive':margin>0,'global_margin_Ha':str(gamma),'terminal_margin_Ha':str(margin),
        'terminal_margin_mHa':float(margin*1000),'response_penalties_Ha':list(map(str,penalties)),
        'proof':'Each Schur complement of a gamma-positive operator is gamma-positive; subtract each response allowance exactly once in its retained identity metric.'}


def run():
    start=time.monotonic();rows=[]
    for name in ('h6','h8'):
        data,tail,_,oldfirst=joint.load_case(name);L,proof=lower(data,name)
        old_b=F(oldfirst['target_Ha']);entry={'case':name,'inherited_lower_Ha':str(L),'inherited_proof':proof}
        if name=='h6':
            old=json.loads((joint.OUT/'h6_joint_response.json').read_text());oldcheck=joint.check(data,tail,oldfirst,old)
            eta1=F(program.check(data,tail,oldfirst)['exact_program_residual_penalty_Ha'])
            eta2=F(oldcheck['second_residual_penalty_Ha'])
            entry['old_terminal']=scalar_terminal_margin(L,old_b,[eta1,eta2])
            entry['old_terminal']['actual_operator']='K2−eta2 I, where K2 is formed from L1=K1−eta1 I; eta1 is not subtracted again from K2.'
        else:
            meta=json.loads((gp.ROOT/'results/certificate_scaling/cubic_precision/intervals/final_h8.json').read_text())
            ref=json.loads((gp.ROOT/meta['reference']).read_text());U,stats=upper(data,ref['independent_upper'])
            value=U-old_b
            if value>=0:raise ValueError('H8 old-target negative witness not established')
            entry['old_target_obstruction']={'target_Ha':str(old_b),'exact_variational_energy_Ha':str(U),
                'negative_Rayleigh_quotient_of_H_minus_b_Ha':str(value),'negative_expectation_mHa':float(1000*value),
                'old_upper_minus_better_upper_mHa':float(1000*(old_b+F(1,1000)-U)),
                'witness_path':meta['reference'],'witness_states':len(ref['independent_upper']['states']),
                'upper_replay':stats,'scope':'Exact physical negative witness for the old full-H target. It does not prove a negative direction inside the eliminated joint sector.'}
        new=json.loads((gp.OUT/f'{name}_global.json').read_text());r=gp.check(data,tail,new)
        if r['status']!='certified_response':raise ValueError('New response unavailable')
        entry['global_response']=r
        entry['new_terminal']=scalar_terminal_margin(L,F(new['target_Ha']),[F(r['residual_penalty_Ha'])])
        rows.append(entry)
    forbidden=[x for x in ('numpy','scipy','cvxpy','pyscf') if x in sys.modules]
    if forbidden:raise AssertionError('Numerical import in exact reference diagnosis')
    result={'cases':rows,'wall_seconds':time.monotonic()-start,'numerical_packages_loaded':forbidden,
        'scope':'Reference-assisted positivity and target diagnosis using independently replayed existing operator SOS lower proofs. No new small-terminal proof is inferred.'}
    (gp.OUT/'reference_diagnosis.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'seconds':result['wall_seconds'],'terminal_margins_mHa':{r['case']:r['new_terminal'].get('terminal_margin_mHa') for r in rows}},indent=2))


if __name__=='__main__':run()
