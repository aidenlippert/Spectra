"""Independent bounded physical-sector controls; excluded from discovery."""
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json
import time

from experiments.marginal_symbolic import add,scale
from research.collective_interference_20260913.validation import sparse_matrix
from research.certificate_scaling.streaming_reference_upper import upper
from research.molecular_collective_20260913.core import extract,retained_polynomial,tail_replay

ROOT=Path(__file__).resolve().parents[2]


def run(campaign,out):
    import numpy as np
    from scipy.sparse.linalg import eigsh
    summary=json.loads((campaign/'summary.json').read_text());rows=[];start=time.monotonic()
    def extreme(matrix,which):
        vals,vecs=eigsh(matrix,k=1,which=which,tol=1e-11,v0=np.random.default_rng(20260913).normal(size=matrix.shape[0]))
        residual=float(np.linalg.norm(matrix@vecs[:,0]-vals[0]*vecs[:,0]))
        if residual>1e-8:raise AssertionError('Control eigenvector residual')
        return float(vals[0]),residual
    for system in summary['systems'][:2]:
        t=time.monotonic();name=system['name'];directory=campaign/name
        data=json.loads((directory/'fixture.json').read_text())
        if data['modes']>12:raise ValueError('Full-sector validation capped at twelve modes')
        cert=json.loads((directory/f"rank_{system['best_rank']}/tail.json").read_text())
        p=extract(data);retained=retained_polynomial(p,cert);receipt=tail_replay(data,cert)
        matrices={};actions=0
        for label,h in [('original',p['h']),('retained',retained),('difference',add(p['h'],scale(retained,-1)))]:
            matrix,basis,work=sparse_matrix(h,p['modes'],p['particles']);actions+=work;matrices[label]=matrix
        lo,reslo=extreme(matrices['difference'],'SA');hi,reshi=extreme(matrices['difference'],'LA')
        original,resoriginal=extreme(matrices['original'],'SA');small,ressmall=extreme(matrices['retained'],'SA')
        certified_lo=float(F(receipt['lower_operator_shift_Ha']));certified_hi=float(F(receipt['upper_operator_shift_Ha']))
        if not certified_lo-1e-8<=lo<=hi<=certified_hi+1e-8:raise AssertionError('Full operator spectrum escaped certified envelope')
        if not certified_lo-1e-8<=original-small<=certified_hi+1e-8:raise AssertionError('Ground-energy difference escaped certificate')
        refpath=ROOT/f'results/certificate_scaling/active_space_ladder_references_aligned/{name}/upper.json'
        reference=json.loads(refpath.read_text());ref,refstats=upper(data,reference['independent_upper'])
        if ref!=F(reference['upper']):raise ValueError('Reference upper mismatch')
        rows.append({'system':name,'rank':system['best_rank'],'full_sector_dimension':len(basis),
                     'full_word_state_actions':actions,'numeric_difference_spectrum_extrema_Ha':[lo,hi],
                     'numeric_original_ground_Ha':original,'numeric_retained_ground_Ha':small,
                     'numeric_ground_shift_mHa':1000*(original-small),
                     'eigenvector_residuals':[reslo,reshi,resoriginal,ressmall],
                     'certified_operator_interval_Ha':[receipt['lower_operator_shift_Ha'],receipt['upper_operator_shift_Ha']],
                     'reference_upper':refstats,'reference_path':str(refpath.relative_to(ROOT)),
                     'reference_sha256':hashlib.sha256(refpath.read_bytes()).hexdigest(),
                     'wall_seconds':time.monotonic()-t})
        print(name,'operator interval',lo,hi,'ground shift mHa',1000*(original-small),flush=True)
    result={'rows':rows,'wall_seconds':time.monotonic()-start,'scope':'Full-sector numerical validation only. Reference rational upper replay exact; no physical full-sector work is part of pattern discovery or accepting tail replay.'}
    out.write_text(json.dumps(result,indent=2)+'\n');return result


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--campaign',type=Path,required=True);parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args();run(args.campaign,args.out)
