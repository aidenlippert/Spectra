"""Fresh standard-library replay of the compact lower/family witnesses."""
import argparse
import copy
from fractions import Fraction as F
import json
from pathlib import Path
from research.local_response_20260916.local_exact import verify,verify_dual
from research.local_response_20260916.projector_correction import verify_correction


def replay(certificate_path,correction_path,out,upper_path=None):
    if out.exists():raise FileExistsError(out)
    out.mkdir(parents=True)
    certificate=json.loads(certificate_path.read_text())
    correction=json.loads(correction_path.read_text())
    receipt={'primal':verify(certificate),'dual':verify_dual(certificate),
             'collective_correction':verify_correction(certificate,correction)}
    if upper_path is not None:
        from research.local_response_20260916.reference import check_upper
        upper=check_upper(json.loads(upper_path.read_text()))
        lo=F(receipt['collective_correction']['energy_lower_over_t'])
        hi=F(upper['energy_upper_over_t'])
        if lo>hi:raise AssertionError('Lower exceeds independently checked upper')
        receipt['physical_upper']=upper
        receipt['complete_interval']={'lower_over_t':str(lo),'upper_over_t':str(hi),
            'width_over_t':str(hi-lo),'width_float':float(hi-lo),
            'target_width_over_t':'1/1000','meets_target':hi-lo<=F(1,1000),
            'enumeration_free_complete_pipeline':False}
    rejected=[]
    bad=copy.deepcopy(certificate);bad['messages'][0][1][4]='13'
    tests=[('asymmetric_boundary',lambda:verify(bad))]
    for name,operation in tests:
        try:operation()
        except ValueError:rejected.append(name)
        else:raise AssertionError(name+' was accepted')
    bad=copy.deepcopy(certificate)
    a=bad['dual_marginals'][0]['0,0'];a[0][0]=str(F(a[0][0])+1)
    try:verify_dual(bad)
    except ValueError:rejected.append('changed_marginal_trace')
    else:raise AssertionError('Corrupt dual accepted')
    bad=copy.deepcopy(correction);bad['correction_over_t']='1'
    try:verify_correction(certificate,bad)
    except ValueError:rejected.append('overclaimed_collective_gain')
    else:raise AssertionError('Overclaimed correction accepted')
    bad=copy.deepcopy(correction);bad['local_profiles'][0]['integer_vector']=[0]*256
    try:verify_correction(certificate,bad)
    except ValueError:rejected.append('zero_local_vector')
    else:raise AssertionError('Zero projector accepted')
    receipt['corruptions_rejected']=rejected
    (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--certificate',type=Path,required=True)
    p.add_argument('--correction',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--upper',type=Path)
    a=p.parse_args();replay(a.certificate,a.correction,a.out,a.upper)
