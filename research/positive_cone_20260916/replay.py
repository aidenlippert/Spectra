"""Independent process replay, standard library only, no discovery solver."""
import argparse
import copy
import json
from pathlib import Path
from research.positive_cone_20260916.cone import verify


def replay(certificate,out):
    if out.exists():raise FileExistsError(out)
    payload=json.loads(certificate.read_text())
    receipt=verify(payload)
    rejected=[]
    for change in ('wrong_U','wrong_charge','asymmetric_C','singular_C','overclaimed_lower'):
        bad=copy.deepcopy(payload)
        if change=='wrong_U':bad['U']='9'
        elif change=='wrong_charge':bad['target_spin_populations']=[3,5]
        elif change=='asymmetric_C':bad['integer_C'][0][1]+=1
        elif change=='singular_C':bad['integer_C']=[[0]*70 for _ in range(70)]
        else:bad['lower_over_t']='-3'
        try:verify(bad)
        except ValueError:rejected.append(change)
        else:raise AssertionError('Corrupted certificate accepted: '+change)
    receipt['corrupted_certificates_rejected']=rejected
    out.mkdir(parents=True)
    (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--certificate',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();replay(a.certificate,a.out)
