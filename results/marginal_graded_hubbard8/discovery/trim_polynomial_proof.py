"""Trim tiny positive atoms using a conservative coefficient-L1 allowance.

Constant positive atoms first become an explicit bound. Discarded atoms spend
at most one quarter of the certified post-residual positivity margin. The
unchanged exact verifier must accept the resulting certificate.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse,json,hashlib,time
from experiments.marginal_polynomial_metric import replay


def trim(source,out):
    source=Path(source);raw=source.read_bytes();certificate=json.loads(raw);started=time.monotonic()
    baseline=replay(certificate);changes={};sites=certificate['modes']//2
    for part in ('weight','numerator'):
        proof=certificate[part+'_proof'];den=proof['denominator']
        constants=sum(a['weight'] for a in proof['positive_indicators'] if a['required']==0)
        proof['positive_indicators']=[a for a in proof['positive_indicators'] if a['required']!=0]
        proof['bound']+=constants
        margin=F(baseline[part+'_positivity']['lower'])+F(constants,den)
        allowance=int(margin*den/4);spent=0;options=[];removed=set()
        for family in ('positive_indicators','charge_indicators'):
            for index,atom in enumerate(proof[family]):
                cost=atom['weight']*(1<<(atom['required']^atom['occupied']).bit_count())*(sites+1 if family=='charge_indicators' else 1)
                options.append((cost,family,index))
        for cost,family,index in sorted(options):
            if spent+cost>allowance:break
            spent+=cost;removed.add((family,index))
        before={family:len(proof[family]) for family in ('positive_indicators','charge_indicators')}
        for family in before:proof[family]=[a for index,a in enumerate(proof[family]) if (family,index) not in removed]
        changes[part]={'constant_promoted':str(F(constants,den)),'margin_before_trimming':str(margin),'conservative_discard_l1':str(F(spent,den)),'atoms_discarded':len(removed),'before':before,'after':{family:len(proof[family]) for family in before}}
    receipt=replay(certificate);out=Path(out);out.mkdir(parents=True,exist_ok=True)
    diagnostic={'source_sha256':hashlib.sha256(raw).hexdigest(),'changes':changes,'seconds':time.monotonic()-started,'scope':'Unchanged Hamiltonian, metric and gamma. Constant-atom promotion and conservative coefficient-L1 trimming only; number identities retained and full original exact replay required.'}
    for name,value in [('certificate',certificate),('receipt',receipt),('trimming',diagnostic)]:
        (out/(name+'.json')).write_text(json.dumps(value,indent=2)+'\n')
    print(json.dumps(diagnostic),flush=True);return receipt

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();trim(a.source,a.out)
