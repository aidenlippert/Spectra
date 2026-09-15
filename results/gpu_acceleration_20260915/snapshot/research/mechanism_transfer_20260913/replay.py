"""Fresh exact transfer bounds with the frozen recipe and references bound."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import sys
import time

from research.mechanism_transfer_20260913.core import replay
from research.certificate_scaling.streaming_reference_upper import upper

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'results/mechanism_transfer_20260913'


def run():
    start=time.monotonic();records=[];cases=[]
    rule_hash=hashlib.sha256((OUT/'frozen_rule.json').read_bytes()).hexdigest()
    for directory in sorted((OUT/'campaign').iterdir()):
        if not directory.is_dir() or not (directory/'summary.json').exists():continue
        meta=json.loads((directory/'summary.json').read_text())
        if meta['frozen_rule_sha256']!=rule_hash:raise ValueError('Recipe binding changed')
        data=json.loads((directory/'fixture.json').read_text());tail=json.loads((directory/'tail.json').read_text())
        reference=json.loads((directory/'upper.json').read_text())
        U,ur=upper(data,reference['independent_upper'])
        if U!=F(meta['upper_Ha']) or U!=F(reference['upper']):raise ValueError('Upper reference changed')
        span=json.loads((directory/'generated_span.json').read_text())
        case_rows=[]
        for row in meta['trials']:
            if not row['accepted']:continue
            path=directory/row['arm']/'certificate.json';cert=json.loads(path.read_text())
            if row['arm']!='quadratic':
                ids=sorted(set(e['group'] for e in span))
                if len(ids)!=len(cert['anti_blocks']):raise ValueError('Recipe block count differs')
                for gid,b in zip(ids,cert['anti_blocks']):
                    if b['directions']!=[e['vector'] for e in span if e['group']==gid]:
                        raise ValueError('Accepted directions differ from the frozen recipe output')
            accepted=replay(data,tail,cert);L=F(accepted['original_lower_Ha'])
            if L!=F(row['lower_Ha']) or U-L!=F(row['width_Ha']):raise ValueError('Transfer interval changed')
            if U<L:raise AssertionError('Invalid interval ordering')
            item={'case':meta['case'],'directory':str(directory.relative_to(ROOT)),
                'arm':row['arm'],'certificate_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
                'lower_Ha':str(L),'upper_Ha':str(U),'width_Ha':str(U-L),'width_mHa':float(1000*(U-L)),
                'directions':row['directions'],'inside_budget':row['inside_budget'],'replay':accepted}
            records.append(item);case_rows.append(item)
            print(json.dumps({k:v for k,v in item.items() if k!='replay'}),flush=True)
        accepted_rows=[r for r in case_rows if r['inside_budget']]
        best=min(accepted_rows,key=lambda r:F(r['width_Ha']))
        by_arm={r['arm']:r for r in case_rows};stats={'case':meta['case'],'upper_replay':ur,'best_arm':best['arm'],
            'best_width_mHa':best['width_mHa']}
        if all(a in by_arm for a in ('quadratic','coupled','separate')):
            stats.update(coupled_gain_over_quadratic_mHa=by_arm['quadratic']['width_mHa']-by_arm['coupled']['width_mHa'],
                coupling_gain_over_separate_mHa=by_arm['separate']['width_mHa']-by_arm['coupled']['width_mHa'])
        cases.append(stats)
    forbidden=[name for name in ('numpy','scipy','cvxpy','pyscf') if name in sys.modules]
    if forbidden:raise AssertionError('Numerical package loaded during accepting transfer replay')
    result={'rows':records,'cases':cases,'wall_seconds':time.monotonic()-start,'numerical_packages_loaded':forbidden,
        'scope':'Exact finite-basis interval proofs for each case. No held-out family ceiling or rigorous FCI ground lower is claimed.'}
    (OUT/'fresh_replay.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':run()
