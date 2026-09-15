"""Measure exact certificate size/error tradeoffs from sparse factor rounding.

Unlike a dense spectral rotation, this preserves the original factor basis.
Every result is accepted only through the existing rational CAR checker.
"""
import argparse
import copy
from fractions import Fraction as F
import json
from pathlib import Path
import sys
import time
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from experiments.marginal_symbolic import verify


def rounded_integer(n, scale, old_scale):
    return round(F(n*scale,old_scale))


def run(fixture, output):
    source=json.loads(fixture.read_text())
    base=verify(source)
    output.mkdir(parents=True,exist_ok=True)
    results=[]
    settings=[(int(source['denominator']),F(0))]+[(10**p,F(0)) for p in (3,4,5,6,7)]
    settings += [(10**6,F(1,10**p)) for p in (3,4,5)]
    for den,threshold in settings:
        start=time.monotonic()
        candidate=copy.deepcopy(source)
        candidate['denominator']=den
        for block in candidate['blocks']:
            rows=[]
            for original in block['factor']:
                row=[rounded_integer(x,den,int(source['denominator']))
                     if abs(F(x,int(source['denominator'])))>=threshold else 0 for x in original]
                if any(row): rows.append(row)
            used=[i for i in range(len(block['words'])) if any(row[i] for row in rows)]
            block['words']=[block['words'][i] for i in used]
            block['factor']=[[row[i] for i in used] for row in rows]
        candidate['blocks']=[b for b in candidate['blocks'] if b['factor']]
        receipt=verify(candidate)
        name=f'd{den}_t{threshold.numerator}-{threshold.denominator}.json'
        path=output/name
        path.write_text(json.dumps(candidate,separators=(',',':'))+'\n')
        receipt.update(certificate=name,bytes=path.stat().st_size,
                       denominator=den,threshold=str(threshold),
                       loss_vs_original=str(F(base['lower'])-F(receipt['lower'])),
                       seconds=time.monotonic()-start)
        results.append(receipt)
    (output/'receipts.json').write_text(json.dumps(results,indent=2)+'\n')
    return results


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--fixture',type=Path,required=True)
    p.add_argument('--outputdir',type=Path,required=True)
    a=p.parse_args()
    for r in run(a.fixture,a.outputdir):
        print(r['certificate'],r['bytes'],float(F(r['loss_vs_original'])))
