import json
from pathlib import Path
import numpy as np
from research.global_response_20260913.upper_model.pair_rotation import ansatz, terms, action, energy

ROOT=Path(__file__).resolve().parents[3]
def main():
    d=json.loads((ROOT/'results/certificate_scaling/active_space_ladder/h4/fixture.json').read_text())
    v=ansatz(d['modes'],d['particles'],[.1,-.2,.03,.07]); assert abs(sum(x*x for x in v.values())-1)<1e-12
    ts=terms(d); hv=action(v,ts,d['modes']); direct=0.
    for s,x in v.items(): direct += x*hv.get(s,0.)
    assert abs(direct-energy(v,ts,d['modes']))<1e-13
    assert len(v)<=2**(d['particles'])
    print(json.dumps({'passed':3,'support':len(v),'energy':direct}))
if __name__=='__main__': main()
