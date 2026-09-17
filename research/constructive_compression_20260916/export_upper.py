"""Round a numerical charge MPS; rounding is proposed, never trusted as a bound."""
import argparse,json
from pathlib import Path
from research.constructive_compression_20260916.model import fixture
from research.molecular_collective_20260913.core import digest


def export(proposal,rungs=4,bits=32):
    data=fixture(rungs);den=1<<bits;layers=[]
    for a in proposal['arrays']:
        layer=[]
        for l,row in enumerate(a):
            for s,rr in enumerate(row):
                for r,value in enumerate(rr):
                    v=round(value*den)
                    if v:layer.append([l,s,r,v])
        layers.append(layer)
    return data,{'kind':'integer_charge_mps_v1','fixture_sha256':digest(data),
      'modes':data['modes'],'particles':data['particles'],'spin_counts':data['spin_counts'],
      'denominator':den,'bond_charges':proposal['charges'],'tensors':layers}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('proposal',type=Path);p.add_argument('out',type=Path)
    p.add_argument('--rungs',type=int,default=4);p.add_argument('--bits',type=int,default=32);a=p.parse_args()
    data,cert=export(json.loads(a.proposal.read_text()),a.rungs,a.bits)
    a.out.mkdir(parents=True,exist_ok=False)
    (a.out/'fixture.json').write_text(json.dumps(data,indent=2)+'\n')
    (a.out/'state.json').write_text(json.dumps(cert)+'\n')
    print(json.dumps({'status':'proposed_exact_mps','entries':sum(map(len,cert['tensors'])),
         'maximum_bond':max(map(len,cert['bond_charges']))}))
