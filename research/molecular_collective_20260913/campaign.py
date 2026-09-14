"""Measured molecular pattern frontiers, supporting bounds, and exact ceilings."""
from fractions import Fraction as F
from pathlib import Path
import argparse
import json
import time

from research.molecular_collective_20260913.core import (
    extract,propose_tail,tail_replay,propose_support,support_replay,propose_rank,rank_replay)

ROOT=Path(__file__).resolve().parents[2]


def write(path,data):
    raw=json.dumps(data,separators=(',',':'))+'\n';path.write_text(raw);return len(raw.encode())


def run(out):
    start=time.monotonic();out.mkdir(parents=True,exist_ok=False);rows=[];systems=[]
    for name in ('h4','h6','h8','h10'):
        t=time.monotonic();source=ROOT/f'results/certificate_scaling/active_space_ladder/{name}/fixture.json'
        data=json.loads(source.read_text());directory=out/name;directory.mkdir();write(directory/'fixture.json',data)
        p=extract(data);d=len(p['basis']);s=p['spatial'];target_rank=2*s-2
        ranks=list(range(d+1)) if s<=6 else [0,s,target_rank-1,target_rank,target_rank+1,d]
        local=[]
        for r in ranks:
            tt=time.monotonic();tail,discovery=propose_tail(data,p,r);receipt=tail_replay(data,tail)
            case=directory/f'rank_{r}';case.mkdir();size=write(case/'tail.json',tail)
            row={'system':name,'requested_rank':r,'tail':receipt,'discovery':discovery,
                 'tail_certificate_bytes':size,'case_seconds':time.monotonic()-tt}
            write(case/'receipt.json',row);rows.append(row);local.append(row)
            print(name,r,receipt['tail_interval_width_mHa'],flush=True)
        passing=[r for r in local if F(r['tail']['tail_interval_width_Ha'])<=F(16,10000)]
        best=min(passing,key=lambda r:r['tail']['factors']);r=best['requested_rank']
        witness=propose_rank(data,p,best['tail']['factors']);rank_receipt=rank_replay(data,witness)
        write(directory/'rank_obstruction.json',witness);write(directory/'rank_receipt.json',rank_receipt)
        support_rows=[]
        if s<=6:
            for retained in sorted({0,s,r,d}):
                tt=time.monotonic();case=directory/f'rank_{retained}';tail=json.loads((case/'tail.json').read_text())
                cert,cost=propose_support(data,p,tail);receipt=support_replay(data,tail,cert)
                size=write(case/'support.json',cert)
                row={'rank':retained,'support':receipt,'discovery':cost,'certificate_bytes':size,'case_seconds':time.monotonic()-tt}
                write(case/'support_receipt.json',row);support_rows.append(row)
        # Number centering is an exact elimination, checked against an uncentered control.
        pp=extract(data,False);tail,cost=propose_tail(data,pp,r);receipt=tail_replay(data,tail)
        write(directory/'uncentered_tail.json',tail);write(directory/'uncentered_receipt.json',receipt)
        system={'name':name,'fixture_source':str(source.relative_to(ROOT)),'hamiltonian_terms':len(p['h']),
                'spatial_orbitals':s,'full_coefficient_dimension':d,'initial_extraction_seconds':p['extract_seconds'],
                'best_rank':r,'best_tail':best,'rank_obstruction':rank_receipt,'supporting_bounds':support_rows,
                'uncentered_same_rank_tail_mHa':receipt['tail_interval_width_mHa'],
                'system_wall_seconds':time.monotonic()-t}
        systems.append(system);write(directory/'system.json',system)
    result={'systems':systems,'tail_cases':len(rows),'wall_seconds':time.monotonic()-start,
            'scope':'Frozen finite-basis rational electronic H; no physical-sector enumeration in campaign. Numerical source fixture/integral generation and separate controls excluded.'}
    write(out/'summary.json',result);return result


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args();run(args.out)
