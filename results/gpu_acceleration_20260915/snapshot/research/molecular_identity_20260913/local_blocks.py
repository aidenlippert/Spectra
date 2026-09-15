"""Bounded orbital-local cubic positivity, using the existing exact exporter."""
from itertools import combinations
from pathlib import Path
import argparse
import hashlib
import json
import time

from experiments.marginal_symbolic import canonical, mono, decode
from experiments.marginal_hunt_car import adj
from research.molecular_identity_20260913.scan import cubic_words


def local_groups(supports, signature):
    groups=[];seen=set()
    for support in supports:
        if len(set(support))!=len(support) or list(support)!=sorted(support):
            raise ValueError('Distinct ascending orbital support required')
        for charge in (-1,-3):
            minus=[mono(w) for w in cubic_words(support,charge)]
            for sign,polys in ((-1,minus),(1,[canonical(adj(p)) for p in minus])):
                parts={}
                for p in polys:
                    key=signature(next(iter(p)))
                    parts.setdefault(key,[]).append(p)
                for key,part in sorted(parts.items()):
                    # Scalar linear constraints already occur in the baseline.
                    if len(part)==1 and max(map(len,part[0]))==1:continue
                    frozen=tuple(sorted(tuple(sorted(p.items())) for p in part))
                    if frozen in seen:continue
                    seen.add(frozen)
                    groups.append({'name':f'local{tuple(support)}:{sign*abs(charge)}:{key}',
                                   'polynomials':part})
    return groups


def run(fixture,out,width=3,solver_seconds=60):
    from research.certificate_scaling.adaptive_block_discovery import partition
    from research.certificate_scaling.commutator_dictionary import run as solve
    started=time.monotonic();raw=fixture.read_bytes();f=json.loads(raw)
    m=f['modes'];h=decode(f['hamiltonian'],m,4)
    if width not in (3,4):raise ValueError('Supported local widths are three and four')
    _,signature,_=partition(h,m,'quadratic',True)
    supports=list(combinations(range(m),width))
    groups=local_groups(supports,signature)
    selected=time.monotonic()
    receipt=solve(h,m,f['particles'],out,solver_seconds=solver_seconds,
                  creator_channels=True,map_backend='contraction',additional_groups=groups)
    result={'policy':'all_orbital_subsets_cubic_symmetry_blocks','width':width,
            'support_count':len(supports),'additional_block_count':len(groups),
            'selection_seconds':selected-started,'total_wall_seconds':time.monotonic()-started,
            'fixture_sha256':hashlib.sha256(raw).hexdigest(),
            'many_body_space_enumerated':False,'source_dual_used':False,
            'source_factors_used':False,'source_upper_used':False,
            'solver_receipt':receipt,'scope':'Finite local positivity restriction; no local particle number is assumed.'}
    (out/'experiment.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--fixture',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True);p.add_argument('--width',type=int,default=3)
    p.add_argument('--solver-seconds',type=float,default=60)
    a=p.parse_args();run(a.fixture,a.out,a.width,a.solver_seconds)
