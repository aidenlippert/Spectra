"""Replayable mathematical examples, not experiments on real materials."""
from fractions import Fraction as F
from pathlib import Path
from dataclasses import asdict
from hashlib import sha256
from time import perf_counter
import json
from .mission_policy_cover import construct,certify,Sample
ROOT=Path(__file__).resolve().parents[1]


def reaction(policy):
    x=F(0);history=[x]
    for u in policy:
        x=x+u*(1-x)/2-x/4;history.append(x)
    return (sum(policy),F(9,20)-x),history


def thermal(policy):
    hot=cold=F(0);history=[(hot,cold)]
    for u in policy:
        hot,cold=hot/2+cold/4+u/2,hot/4+cold/2;history.append((hot,cold))
    return (sum(policy),F(1,10)-cold,hot-F(3,10)),history


def pack(v):
    if isinstance(v,F):return str(v)
    if isinstance(v,Sample):return pack(asdict(v))
    if isinstance(v,dict):return {k:pack(x) for k,x in v.items()}
    if isinstance(v,(tuple,list)):return [pack(x) for x in v]
    return v


def run():
    rows=[];epsilon=F(1,5)
    for name,physical,L,optimum in [('two_state_reaction',reaction,[2,F(7,8)],F(9,10)),
                                  ('two_node_thermal',thermal,[2,F(1,8),F(3,4)],F(4,5))]:
        def oracle(point,requested):
            values,_=physical(point)
            return tuple((v,v) for v in values)
        start=perf_counter();ans=construct(oracle,L,2,epsilon,72);constructed=perf_counter()
        receipt=ans['receipt']
        replay=certify(ans['paths'],ans['samples'],L,2,epsilon);checked=perf_counter()
        if receipt!=replay or receipt['status']!='conditional_certificate':raise AssertionError('replay')
        values,trajectory=physical(receipt['policy'])
        if any(g>0 for g in values[1:]) or not receipt['lower']<=optimum<=receipt['upper']:
            raise AssertionError('independent physical-equation/optimum check')
        rows.append(dict(name=name,dimension=2,lipschitz=L,epsilon=epsilon,known_optimum=optimum,
            result=ans,selected_trajectory=trajectory,selected_values=values,
            construction_seconds=constructed-start,rechecking_seconds=checked-constructed,
            evidence='exact supplied process equations; no physical measurements'))
    files=['experiments/mission_policy_cover.py','experiments/mission_examples.py',
           'tests/test_mission_policy_cover.py','research/MISSION.md','research/mission/CONSTRUCTIVE_BRIDGE.md']
    report=dict(schema='mission-math-replay-1',scope='mathematical verification only',rows=rows,
        source_sha256={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in files},
        physical_validation=False,autonomous_discovery=False,universal_synthesis=False)
    destination=ROOT/'results/mission';destination.mkdir(exist_ok=True)
    (destination/'mathematical_examples.json').write_text(json.dumps(pack(report),indent=2)+'\n')
    print(json.dumps(pack([dict(name=r['name'],policy=r['result']['receipt']['policy'],
        lower=r['result']['receipt']['lower'],upper=r['result']['receipt']['upper'],
        observations=r['result']['observation_calls'],cells=r['result']['leaf_cells']) for r in rows]),indent=2))


if __name__=='__main__':run()
