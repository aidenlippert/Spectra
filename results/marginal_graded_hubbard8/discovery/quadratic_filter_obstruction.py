"""Exact counterexample to reusing linear-filter scalar closure at degree two."""
from pathlib import Path
from fractions import Fraction as F
import sys,json,hashlib
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT))
from experiments.marginal_local_hubbard_block import _actions
from experiments.marginal_transfer_verify import apply_word

def main():
    h=_actions(2,0,1)
    def move(matrix,v):
        result={}
        for s,a in v.items():
            for t,b in matrix[s].items():result[t]=result.get(t,0)+a*b
        return {s:a for s,a in result.items() if a}
    fourth={s:{s:1} for s in range(16)}
    for _ in range(4):fourth={s:move(h,v) for s,v in fourth.items()}
    # R(d)=Tr_B[rho_B h^4]; the dependence is affine in d.
    # Testing two distinct d values establishes every coefficient.
    for d in (F(0),F(1,2)):
        probabilities=[d,F(1,2)-d,F(1,2)-d,d]
        for s in range(4):
            C=int(s in (0,3))
            for t in range(4):
                actual=sum(probabilities[b]*fourth[s+4*b].get(t+4*b,0) for b in range(4))
                expected=(4-6*d+(12*d-3)*C)*int(s==t)
                if actual!=expected:raise ValueError('Fourth-moment partial-operator identity failed')
    A={3:1,12:1,9:1,6:-1};B={9:1,6:-1}
    state={a+(b<<4):x*y for a,x in A.items() for b,y in B.items()}
    terms=[]
    for spin in (0,1):
        a,b=2+spin,4+spin
        terms.extend([((1,a),(0,b)),((1,b),(0,a))])
    matrix={}
    for s in range(256):
        matrix[s]={}
        for word in terms:
            image=apply_word(word,s)
            if image:
                t,phase=image;matrix[s][t]=matrix[s].get(t,0)-phase
    twice=move(matrix,move(matrix,state))
    filtered={s:state.get(s,0)+twice.get(s,0) for s in state.keys()|twice.keys()}
    norm=lambda v:sum(a*a for a in v.values())
    doublon=lambda v:F(sum(a*a*int((s&3)==3) for s,a in v.items()),norm(v))
    if norm(state)!=8 or norm(filtered)!=44 or doublon(state)!=F(1,4) or doublon(filtered)!=F(2,11):
        raise ValueError('Remote doublon counterexample failed')
    # Both contacts initially have each spin occupation one half.
    for site in (1,2):
        for spin in (0,1):
            if F(sum(a*a*((s>>(2*site+spin))&1) for s,a in state.items()),norm(state))!=F(1,2):
                raise ValueError('Half occupation hypothesis failed')
    result={'accepted':True,'fourth_moment_partial_operator':'(4-6*d)I+(12*d-3)C_A',
        'C_A':'1-n_A_up-n_A_down+2*D_A','d':'contact doublon probability of B',
        'scalar_for_arbitrary_A_iff_d':'1/4','filter':'I+h^2',
        'A_sparse_amplitudes':A,'B_sparse_amplitudes':B,'before_norm':norm(state),'after_norm':norm(filtered),
        'remote_doublon_before':str(doublon(state)),'remote_doublon_after':str(doublon(filtered)),
        'scope':'Exact failure of general remote-observable scalar closure for quadratic hopping filters, even with fixed block spin populations and half contact occupations. Does not rule out finite boundary transfer contraction.'}
    paths=['results/marginal_graded_hubbard8/discovery/quadratic_filter_obstruction.py',
           'experiments/marginal_local_hubbard_block.py','experiments/marginal_transfer_verify.py']
    result['source_sha256']={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}
    dest=ROOT/'results/marginal_graded_hubbard8/adapted_block/quadratic_filter_obstruction.json'
    dest.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)

if __name__=='__main__':main()
