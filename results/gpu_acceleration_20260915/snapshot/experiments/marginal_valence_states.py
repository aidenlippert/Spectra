"""Compact valence (one-electron-per-site) states and singlet witness."""
from itertools import combinations
from math import comb

def valence_states(modes, particles, cap=32):
    if type(modes) is not int or type(particles) is not int or not 4<=modes<=64 or modes%4 or particles!=modes//2:
        raise ValueError('Closed half-filled even-site sector required')
    m=modes//2
    count=comb(m,particles//2)
    if type(cap) is not int or cap<1 or cap>32 or count>cap: raise ValueError('Valence space exceeds cap')
    return sorted(sum(1<<(2*i) for i in up)|sum(1<<(2*i+1) for i in range(m) if i not in up) for up in combinations(range(m),particles//2))

def dimer_singlet_witness(modes, particles):
    allowed=set(valence_states(modes,particles))
    states={0:1}; m=modes//2
    for pair in range(0,m,2):
        nxt={}
        for s,a in states.items():
            for x,sgn in (((2*pair,2*(pair+1)+1),1),((2*pair+1,2*(pair+1)),-1)):
                z=s|sum(1<<q for q in x); nxt[z]=nxt.get(z,0)+a*sgn
        states=nxt
    if not set(states)<=allowed: raise ValueError('Singlet construction escaped valence sector')
    return {'states':sorted(states),'amplitudes':[states[s] for s in sorted(states)]}

def spin_raise_residual(witness,modes):
    from experiments.marginal_transfer_verify import apply_word
    out={}
    for s,a in zip(witness['states'],witness['amplitudes']):
        for i in range(modes//2):
            z=apply_word(((1,2*i),(0,2*i+1)),s)
            if z: t,sgn=z; out[t]=out.get(t,0)+a*sgn
    return {s:a for s,a in out.items() if a}
