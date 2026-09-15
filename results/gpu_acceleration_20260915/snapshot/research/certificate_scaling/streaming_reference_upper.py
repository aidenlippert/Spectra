"""Bounded exact upper replay for larger validation-only determinant witnesses.

No action cache or K-by-K retained matrix; O(K * Hamiltonian terms) arithmetic.
The existing small-witness oracle and its limits are unchanged.
"""
from fractions import Fraction as F
from math import lcm
from pathlib import Path
import sys,time
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from experiments.marginal_determinant_tree import DeterminantOracle


def compile_term(word,coefficient):
    required=occupied=flip=parity=0
    for creation,mode in reversed(word):
        bit=1<<mode;need=(1-creation)^bool(flip&bit)
        if required&bit and bool(occupied&bit)!=need:return None
        required|=bit
        if need:occupied|=bit
        if (flip&(bit-1)).bit_count()%2:coefficient=-coefficient
        parity^=bit-1;flip^=bit
    return required,occupied,flip,parity,coefficient


def upper(certificate,witness):
    start=time.monotonic();oracle=DeterminantOracle(certificate)
    if type(witness) is not dict:raise ValueError('Integer validation witness required')
    states,amps=witness.get('states'),witness.get('amplitudes')
    if (type(states) is not list or not 1<=len(states)<=65536 or type(amps) is not list or len(amps)!=len(states)
        or len(set(states))!=len(states) or any(not oracle.valid_state(s) for s in states)
        or any(type(a) is not int for a in amps) or not any(amps)):
        raise ValueError('Invalid bounded integer validation witness')
    if len(oracle.h)>20000:raise ValueError('Validation Hamiltonian term budget exceeded')
    denominator=lcm(*(c.denominator for c in oracle.h.values()))
    terms=[compile_term(w,int(c*denominator)) for w,c in oracle.h.items()]
    terms=[t for t in terms if t is not None]
    vector={s:a for s,a in zip(states,amps) if a};energy=0;contributing=0
    for state,x in vector.items():
        for required,occupied,flip,parity,c in terms:
            if state&required!=occupied:continue
            y=vector.get(state^flip,0)
            if y:
                sign=-1 if (state&parity).bit_count()%2 else 1
                energy+=sign*c*x*y;contributing+=1
    norm=sum(a*a for a in vector.values());value=F(energy,denominator*norm)
    return value,{'method':'bounded_streaming_integer_reference_upper_v1','witness_states':len(vector),
                  'hamiltonian_terms':len(terms),'word_state_checks':len(vector)*len(terms),
                  'contributing_word_state_pairs':contributing,'action_cache_states':len(oracle.cache),
                  'replay_seconds':time.monotonic()-start,'upper':str(value),
                  'scope':'validation upper witness only; exponentially growing FCI discovery cost is not removed'}
