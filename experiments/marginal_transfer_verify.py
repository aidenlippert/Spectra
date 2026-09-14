"""Dependency-free replay of connected transfer certificates and upper witnesses."""
from fractions import Fraction as F
import json
from pathlib import Path

from experiments.marginal_symbolic import decode, verify

def apply_word(word, state):
    amp=1; bits=state
    for creation,i in reversed(word):
        amp *= -1 if ((bits & ((1<<i)-1)).bit_count() & 1) else 1
        if creation:
            if bits>>i&1: return None
            bits |= 1<<i
        else:
            if not (bits>>i&1): return None
            bits &= ~(1<<i)
    return bits,amp

def replay(certificate):
    modes=certificate.get("modes"); particles=certificate.get("particles")
    if type(modes) is not int or type(particles) is not int or not 0<=particles<=modes:
        raise ValueError("invalid mode or particle count")
    states=[s for s in range(1<<modes) if s.bit_count()==particles]
    witness=certificate.get("independent_upper")
    if not isinstance(witness,dict): raise ValueError("missing independent upper witness")
    amps=witness.get("amplitudes")
    if type(amps) is not list or len(amps)!=len(states) or any(type(a) is not int for a in amps) or not any(amps):
        raise ValueError("upper amplitudes must be nonzero integers in ascending bitstring order")
    h=decode(certificate["hamiltonian"],modes,4)
    norm=sum(a*a for a in amps)
    energy=F(0)
    for word,coefficient in h.items():
        for source,a in zip(states,amps):
            result=apply_word(word,source)
            if result:
                destination,sign=result
                energy += coefficient*a*amps[states.index(destination)]*sign
    upper=F(energy,norm)
    lower=F(verify(certificate)["lower"])
    width=upper-lower
    if width<0: raise ValueError("inconsistent lower/upper interval")
    return {"lower":str(lower),"upper":str(upper),"norm":str(norm),"width":str(width),"width_float":float(width)}

if __name__=="__main__":
    import argparse
    p=argparse.ArgumentParser(); p.add_argument("certificate",type=Path); a=p.parse_args()
    print(json.dumps(replay(json.loads(a.certificate.read_text())),indent=2))
