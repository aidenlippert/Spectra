"""Generate and exactly verify a connected noncommuting frustration-free chain.

Two fermionic modes (up, down) live at each site.  For each bond,

  B_i = a_down,i a_up,i+1 - a_up,i a_down,i+1

and H=sum B_i^dagger B_i.  The all-up determinant at N=L is annihilated by
every B_i, giving a matched lower/upper energy zero.  Adjacent B_i overlap,
so this is a noncommuting control family rather than disjoint blocks.
"""
import argparse, json
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from experiments.marginal_symbolic import encode, verify, word_product
from experiments.marginal_determinant_tree import DeterminantOracle


def dagger(word):
    return tuple((1-c, i) for c, i in reversed(word))


def bond_words(site, sites):
    u, d = 2 * site, 2 * site + 1
    v, e = 2 * ((site + 1) % sites), 2 * ((site + 1) % sites) + 1
    return [(((0, d), (0, v)), 1), (((0, u), (0, e)), -1)]


def build(sites, hinput=None):
    factors=[]; h={}
    for i in range(sites):
        terms=bond_words(i, sites); words=[w for w,c in terms]
        for left, a in terms:
            for right, b in terms:
                for word, sign in word_product(dagger(left), right):
                    h[word]=h.get(word,0)+a*b*sign
        factors.append({"name":f"bond-{i}","words":[list(map(list,w)) for w in words],"factor":[[1,-1]]})
    h={w:c for w,c in h.items() if c}
    if hinput is None: hinput=h
    recognize(hinput, sites, sites)
    cert={"modes":2*sites,"particles":sites,"hamiltonian":encode(hinput),
          "number_multiplier":[],"b":"0","denominator":1,"blocks":factors}
    receipt=verify(cert)
    receipt.update({"sites":sites,"factors":sites,"factor_storage":2*sites,
                    "matched_upper":0,"witness":"all-up determinant",
                    "connected_overlap":True,"source_factor_used":False})
    return cert,receipt


def make_hamiltonian(sites):
    """Public H-only source; callers may mutate its returned polynomial."""
    h={}
    for i in range(sites):
        for left, a in bond_words(i, sites):
            for right, b in bond_words(i, sites):
                for word, sign in word_product(dagger(left), right):
                    h[word]=h.get(word,0)+a*b*sign
    return {w:c for w,c in h.items() if c}


def recognize(h, sites, particles):
    """Recognize only the exact H,M,N pattern; reject a coefficient mutation."""
    expected=make_hamiltonian(sites)
    if particles != sites or h != expected:
        raise ValueError("H does not match the recognized overlapping-bond family")
    return True


def determinant_upper(sites):
    """Exact occupation oracle: all-up determinant is killed by every B_i."""
    occupied={2*i for i in range(sites)}
    for i in range(sites):
        # Every term in B_i contains a down annihilator, absent from the state.
        if any(mode in occupied for _, mode in bond_words(i, sites)[0][0] if mode % 2):
            return False
    return 0


def main():
    p=argparse.ArgumentParser(); p.add_argument("--sites",type=int,default=4); p.add_argument("--outputdir",type=Path,default=Path("results/certificate_scaling/structural_chain")); a=p.parse_args()
    h=make_hamiltonian(a.sites)
    try:
        recognize(h,a.sites,a.sites)
        mutated=dict(h); key=next(iter(mutated)); mutated[key]+=1
        try: recognize(mutated,a.sites,a.sites); mutation_rejected=False
        except ValueError: mutation_rejected=True
        if not mutation_rejected: raise RuntimeError("recognizer mutation test failed")
        def energy(i):
            out={}
            for left,a0 in bond_words(i,a.sites):
                for right,b0 in bond_words(i,a.sites):
                    for w,s in word_product(dagger(left),right): out[w]=out.get(w,0)+a0*b0*s
            return {w:c for w,c in out.items() if c}
        e0,e1=energy(0),energy(1); comm={}
        for left,a0 in e0.items():
            for right,b0 in e1.items():
                for word,sign in word_product(left,right): comm[word]=comm.get(word,0)+a0*b0*sign
                for word,sign in word_product(right,left): comm[word]=comm.get(word,0)-b0*a0*sign
        comm={w:c for w,c in comm.items() if c}
        cert,receipt=build(a.sites,h)
        upper=DeterminantOracle(cert).upper({"states":[sum(1<<(2*i) for i in range(a.sites))],"amplitudes":[1]})
        if upper != 0: raise RuntimeError("determinant upper oracle was not zero")
        receipt.update({"mutation_rejected":mutation_rejected,"commutator_terms":len(comm),"determinant_upper":str(upper)})
    except Exception as exc:
        raise SystemExit(str(exc))
    a.outputdir.mkdir(parents=True,exist_ok=True)
    (a.outputdir/f"L{a.sites}.json").write_text(json.dumps(cert)+'\n')
    (a.outputdir/f"L{a.sites}_receipt.json").write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2))
if __name__=="__main__": main()
