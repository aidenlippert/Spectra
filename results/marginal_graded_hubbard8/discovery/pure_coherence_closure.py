"""Exact full spin-word stationary law and classical Markov-extension witness."""
from pathlib import Path
from fractions import Fraction as F
from math import lcm
import argparse,hashlib,json,sys
from full_overlap_density import ROOT,orbit
sys.set_int_max_str_digits(10000)

def main():
    parser=argparse.ArgumentParser();parser.add_argument('directory',type=Path)
    directory=parser.parse_args().directory.resolve();cp=directory/'range_two_family_limit_certificate.json';rp=directory/'range_two_family_limit_replay.json'
    c=json.loads(cp.read_text());r=json.loads(rp.read_text())
    if not r.get('accepted') or not r.get('full_spin_word') or not r.get('pure_coherence'):raise ValueError('Accepted full spin-word family required')
    for name,h in r['source_sha256'].items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=h:raise ValueError('Stale family source')
    sources=[]
    for item in c['mixture']:
        vector={int(s):a for s,a in item['vector'].items() if a};norm=sum(a*a for a in vector.values())
        sources.append((F(item['weight'])/(8*norm),vector))
    denominator=lcm(*(weight.denominator for weight,vector in sources))
    if denominator.bit_length()>12000:raise ValueError('Common denominator exceeds12000-bit budget')
    law=[0]*4096
    for weight,vector in sources:
        integer=weight*denominator
        if integer.denominator!=1:raise ValueError('Noninteger common-denominator multiplier')
        for image,multiplicity in orbit(vector).items():
            for state,a in image:law[state]+=int(integer)*multiplicity*a*a
    if min(law)<0 or sum(law)!=denominator:raise ValueError('Invalid positive normalized six-word law')
    left=[0]*1024;right=[0]*1024
    for state,mass in enumerate(law):left[state&1023]+=mass;right[state>>2]+=mass
    if left!=right:raise ValueError('Full spin-resolved five-word overlap did not close')
    if len(r['family_replay']['spin_word_moments'])!=120 or any(map(F,r['family_replay']['spin_word_moments'].values())):raise ValueError('Incomplete spin-word basis closure')
    if len(r['family_replay']['pure_coherence_moments'])!=2 or any(map(F,r['family_replay']['pure_coherence_moments'].values())):raise ValueError('Incomplete pure-coherence closure')
    # Each transition carries stationary flow mass/D; cancellation of its row denominator
    # proves both row normalization and stationary inflow with exact integer sums above.
    active=sum(bool(v) for v in left);edges=sum(bool(v) for v in law)
    for state,mass in enumerate(law):
        if mass and (not left[state&1023] or not left[state>>2]):raise ValueError('Positive flow enters a zero-mass word')
    digest=lambda values:hashlib.sha256(''.join(str(v)+'\n' for v in values).encode()).hexdigest()
    files={Path(__file__).resolve(),cp,rp,Path(sys.modules['full_overlap_density'].__file__).resolve()}
    for module in tuple(sys.modules.values()):
        name=getattr(module,'__file__',None)
        if name and str(Path(name).resolve()).startswith(str(ROOT/'experiments')+'/'):files.add(Path(name).resolve())
    result={'accepted':True,'full_spin_word_overlap_exactly_zero':True,'spin_word_zero_moments':120,'pure_coherence_zero_moments':2,
            'alphabet':['empty','up','down','double'],'mixture_sources':len(c['mixture']),
            'common_denominator_bits':denominator.bit_length(),'common_denominator':str(denominator),
            'six_site_diagonal_integer_sha256':digest(law),'five_site_stationary_integer_sha256':digest(left),
            'markov_states_with_positive_mass':active,'markov_edges_with_positive_flow':edges,
            'classical_markov_extension':True,
            'construction':'For positive five-word prefix mass q(u), append a with probability p6(u,a)/q(u). Prefix and suffix sums coincide exactly, proving normalized rows and stationary inflow q. For zero-mass prefixes append empty deterministically; zero stationary mass makes these rows irrelevant. This gives an order-five stationary classical process reproducing the full six-site spin-resolved diagonal law.',
            'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(files)},
            'scope':'All spin-resolved diagonal overlaps close for the positive symmetry-averaged local mixture, and its diagonal occupation law has a stationary classical extension. This does not extend its offdiagonal quantum coherences, impose a fixed global particle number, or establish general quantum representability.'}
    (directory/'pure_coherence_closure.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'accepted':True,'full_diagonal_overlap_zero':True,'classical_markov_extension':True,'active_states':active,'positive_edges':edges}))

if __name__=='__main__':main()
