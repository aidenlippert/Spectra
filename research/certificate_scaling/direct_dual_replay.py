"""Stdlib exact replay of a DD quadratic-cone obstruction, independent of LP."""
from fractions import Fraction as F
from itertools import combinations
import hashlib
import json
from pathlib import Path
import sys
import time
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from experiments.marginal_symbolic import add, decode, multiplier_basis, number_shift, product, verify, word_product

ROOT=Path(__file__).resolve().parents[2]

def replay(dual, source):
    start=time.monotonic();m,n=source['modes'],source['particles']
    y={tuple(map(tuple,row['word'])):F(row['value']) for row in dual['dual_words']}
    if len(y)!=len(dual['dual_words']):raise ValueError('Duplicate dual coordinates')
    if y.get(())!=1 or max(map(abs,y.values()))>1:raise ValueError('Dual normalization failed')
    def pairing(poly):return sum((c*y.get(w,F(0)) for w,c in poly.items()),F(0))
    shift=number_shift(m,n);ideals=multiplier_basis(m,max_body=2)
    for x in ideals:
        if pairing(product(shift,x)):raise ValueError('Dual ideal constraint failed')
    dagger=lambda w:tuple((1-c,i) for c,i in reversed(w))
    linear=[((0,i),) for i in range(m)]
    pairs=[((0,j),(0,i)) for i,j in combinations(range(m),2)]
    words=linear+[dagger(w) for w in linear]+pairs+[dagger(w) for w in pairs]+[((1,i),(0,j)) for i in range(m) for j in range(m)]
    diagonal=[pairing(dict(word_product(dagger(w),w))) for w in words]
    if min(diagonal)<0:raise ValueError('Negative monomial square')
    count=len(words);negative_minor=None
    for i,u in enumerate(words):
        for j in range(i+1,len(words)):
            v=words[j]
            cross=pairing(add(dict(word_product(dagger(u),v)),dict(word_product(dagger(v),u))))
            if diagonal[i]+diagonal[j]-abs(cross)<0:raise ValueError('Signed pair square constraint failed')
            count+=2
            minor=diagonal[i]*diagonal[j]-(cross/2)**2
            same_charge=sum(2*c-1 for c,_ in u)==sum(2*c-1 for c,_ in v)
            if minor<0 and same_charge and negative_minor is None:
                negative_minor={'u':u,'v':v,'diagonal_u':str(diagonal[i]),'diagonal_v':str(diagonal[j]),'cross':str(cross),'determinant':str(minor)}
    U=pairing(decode(source['hamiltonian'],m,4))
    L=F(verify(source)['lower'])
    return {'valid':True,'restricted_cone':'all quadratic word singleton and equal-magnitude signed pair squares, X body <= 2, full coefficient-l1 residual',
            'dual_ceiling':str(U),'dual_ceiling_float':float(U),'independent_ground_lower':str(L),
            'unavoidable_gap_at_least':str(L-U),'unavoidable_gap_at_least_float':float(L-U),
            'square_constraints_checked':count,'ideal_constraints_checked':len(ideals),
            'includes_cross_charge_atoms':True,'arbitrary_pair_PSD_counterexample':negative_minor,
            'seconds':time.monotonic()-start,'uses_external_numeric_library':False}

if __name__=='__main__':
    p=ROOT/'results/certificate_scaling/direct_dual_obstruction'
    d=p/'receipt.json';s=ROOT/'results/marginal_molecule_stress/h4_square_degree3_certificate.json'
    out=replay(json.loads(d.read_text()),json.loads(s.read_text()))
    out['dual_sha256']=hashlib.sha256(d.read_bytes()).hexdigest()
    out['source_sha256']=hashlib.sha256(s.read_bytes()).hexdigest()
    (p/'independent_replay.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ('dual_ceiling','independent_ground_lower','unavoidable_gap_at_least')},indent=2))
