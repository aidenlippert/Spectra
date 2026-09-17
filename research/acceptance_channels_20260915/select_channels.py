"""Bounded four-term channel proposals from an approximate molecular dual.

Only already supplied degree <=4 moments are queried. The approximate dual is
not an accepted obstruction. Added channels must earn a new exact energy bound.
"""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import time
from itertools import combinations
import numpy as np
from scipy import sparse
from experiments.marginal_symbolic import add, adj, canonical, encode, mono, product, word_product
from research.interacting_scaling_20260915.dictionary import representative


def channel_polynomial(words, integers, denominator):
    B = add(*(mono(tuple(map(tuple, w)), Fraction(v, denominator)) for w, v in zip(words, integers)))
    P = add(product(adj(B), B), product(B, adj(B)))
    if max(map(len, P), default=0) > 4:
        raise ValueError('Channel requires unavailable higher moments')
    return canonical(B), P


def t2_terms(p, q, r, s, t, u):
    """Quartic normal-order identity for {a†p bq br, b†u b†t as}."""
    terms = []
    if p == s: terms.append((((1,2*u+1),(1,2*t+1),(0,2*q+1),(0,2*r+1)),1))
    constant = int(r==u and q==t)-int(q==u and r==t)
    if constant: terms.append((((1,2*p),(0,2*s)),constant))
    for condition, creator, annihilator, sign in ((r==u,t,q,-1),(q==u,t,r,1),(r==t,u,q,1),(q==t,u,r,-1)):
        if condition:
            terms.append((((1,2*p),(1,2*creator+1),(0,2*annihilator+1),(0,2*s)),sign))
    return terms


def run(case, checkpoint, output, maximum=8, mode='spectator', max_terms=4, previous=None):
    start = time.monotonic()
    prepared = case/'prepared'
    frame = json.loads((prepared/'frame.json').read_text())
    T = sparse.load_npz(prepared/'twirl.npz')
    selected = np.load(prepared/'selected.npy')
    y = np.load(checkpoint)['y']
    values = -T[selected].T@(np.load(prepared/'scale.npy')*y)/np.load(prepared/'weights.npy')
    moments = {tuple(map(tuple, w)): float(v) for w, v in zip(frame['rows'], values) if len(w) <= 4}
    original_normalization = moments[()]
    if not np.isfinite(original_normalization) or original_normalization <= 0:
        raise ValueError('The proposal has no positive finite normalization')
    # This rescales an explicitly approximate proposal, not an exact witness.
    moments = {w:v/original_normalization for w,v in moments.items()}
    def expectation(word):
        total = 0.
        for normal, sign in word_product((), word):
            key = representative(normal)
            if key not in moments:
                raise ValueError(('No supplied moment', key))
            total += sign*moments[key]
        return total
    spaces = []
    for b in frame['blocks']:
        if b.get('kind') == 'spin_three_half_highest':
            g = frame['groups'][b['physical_group']]
            if g.get('signature') == [-1, 3]:
                space = set()
                for w in g['words']: space.update(canonical(mono(tuple(map(tuple, w)))))
                spaces.append(space)
    s = frame['modes']//2
    candidates = {}
    diagnostics = []
    for spectator in (range(s) if mode == 'spectator' else [None]):
        if spectator is None:
            triples = [(p,q,r) for p in range(s) for q,r in combinations(range(s),2)]
            if len(triples)>1024: raise ValueError('Global diagnostic operator-matrix envelope exceeded')
            words = [((1,2*p),(0,2*q+1),(0,2*r+1)) for p,q,r in triples]
            matrix = np.empty((len(words),len(words)))
            for i,(p,q,r) in enumerate(triples):
                for k in range(i,len(words)):
                    a,b,c = triples[k]
                    value = sum(sign*expectation(w) for w,sign in t2_terms(p,q,r,a,b,c))
                    matrix[i,k] = matrix[k,i] = value
        else:
            j = 2*spectator+1
            pairs = [(p, q) for p in range(s) for q in range(s) if q != spectator]
            words = [((1, 2*p), (0, 2*q+1), (0, j)) for p, q in pairs]
            matrix = np.empty((len(pairs), len(pairs)))
            for i, (p, q) in enumerate(pairs):
                for k in range(i, len(pairs)):
                    r, t = pairs[k]
                    value = expectation(((1,2*p),(0,2*q+1),(1,2*t+1),(0,2*r)))
                    if p == r:
                        value += expectation(((1,j),(0,j),(1,2*q+1),(0,2*t+1)))
                    if q == t:
                        value -= expectation(((1,j),(0,j),(1,2*r),(0,2*p)))
                    matrix[i,k] = matrix[k,i] = value
        ev, vectors = np.linalg.eigh(matrix)
        diagnostics.append({'spectator_spatial': spectator, 'matrix_dimension': len(words),
            'minimum_unverified_dual_eigenvalue': float(ev[0])})
        if mode == 'dense':
            inherited=np.zeros((len(words),0),dtype=np.int64)
            parent_hash=None
            if previous:
                parent=json.loads(previous.read_text())
                old=parent['channels'][0]
                if old['words']!=json.loads(json.dumps(words)) or old['denominator']!=10**7:
                    raise ValueError('Parent collective dictionary differs')
                if parent['frame_sha256']!=hashlib.sha256((prepared/'frame.json').read_bytes()).hexdigest():
                    raise ValueError('Parent base family differs')
                inherited=np.asarray(old['basis_integers'],dtype=np.int64)
                if inherited.shape[1]>8:raise ValueError('Nested selector is capped at twelve directions')
                U=np.linalg.qr(inherited/10**7,mode='reduced')[0]
                reduced=matrix-U@(U.T@matrix)-(matrix@U)@U.T+U@(U.T@matrix@U)@U.T
                ev,vectors=np.linalg.eigh((reduced+reduced.T)/2)
                diagnostics[-1]['minimum_complement_eigenvalue']=float(ev[0])
                parent_hash=hashlib.sha256(previous.read_bytes()).hexdigest()
            rank = min(4, int(np.count_nonzero(ev < -1e-7)))
            if rank:
                integers = np.column_stack((inherited,np.rint(vectors[:,:rank]*10**7).astype(np.int64)))
                candidates['dense'] = {'words':words,'basis_integers':integers.tolist(),
                    'denominator':10**7,'dimension':integers.shape[1],'operator_terms':len(words),
                    'parent_dimension':inherited.shape[1],'parent_file_sha256':parent_hash,
                    'unverified_dual_score':float(ev[0]),'unverified_eigenvalues':ev[:rank].tolist(),
                    'not_contained_in_a_retained_highest_weight_block':True}
            print(json.dumps(diagnostics[-1]),flush=True)
            continue
        for vector in vectors[:, :min(12, len(words))].T:
            for size in [v for v in (4,8,16,32,64) if v<=max_terms]:
                indices = np.argsort(abs(vector))[-size:]
                sub = matrix[np.ix_(indices, indices)]
                small_ev, small_vec = np.linalg.eigh(sub)
                if small_ev[0] >= -1e-7:
                    continue
                coeff = small_vec[:, 0]
                if coeff[np.argmax(abs(coeff))] < 0: coeff = -coeff
                integers = np.rint(coeff*10**7).astype(np.int64).tolist()
                selected_words = [words[i] for i in indices]
                B, P = channel_polynomial(selected_words, integers, 10**7)
                if any(set(B) <= space for space in spaces):
                    continue
                score = sum(float(c)*moments[representative(w)] for w,c in P.items())
                predicted = float((np.array(integers)/10**7)@sub@(np.array(integers)/10**7))
                if abs(score-predicted) > 1e-8:
                    raise ValueError('Structured channel matrix failed literal CAR comparison')
                key = json.dumps(sorted(B))
                candidate = {'words': selected_words, 'integers': integers, 'denominator': 10**7,
                    'operator_terms':size,'unverified_dual_score': score, 'P_terms': len(P),
                    'spatial_support': sorted({i//2 for w in B for _,i in w}),
                    'not_contained_in_a_retained_highest_weight_block': True}
                if key not in candidates or score < candidates[key]['unverified_dual_score']:
                    candidates[key] = candidate
                break
        print(json.dumps(diagnostics[-1]), flush=True)
    chosen = sorted(candidates.values(), key=lambda c:c['unverified_dual_score'])[:maximum]
    result = {'kind':'bounded_spectator_channel_proposals', 'mode':mode, 'maximum_terms':max_terms,'source_case':str(case),
        'source_checkpoint':str(checkpoint), 'checkpoint_sha256':hashlib.sha256(checkpoint.read_bytes()).hexdigest(),
        'frame_sha256':hashlib.sha256((prepared/'frame.json').read_bytes()).hexdigest(),
        'seconds':time.monotonic()-start, 'matrices':diagnostics,
        'original_approximate_dual_normalization':original_normalization,
        'positive_normalization_rescaling_applied':True,
        'distinct_noncontained_candidates':len(candidates), 'channels':chosen,
        'source_dual_exactly_feasible':False, 'family_obstruction_claimed':False,
        'new_energy_improvement_claimed':False, 'unavailable_higher_moments_used':False,
        'complete_global_cubic_map_built':False}
    with output.open('x') as stream: json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps({'selected':len(chosen),'seconds':result['seconds']}),flush=True)


if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('case',type=Path);parser.add_argument('checkpoint',type=Path)
    parser.add_argument('output',type=Path);parser.add_argument('--maximum',type=int,default=8)
    parser.add_argument('--mode',choices=('spectator','global','dense'),default='spectator')
    parser.add_argument('--max-terms',type=int,choices=(4,8,16,32,64),default=4)
    parser.add_argument('--previous',type=Path)
    args=parser.parse_args()
    if args.previous and args.mode!='dense':parser.error('--previous requires dense selection')
    run(args.case.resolve(),args.checkpoint.resolve(),args.output.resolve(),args.maximum,args.mode,args.max_terms,
        args.previous.resolve() if args.previous else None)
