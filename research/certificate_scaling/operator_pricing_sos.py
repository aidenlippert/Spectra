"""Actual CAR SOS dictionary pricing on small Hubbard Hamiltonians.

Each candidate is a positive Gram block from ``marginal_coefficient``.  A
restricted SDP is solved for matched fixed, random, and greedy block sets.
The residual is reconstructed independently from coefficient maps; it is
reported as a numerical diagnostic and is not promoted to an exact proof.
"""
from fractions import Fraction as F
import argparse, json, time
from pathlib import Path
import numpy as np

from experiments.marginal_coefficient import (dictionaries, coefficient_rows,
    gram_map, solve_coefficients, hopping_model, export)
from experiments.marginal_symbolic import multiplier_basis, number_shift, product


def residual(h, modes, particles, blocks, solution):
    rows = coefficient_rows(modes); lookup = {w:i for i,w in enumerate(rows)}
    out = np.zeros(len(rows)); out[lookup[()]] = solution['b']
    basis = solution['basis']; shift = number_shift(modes, particles)
    for q, val in zip(basis, solution['x']):
        for w, c in product(shift, q).items():
            if w in lookup: out[lookup[w]] += float(c) * float(val)
    for block, gram in zip(blocks, solution['grams']):
        g = gram_map(block['words'], lookup)
        out += np.asarray(g @ np.asarray(gram).reshape(-1, order='C')).ravel()
    rhs = np.array([float(h.get(w, 0)) for w in rows])
    return float(np.max(np.abs(out-rhs))), float(np.linalg.norm(out-rhs))


def solve_one(h, modes, blocks):
    started = time.perf_counter()
    try:
        sol = solve_coefficients(h, modes, modes//2, blocks)
    except RuntimeError as exc:
        return {'b': None, 'residual_max': None, 'residual_l2': None,
                'gram_rank': [], 'support_blocks': len(blocks),
                'wall_seconds': time.perf_counter()-started,
                'solver_status': str(exc), 'coefficient_rows': None}
    mx, l2 = residual(h, modes, modes//2, blocks, sol)
    return {'b': sol['b'], 'residual_max': mx, 'residual_l2': l2,
            'gram_rank': [int(np.linalg.matrix_rank((q+q.T)/2, tol=1e-8)) for q in sol['grams']],
            'support_blocks': len(blocks), 'wall_seconds': time.perf_counter()-started,
            'solver_status': sol['status'], 'coefficient_rows': sol['coefficient_rows']}


def one_case(modes, budget, seed=901, asymmetric=False):
    h = hopping_model(modes, F(1, 5), asymmetric=asymmetric)
    all_blocks = dictionaries(modes, 'mixed')
    if budget > len(all_blocks): raise ValueError('budget exceeds dictionary')
    rng = np.random.default_rng(seed); fixed = list(range(budget)); random = list(rng.permutation(len(all_blocks))[:budget])
    # Greedy pricing uses the restricted SDP objective as its score.
    chosen=[]; adaptive_trace=[]; pricing_log=[]
    for _ in range(budget):
        candidates=[]
        round_started = time.perf_counter()
        for i in range(len(all_blocks)):
            if i in chosen: continue
            trial=chosen+[i]
            try:
                scored=solve_one(h,modes,[all_blocks[j] for j in trial])
                if scored['b'] is not None: candidates.append((scored['b'], i))
                else: pricing_log.append({'candidate': i, 'status': scored['solver_status']})
            except (ValueError, RuntimeError) as exc:
                pricing_log.append({'candidate': i, 'status': 'error', 'error': str(exc)})
            except Exception:
                raise
        if not candidates: break
        _, pick=max(candidates); chosen.append(pick)
        adaptive_trace.append(solve_one(h,modes,[all_blocks[j] for j in chosen]))
        pricing_log.append({'round': len(chosen), 'candidate_count': len(all_blocks)-len(chosen)+1,
                            'round_wall_seconds': time.perf_counter()-round_started})
    receipt={'scope':'CAR SOS restricted SDP; residual is independently reconstructed numerical diagnostic',
             'modes':modes,'particles':modes//2,'dictionary_blocks':len(all_blocks),'budget':budget,
             'asymmetric':asymmetric,'strategies':{}}
    for name, order in [('adaptive',chosen),('fixed',fixed),('random',random)]:
        trace=adaptive_trace if name=='adaptive' else [solve_one(h,modes,[all_blocks[j] for j in order[:k]]) for k in range(1,len(order)+1)]
        entry={'order':[all_blocks[i]['name'] for i in order], 'trace':trace}
        if name == 'adaptive': entry['pricing_log'] = pricing_log
        if name == 'adaptive' and order:
            try:
                solution = solve_coefficients(h, modes, modes//2, [all_blocks[i] for i in order])
                certificate, exact = export(h, modes, modes//2, [all_blocks[i] for i in order], solution)
                entry['exact_replay'] = exact
                entry['exact_replay_accepted'] = True
                entry['certificate_bytes'] = len(json.dumps(certificate).encode())
                cert_out=(Path(__file__).resolve().parents[2]/'results/certificate_scaling/operator_pricing'/
                          f"m{modes}_b{budget}_{'asym' if asymmetric else 'sym'}_adaptive_certificate.json")
                cert_out.write_text(json.dumps(certificate, indent=2)+'\n')
            except Exception as exc:
                entry['exact_replay_accepted'] = False
                entry['exact_replay_error'] = str(exc)
        receipt['strategies'][name]=entry
    return receipt


def run(modes=4, budget=2, seed=901, campaign=False):
    if not campaign:
        receipt = one_case(modes, budget, seed)
        cases = [receipt]
    else:
        # Keep local campaign bounded; modes 8 is intentionally opt-in because
        # restricted SDP dimensions grow rapidly with the quartic row space.
        cases=[]
        for m in (4, 6, 8):
            for b in (1, 2, 3):
                for asym in (False, True):
                    if asym and m != 6: continue
                    if m == 8 and b > 2: continue
                    cases.append(one_case(m, b, seed + m + b, asym))
        receipt={'scope':'CAR SOS campaign; each proposal has independent residual and exact replay attempt', 'cases':cases}
    out=Path(__file__).resolve().parents[2]/'results/certificate_scaling/operator_pricing'; out.mkdir(parents=True,exist_ok=True)
    (out/('sos_campaign.json' if campaign else 'sos_receipt.json')).write_text(json.dumps(receipt,indent=2)+'\n')
    return receipt


if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('--modes',type=int,default=4); p.add_argument('--budget',type=int,default=2); p.add_argument('--campaign',action='store_true'); a=p.parse_args(); print(json.dumps(run(a.modes,a.budget,campaign=a.campaign),indent=2))
