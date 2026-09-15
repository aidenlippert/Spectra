"""Development probe: fixed-grid segmented sparse Taylor certificates.

This is a conventional time-stepping baseline.  Every segment is included in
one certificate, so endpoint jumps and all residual integrals are charged by
the unchanged v7 checker.
"""
from fractions import Fraction as F
from pathlib import Path
from time import perf_counter
import json
from .v7_headroom import model, TOL
from .v7_certificate import Generator, derive_certificate, check_certificate, evaluate, Piece
from .v7_reducers import full_taylor
from .v7_adaptive_taylor import adaptive_taylor

ROOT = Path(__file__).resolve().parents[1]

def attempt(n, gamma, horizon, segments, order, family="xxz", deadline=20):
    began = perf_counter(); h, initial = model(n, family)
    duration = F(horizon) / segments; current = dict(initial); pieces=[]
    work=[]
    try:
        gen = Generator(h, gamma, n, max_terms=512)
        for _ in range(segments):
            piece, w = full_taylor(gen, current, duration, order)
            pieces.append(piece); work.append(w)
            current = evaluate(piece.coefficients, duration)
            if perf_counter()-began > deadline: raise RuntimeError("calculation time budget")
        best = None; tried=[]
        for integration in ("power", "bernstein"):
            for grouping in ("l1", "firstfit", "weighted"):
                wg = Generator(h, gamma, n, max_terms=512)
                t0=perf_counter(); witness=derive_certificate(wg, initial, pieces, grouping, integration)
                wt=perf_counter()-t0
                cg = Generator(h, gamma, n, max_terms=512)
                t0=perf_counter(); result=check_certificate(cg, initial, pieces, witness, TOL, expected_time=F(horizon))
                ct=perf_counter()-t0
                row=dict(grouping=grouping,integration_basis=integration,status=result['status'],bound=result.get('bound'),
                         witness_seconds=wt,checking_seconds=ct,polynomial_entries=sum(len(c) for p in pieces for c in p.coefficients),
                         peak_support=max(len(c) for p in pieces for c in p.coefficients),witness_work=witness['construction_cost'],
                         checking_work=result.get('checking_cost'),generator_work=result.get('generator_cost'))
                tried.append(row)
                if result['status']=='certified': best=row; break
            if best: break
        return dict(status='certified' if best else 'over_tolerance', best=best, attempts=tried,
                    construction_work=work, elapsed=perf_counter()-began, segments=segments, order=order)
    except (ValueError, RuntimeError) as exc:
        return dict(status='refused', reason=str(exc), elapsed=perf_counter()-began,
                    segments=segments, order=order, construction_work=work,
                    generator_work=dict(gen.cost) if 'gen' in locals() else {})

def adaptive_baseline(n, gamma, horizon, family):
    began=perf_counter(); h, initial=model(n,family); gen=Generator(h,gamma,n,max_terms=512)
    try:
        piece,witness,proposal=adaptive_taylor(gen,initial,F(horizon),TOL,max_order=24)
        cg=Generator(h,gamma,n,max_terms=512); result=check_certificate(cg,initial,[piece],witness,TOL,expected_time=F(horizon))
        return dict(status=result['status'],bound=result.get('bound'),elapsed=perf_counter()-began,
                    order=len(piece.coefficients)-1,polynomial_entries=sum(map(len,piece.coefficients)),
                    construction_work=proposal,checking_work=result.get('checking_cost'),generator_work=result.get('generator_cost'))
    except (ValueError,RuntimeError) as exc:
        return dict(status='refused',reason=str(exc),elapsed=perf_counter()-began,generator_work=dict(gen.cost))

def run():
    rows=[]
    for n in (3,4,6):
      for family in ('xxz','mixed'):
       for gamma in (F(0),F(1,5),F(2)):
        for horizon in (F(1,5),F(1,2)):
         base=adaptive_baseline(n,gamma,horizon,family)
         for segments in (1,2,4):
          arm=[]; accepted=False
          for order in (4,8,12,16,24):
           row=attempt(n,gamma,horizon,segments,order,family); row.update(n=n,family=family,gamma=str(gamma),horizon=str(horizon),baseline=base); rows.append(row); arm.append(row)
           print(n,family,gamma,horizon,segments,order,row['status'],flush=True)
           if row['status']=='certified': accepted=True; break
           # A budget refusal at a higher order cannot be repaired by more order
           # for this fixed segmentation; retain the refusal and close the arm.
           if row['status']=='refused' and 'budget' in row.get('reason','').lower(): break
    out=dict(protocol='v7 fixed task; segmented Taylor development baseline', tolerance=str(TOL), rows=rows,
             claim='conventional time stepping diagnostic; no novelty or heldout evaluation')
    path=ROOT/'results/v8/time_reduction_probe.json'; path.write_text(json.dumps(out,indent=2)); return out
if __name__=='__main__': run()
