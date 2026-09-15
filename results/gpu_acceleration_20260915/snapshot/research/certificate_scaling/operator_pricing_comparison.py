"""Fair small CAR-SOS comparison with exact replay where export succeeds."""
import json, time
from pathlib import Path
import numpy as np
from fractions import Fraction as F
from experiments.marginal_coefficient import dictionaries, solve_coefficients, export, hopping_model, symmetric_upper
from research.certificate_scaling.operator_pricing_sos import solve_one


def run(modes=4, seed=812):
    h=hopping_model(modes,F(1,5),False); blocks=dictionaries(modes,'mixed'); rng=np.random.default_rng(seed)
    upper=F(symmetric_upper(modes,F(1,5))['upper']); out=Path(__file__).resolve().parents[2]/'results/certificate_scaling/operator_pricing_comparison'; out.mkdir(parents=True,exist_ok=True)
    orders={'full':list(range(len(blocks))), 'fixed':list(range(3)), 'random':list(rng.permutation(len(blocks))[:3])}
    # Greedy choices and complete scoring ledger.
    chosen=[]; ledger=[]
    for _ in range(3):
        scored=[]; started=time.perf_counter()
        for i,b in enumerate(blocks):
            if i in chosen: continue
            try:
                s=solve_one(h,modes,[blocks[j] for j in chosen+[i]])
                ledger.append({'candidate':i,'round':len(chosen)+1,'status':s['solver_status'],'b':s['b']})
                if s['b'] is not None: scored.append((s['b'],i))
            except (RuntimeError,ValueError) as e: ledger.append({'candidate':i,'round':len(chosen)+1,'status':'error','error':str(e)})
        if not scored: break
        chosen.append(max(scored)[1]); ledger.append({'round':len(chosen),'scored':len(scored),'seconds':time.perf_counter()-started})
    orders['adaptive']=chosen
    receipt={'modes':modes,'dictionary_blocks':len(blocks),'upper_bound':str(upper),'upper_bound_float':float(upper),'strategies':{},'pricing_ledger':ledger}
    for name,order in orders.items():
        traces=[]
        budgets=[len(order)] if name=='full' else range(1,min(3,len(order))+1)
        for budget in budgets:
            subset=[blocks[i] for i in order[:budget]]; started=time.perf_counter(); s=solve_one(h,modes,subset); rec={'budget':budget,'order':[blocks[i]['name'] for i in order[:budget]],'solve':s}
            if s['b'] is not None:
                try:
                    sol=solve_coefficients(h,modes,modes//2,subset); cert,exact=export(h,modes,modes//2,subset,sol)
                    rec.update({'exact_replay_accepted':True,'exact_replay':exact,
                                'gap_to_upper':float(upper)-float(F(exact['lower'])),
                                'certificate_bytes':len(json.dumps(cert).encode())})
                    (out/f'{name}_b{budget}_certificate.json').write_text(json.dumps(cert,indent=2)+'\n')
                except Exception as e: rec.update({'exact_replay_accepted':False,'exact_replay_error':str(e)})
            traces.append(rec)
        receipt['strategies'][name]={'traces':traces}
    (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n'); return receipt

if __name__=='__main__': print(json.dumps(run(),indent=2))
