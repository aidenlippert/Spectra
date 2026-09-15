"""Direct rank-one SOS proposal from H, with exact verifier gate.

No source Gram, Fock-sector matrix, or full-rank SDP solution is used. The
nonconvex optimizer only proposes factors; `verify` decides the exported
rational certificate.
"""
from __future__ import annotations
import argparse, json, sys, time
from fractions import Fraction as F
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from experiments.marginal_collective import hopping_polynomial
from experiments.marginal_coefficient import coefficient_rows, dictionaries, gram_map
from experiments.marginal_symbolic import encode, verify, decode


def run(modes=4, t=F(1, 5), iterations=120, seed=7, denominator=10**6, fixture=None):
    if fixture:
        source = json.loads(Path(fixture).read_text())
        modes, particles = source["modes"], source["particles"]
        h = decode(source["hamiltonian"], modes, 4)
        # The fixture contributes only H,M,N.  Operator supports are generated
        # from the requested dictionary, never copied from an accepted factor.
        blocks = dictionaries(modes, "quadratic")
    else:
        particles = modes // 2; h = hopping_polynomial(modes, t); blocks = dictionaries(modes, "quadratic")
    rows = coefficient_rows(modes); lookup = {w:i for i,w in enumerate(rows)}
    maps = [gram_map(b["words"], lookup) for b in blocks]
    target = np.array([float(h.get(w, 0)) for w in rows])
    shapes = [len(b["words"]) for b in blocks]
    offsets = np.cumsum([0] + shapes)
    rng = np.random.default_rng(seed)
    x0 = np.zeros(offsets[-1] + 1); x0[:-1] = rng.normal(scale=.02, size=offsets[-1])
    x0[-1] = -1.0
    def residual(x):
        out = -x[-1] * np.eye(1, len(rows), 0).ravel() if False else np.zeros(len(rows))
        out[lookup[()]] += x[-1]
        for k, mat in enumerate(maps):
            v = x[offsets[k]:offsets[k+1]]
            out += mat @ np.outer(v, v).reshape(-1)
        return out - target
    def objective(x):
        r = residual(x)
        return float(-x[-1] + np.sqrt(r*r + 1e-10).sum())
    def jacobian(x):
        r = residual(x); weight = r / np.sqrt(r*r + 1e-10)
        grad = np.zeros_like(x); grad[-1] = -1.0 + weight[lookup[()]]
        for k, mat in enumerate(maps):
            v = x[offsets[k]:offsets[k+1]]
            a = (mat.T @ weight).reshape(len(v), len(v))
            grad[offsets[k]:offsets[k+1]] = (a + a.T) @ v
        return grad
    started = time.monotonic()
    result = minimize(objective, x0, jac=jacobian, method="L-BFGS-B", options={"maxiter": iterations, "ftol": 1e-12})
    factors=[]
    for k,b in enumerate(blocks):
        v=result.x[offsets[k]:offsets[k+1]]
        row=[int(round(float(a)*denominator)) for a in v]
        factors.append({"name":b["name"],"words":b["words"],"factor":[row]})
    cert={"modes":modes,"particles":particles,"hamiltonian":encode(h),
          "number_multiplier":[],"b":str(F(int(round(result.x[-1]*denominator)),denominator)),
          "denominator":denominator,"blocks":factors}
    try:
        receipt=verify(cert); status="accepted"
    except Exception as exc:
        receipt={"error":f"{type(exc).__name__}: {exc}"}; status="rejected"
    baseline = dict(cert); baseline["b"] = "0"; baseline["blocks"] = [dict(b, factor=[[0] * len(b["words"])]) for b in factors]
    baseline_receipt = verify(baseline)
    receipt.update({"status":status,"optimizer_success":bool(result.success),"optimizer_message":result.message,
                    "iterations":result.nit,"proposal_objective":float(result.fun),"proposed_b":float(result.x[-1]),
                    "factor_variables":int(offsets[-1]),"coefficient_rows":len(rows),
                    "elapsed_seconds":time.monotonic()-started,"source_gram_used":False,
                    "fock_sector_matrix_used":False,
                    "zero_sos_baseline_lower":baseline_receipt["lower_float"],
                    "zero_sos_baseline_residual_l1":baseline_receipt["residual_l1_float"]})
    return cert,receipt


def main():
    p=argparse.ArgumentParser(); p.add_argument("--modes",type=int,default=4); p.add_argument("--t",default="1/5"); p.add_argument("--iterations",type=int,default=120); p.add_argument("--seed",type=int,default=7); p.add_argument("--fixture",type=Path); p.add_argument("--outputdir",type=Path,default=Path("results/certificate_scaling/direct_factor_optimization")); a=p.parse_args()
    cert,receipt=run(a.modes,F(a.t),a.iterations,a.seed,fixture=a.fixture); a.outputdir.mkdir(parents=True,exist_ok=True)
    (a.outputdir/"candidate.json").write_text(json.dumps(cert)+'\n'); (a.outputdir/"receipt.json").write_text(json.dumps(receipt,indent=2)+'\n'); print(json.dumps(receipt,indent=2))
if __name__=="__main__": main()
