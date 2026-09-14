"""Exact rational dual certificate for the finite gamma=-10 diagnostic."""
from pathlib import Path
import json, sys
import numpy as np
from scipy.optimize import linprog
from sympy import Matrix, linsolve

sys.path.insert(0, str(Path(__file__).parents[3]))
exec((Path(__file__).parent / "charge_metric_limit.py").read_text().split("receipts=[]")[0])
gamma = -10
M = K - gamma * W
A = np.vstack((np.c_[-M, np.ones(len(patterns))], np.c_[-W, np.zeros(len(patterns))]))
r = linprog(np.r_[np.zeros(n), -1.], A_ub=A, b_ub=np.zeros(2*len(patterns)),
            A_eq=np.array([np.r_[mean, 0.]]), b_eq=[1.],
            bounds=[(None, None)]*(n+1), method="highs-ds")
support = np.where(r.ineqlin.marginals[:len(patterns)] < -1e-8)[0].tolist()
mean_num = (np.array(multiplicities) @ W).astype(int).tolist()
total = int(sum(multiplicities))
Ms = Matrix(M[support, :].astype(int)); mn = Matrix(mean_num)
ae = Ms.T.row_join(-mn/total).col_join(Matrix([[1]*len(support)+[0]]))
ys, a = list(linsolve((ae, Matrix([0]*n+[1]))))[0][:-1], list(linsolve((ae, Matrix([0]*n+[1]))))[0][-1]
outdir = Path(__file__).parents[1] / "quadratic_metric_frontier"; outdir.mkdir(exist_ok=True)
cert = {"scope":"finite 1106-pattern diagnostic", "sites":8, "gamma":gamma,
 "patterns":len(patterns), "orbits":orbits, "support":support,
 "y":[str(x) for x in ys], "z":[], "a":str(a), "mean_num":mean_num, "mean_den":total,
 "W":W.astype(int).tolist(), "K":K.astype(int).tolist()}
(outdir/"exact_certificate.json").write_text(json.dumps(cert, separators=(",",":"))+"\n")
(outdir/"exact_receipt.json").write_text(json.dumps({"success":True,"a":str(a),"support_size":len(support),"numerical_floor":float(r.x[-1])}, indent=2)+"\n")
print(json.dumps({"a":str(a),"support":len(support),"patterns":len(patterns)}))
