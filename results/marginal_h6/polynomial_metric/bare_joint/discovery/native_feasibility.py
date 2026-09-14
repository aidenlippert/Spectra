"""Resume the bare-H coefficient model with native HiGHS and a fixed positive margin."""
from pathlib import Path
from fractions import Fraction as F
from unittest.mock import patch
import json,time
import numpy as np
from scipy.sparse import load_npz
import highspy
from experiments.marginal_joint_coefficient_constructor import input_digest,export
root=Path('results/marginal_h6/polynomial_metric/bare_joint');prepared=root/'search';out=root/'native_feasibility';out.mkdir(exist_ok=True);data=json.loads((root/'hamiltonian.json').read_text());meta=json.loads((prepared/'prepared.json').read_text());assert input_digest(data)==meta['source_sha256'];matrix=load_npz(prepared/'model.npz');rhs=np.load(prepared/'rhs.npy',allow_pickle=False);orbits=meta['orbits'];labels=meta['labels'];floor=F(1,200);n=len(orbits);a=len(labels);assert matrix.shape==(len(rhs),n+2*a+1)
lp=highspy.HighsLp();lp.num_col_=matrix.shape[1];lp.num_row_=matrix.shape[0];lp.col_cost_=np.zeros(matrix.shape[1]);lower=np.zeros(matrix.shape[1]);upper=np.full(matrix.shape[1],highspy.kHighsInf);lower[:n]=-highspy.kHighsInf;lower[-1]=upper[-1]=float(floor);lp.col_lower_=lower;lp.col_upper_=upper;lp.row_lower_=rhs;lp.row_upper_=rhs;lp.a_matrix_.format_=highspy.MatrixFormat.kColwise;lp.a_matrix_.start_=matrix.indptr;lp.a_matrix_.index_=matrix.indices;lp.a_matrix_.value_=matrix.data
solver=highspy.Highs()
for key,value in [('output_flag',True),('threads',1),('solver','simplex'),('simplex_strategy',1),('presolve','on'),('time_limit',180.),('primal_feasibility_tolerance',1e-8),('dual_feasibility_tolerance',1e-8)]:
 if solver.setOptionValue(key,value)!=highspy.HighsStatus.kOk:raise ValueError('Rejected native option '+key)
if solver.passModel(lp)==highspy.HighsStatus.kError:raise ValueError('Native model rejected')
start=time.monotonic();solver.run();status=solver.getModelStatus();solution=solver.getSolution();info=solver.getInfo();receipt={'solver_version':solver.version(),'method':'dual simplex feasibility with fixed numerator floor','floor':str(floor),'status':solver.modelStatusToString(status),'seconds':time.monotonic()-start,'iterations':info.simplex_iteration_count,'source_sha256':input_digest(data),'prepared_source':str(prepared),'shape':list(matrix.shape),'nonzeros':matrix.nnz,'value_valid':solution.value_valid,'scope':'Same bare-H model, no imported metric or atom weights; fixed positive numerator floor rather than bound optimization. Exact export is the acceptance gate.'};accepted=False
if solution.value_valid:
 x=np.asarray(solution.col_value);receipt.update(numerical_floor=float(x[-1]),minimum_atom=float(min(x[n:-1])),max_equation_residual=float(max(abs(matrix@x-rhs))))
 if x[-1]>0:
  proposal={'source_sha256':input_digest(data),'orbits':orbits,'metric_coefficients':x[:n].tolist(),'metric_bound':.001,'numerator_bound':float(x[-1]),'gamma':meta['gamma']}
  for name,values in [('weight',x[n:n+a]),('numerator',x[n+a:n+2*a])]:
   pairs=[(l,float(v)) for l,v in zip(labels,values) if abs(v)>1e-13];proposal[name+'_labels']=[p[0] for p in pairs];proposal[name+'_values']=[p[1] for p in pairs]
  (out/'proposal.json').write_text(json.dumps(proposal,indent=2)+'\n')
  try:
   with patch('experiments.marginal_determinant_tree.DeterminantOracle.action',side_effect=AssertionError('No physical states')),patch('experiments.marginal_spin_constructor.spin_states',side_effect=AssertionError('No state generation')),patch('experiments.marginal_polynomial_metric.complete_number_ideals',side_effect=AssertionError('No full-population lift')):proof=export(data,proposal,out/'proof')
   accepted=True;receipt['exact_receipt']=proof
  except ValueError as error:receipt['exact_export_rejection']=str(error)
receipt['exact_accepted']=accepted;(out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt),flush=True)
