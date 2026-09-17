"""Direct local boundary SDP proposal, followed by exact rational replay.

No full ladder determinant sector or ground-state teacher is used.
"""
import argparse
from fractions import Fraction as F
import json
from math import floor, sqrt
from pathlib import Path
import time
import numpy as np
from scipy import sparse
import clarabel
from research.local_response_20260916.local_exact import (
    zero, sector_labels, local_blocks, partial_trace, verify, verify_dual)


def encode_matrix(a):
    return [[str(x) for x in row] for row in a]


def rounded(x, scale=10**10):
    if not np.isfinite(x):
        raise ValueError('Nonfinite numerical proposal')
    return F(round(float(x)*scale), scale)


def format_base():
    return {'kind':'local_hubbard_boundary_v1','rungs':4,'U':'8','t':'1'}


def repair_dual(densities, mix):
    # First make each local trace exactly one by a uniform correction.
    normalized=[]
    for rho in densities:
        trace=sum(matrix[i][i] for labels,matrix in rho.values()
                  for i in range(len(labels)))
        fixed={}
        for key,(labels,matrix) in rho.items():
            a=[row[:] for row in matrix]
            for i in range(len(a)):
                a[i][i]+=(1-trace)/256
            fixed[key]=(labels,a)
        normalized.append(fixed)
    left=[partial_trace(rho,0) for rho in normalized]
    right=[partial_trace(rho,1) for rho in normalized]
    targets=[left[0]]
    targets += [[[ (right[i][a][b]+left[i+1][a][b])/2
                   for b in range(16)] for a in range(16)] for i in range(2)]
    targets += [right[-1]]
    output=[]
    for j,rho in enumerate(normalized):
        dl=[[targets[j][a][b]-left[j][a][b] for b in range(16)] for a in range(16)]
        dr=[[targets[j+1][a][b]-right[j][a][b] for b in range(16)] for a in range(16)]
        fixed={}
        for key,(labels,matrix) in rho.items():
            a=zero(len(labels))
            for u,x in enumerate(labels):
                for v,y in enumerate(labels):
                    value=matrix[u][v]
                    if x>>4 == y>>4:
                        value+=dl[x&15][y&15]/16
                    if x&15 == y&15:
                        value+=dr[x>>4][y>>4]/16
                    a[u][v]=(1-mix)*value+(mix/256 if u==v else 0)
            fixed[f'{key[0]},{key[1]}']=encode_matrix(a)
        output.append(fixed)
    return output


def run(out):
    start=time.monotonic()
    if out.exists():
        raise FileExistsError(out)
    out.mkdir(parents=True)
    bare=local_blocks([zero(16),zero(16)])
    rung_coordinates=[(a,b) for labels in sector_labels(2).values()
                      for j,b in enumerate(labels) for a in labels[:j+1]]
    coords={pair:i for i,pair in enumerate(rung_coordinates)}
    nv=3+2*len(coords)
    row_indices=[];col_indices=[];values=[];rhs=[];cones=[];layout=[]
    # PSD slack h_i + right J_i - left J_(i-1) - lambda_i I.
    for block_id,blocks in enumerate(bare):
        for key,(labels,matrix) in blocks.items():
            startrow=len(rhs)
            for j,y in enumerate(labels):
                for i,x in enumerate(labels[:j+1]):
                    scale=1 if i==j else sqrt(2)
                    row=len(rhs);rhs.append(scale*float(matrix[i][j]))
                    if i==j:
                        row_indices.append(row);col_indices.append(block_id);values.append(1.)
                    if block_id and x>>4 == y>>4:
                        pair=tuple(sorted((x&15,y&15)))
                        if pair in coords:
                            row_indices.append(row)
                            col_indices.append(3+(block_id-1)*len(coords)+coords[pair])
                            values.append(scale)
                    if block_id<2 and x&15 == y&15:
                        pair=tuple(sorted((x>>4,y>>4)))
                        if pair in coords:
                            row_indices.append(row)
                            col_indices.append(3+block_id*len(coords)+coords[pair])
                            values.append(-scale)
            cones.append(clarabel.PSDTriangleConeT(len(labels)))
            layout.append((block_id,key,labels,startrow,len(rhs)))
    A=sparse.csc_matrix((values,(row_indices,col_indices)),shape=(len(rhs),nv))
    objective=np.zeros(nv);objective[:3]=-1
    settings=clarabel.DefaultSettings()
    settings.verbose=False;settings.max_iter=150
    settings.tol_gap_abs=1e-9;settings.tol_feas=1e-9;settings.tol_gap_rel=1e-9
    settings.time_limit=100
    # Keep submitted cones intact so returned duals have the declared ordering.
    settings.chordal_decomposition_enable=False
    settings.max_threads=1
    solver=clarabel.DefaultSolver(sparse.csc_matrix((nv,nv)),objective,A,
                                 np.array(rhs),cones,settings)
    result=solver.solve()
    discovery_seconds=time.monotonic()-start
    print(json.dumps({'solver_status':str(result.status),'objective':result.obj_val,
                      'seconds':discovery_seconds}),flush=True)
    if str(result.status) not in ('Solved','AlmostSolved'):
        raise RuntimeError('Numerical optimization failed to propose a solution')
    x=np.array(result.x);z=np.array(result.z)
    messages=[]
    for j in range(2):
        matrix=zero(16)
        for (a,b),k in coords.items():
            matrix[a][b]=matrix[b][a]=rounded(x[3+j*len(coords)+k])
        messages.append(matrix)
    candidate=format_base()
    candidate['messages']=[encode_matrix(m) for m in messages]
    shifted=local_blocks(messages)
    minima=[min(float(np.linalg.eigvalsh(np.array(m,dtype=float))[0])
                for labels,m in blocks.values()) for blocks in shifted]
    # Exact acceptance decides; these float eigenvalues only suggest shifts.
    attempt_log=[]
    for margin in (1e-8,1e-7,1e-6,1e-5):
        candidate['lower_shifts']=[str(F(floor((v-margin)*10**10),10**10)) for v in minima]
        try:
            primal=verify(candidate)
            attempt_log.append({'primal_margin':margin,'accepted':True})
            break
        except ValueError as error:
            attempt_log.append({'primal_margin':margin,'accepted':False,'reason':str(error)})
    else:
        raise RuntimeError('Primal rational repair failed')
    densities=[{} for _ in range(3)]
    for bi,key,labels,begin,end in layout:
        matrix=zero(len(labels));position=begin
        for j in range(len(labels)):
            for i in range(j+1):
                value=rounded(z[position]/(1 if i==j else sqrt(2)))
                matrix[i][j]=matrix[j][i]=value;position+=1
        assert position==end
        densities[bi][key]=(labels,matrix)
    for mix in (F(1,10**7),F(1,10**6),F(1,10**5),F(1,10**4),F(1,1000)):
        candidate['dual_marginals']=repair_dual(densities,mix)
        try:
            dual=verify_dual(candidate)
            attempt_log.append({'dual_identity_mixture':str(mix),'accepted':True})
            break
        except ValueError as error:
            attempt_log.append({'dual_identity_mixture':str(mix),'accepted':False,'reason':str(error)})
    else:
        raise RuntimeError('Dual exact repair failed')
    # Matched zero-message baseline; independent rational acceptance as above.
    baseline=format_base();baseline['messages']=[encode_matrix(zero(16)) for _ in range(2)]
    base_minima=[min(float(np.linalg.eigvalsh(np.array(m,dtype=float))[0])
                     for labels,m in blocks.values()) for blocks in bare]
    baseline['lower_shifts']=[str(F(floor((v-1e-7)*10**10),10**10)) for v in base_minima]
    baseline_result=verify(baseline)
    receipt={'status':'primal_and_family_ceiling_exactly_accepted',
             'primal':primal,'dual':dual,'zero_message_baseline':baseline_result,
             'optimization_variables':nv,'boundary_coordinates_per_rung':len(coords),
             'PSD_blocks':len(cones),'coefficient_rows':len(rhs),'map_nonzeros':A.nnz,
             'solver_status':str(result.status),'solver_iterations':result.iterations,
             'discovery_seconds':discovery_seconds,'total_internal_seconds':time.monotonic()-start,
             'repair_attempts':attempt_log,'global_ground_state_teacher_used':False,
             'global_determinants_enumerated':0,'max_local_fock_dimension':256}
    receipt['exact_boundary_gain_over_t']=str(F(primal['energy_lower_over_t'])-F(baseline_result['energy_lower_over_t']))
    receipt['certified_optimization_bracket_width_over_t']=str(F(dual['family_lower_ceiling_over_t'])-F(primal['energy_lower_over_t']))
    (out/'certificate.json').write_text(json.dumps(candidate,separators=(',',':'))+'\n')
    (out/'baseline.json').write_text(json.dumps(baseline,separators=(',',':'))+'\n')
    (out/'construction.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True)
    run(p.parse_args().out)
