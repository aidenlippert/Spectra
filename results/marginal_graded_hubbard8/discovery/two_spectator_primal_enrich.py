"""Bounded primal subspace SDP with full-block spectral enrichment; nonaccepting."""
from pathlib import Path
from fractions import Fraction as F
import argparse
import gzip
import hashlib
import json
import time
import numpy as np
import cvxpy as cp
from scipy.linalg import eigh
from two_spectator_family_trust import (prepare, BASE, ROOT, actions, projected_matrix,
    spin_actions, spectator_actions, pair_actions, two_actions, LABELS, SPECTATORS, PAIRS, TWO)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('source', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    cert_path = args.source / 'profile_joint_r1_2_certificate.json'
    proposal_path = args.source / 'diagonal_family_limit_proposal.json'
    ledger_path = args.source / 'candidate_ledger.json.gz'
    c = json.loads(cert_path.read_text())
    proposed = json.loads(proposal_path.read_text())
    with gzip.open(ledger_path, 'rt') as f:
        ledger = json.load(f)
    if len(ledger['states']) > 13000:
        raise ValueError('Bounded candidate input required')
    canonical = lambda v: tuple(sorted((int(s),a) for s,a in v.items() if a))
    wanted = {canonical(item['vector']) for item in proposed['mixture']}
    selected = {}
    found = set()
    anchors = []
    for i, state in enumerate(ledger['states']):
        key = canonical(state['vector'])
        if state['sector'] is None:
            anchors.append(i)
        if key in wanted:
            found.add(key)
            if state['sector'] is not None:
                selected.setdefault(tuple(state['sector']), []).append(state)
    if found != wanted or len(selected) > 32:
        raise ValueError('Complete bounded selected physical spans required')
    candidates_path = BASE / 'joint_projector/signed_density/symmetry_diagonal_candidates.json'
    candidates = json.loads(candidates_path.read_text())['candidates']
    shapes = [{int(s):v for s,v in candidates[i]['diagonal'].items()} for i in [1,2,4,5,6,7,31,11,20]]
    x, _, mats, *_ = prepare(c, shapes)
    action = [actions()] + [spin_actions({k:1}) for k in LABELS] + [spectator_actions({k:1}) for k in SPECTATORS] + [pair_actions({k:1}) for k in PAIRS] + [two_actions({k:1}) for k in TWO]
    rhs = np.array(list(map(F, ledger['rhs'])), float)
    if len(rhs) != 119:
        raise ValueError('119-row target required')
    anchor_rows = np.array([ledger['rows'][i] for i in anchors], float).T
    anchor_energy = np.array([ledger['energies'][i] for i in anchors], float)
    bases={};prepared={}
    for key,a,ds,ph,q,ts,cols,*_ in mats:
        scale=np.sqrt([sum(v*v for v in col.values()) for col in cols])
        hop=np.array([projected_matrix(cols,item) for item in action],float)/scale[:,None]/scale[None,:]
        prepared[key]=(a-np.einsum('i,ijk->jk',np.array(x[2:4],float),ds),ds,ph,q,ts,hop)
        if key not in selected:continue
        vectors=np.array([np.array(s['coeff'],float)*scale/np.sqrt(float(s['norm'])) for s in selected[key]]).T
        u,singular,_=np.linalg.svd(vectors,full_matrices=False);rank=int(sum(singular>1e-10))
        if not 1<=rank<=48:raise ValueError('Initial compressed PSD dimension cap48 exceeded')
        bases[key]=u[:,:rank]
    history=[]
    for iteration in range(4):
        mass=cp.Variable(len(anchors),nonneg=True);moments=anchor_rows@mass;objective=anchor_energy@mass
        blocks=[];constraints=[]
        for key,v in bases.items():
            a,ds,ph,q,ts,hop=prepared[key];rank=v.shape[1]
            operators=[np.eye(len(a))]+list(ds)+[np.diag(t) for t in ts]+list(hop)+[ph,q]
            coefficients=np.array([(v.T@op@v).ravel(order='C') for op in operators]);fixed=v.T@a@v
            density=cp.Variable((rank,rank),symmetric=True);constraints.append(density>>0)
            moments=moments+coefficients@cp.vec(density,order='C');objective=objective+cp.sum(cp.multiply(fixed,density))
            blocks.append((key,v,density,rank))
        equality=moments[:117]==rhs[:117];inequality=moments[117:]<=rhs[117:]
        constraints += [equality,inequality]
        problem=cp.Problem(cp.Minimize(objective/5),constraints)
        problem.solve(solver='CLARABEL',time_limit=30.,max_iter=80,tol_gap_abs=1e-9,tol_feas=1e-9,verbose=False)
        entry={'iteration':iteration,'status':problem.status,'objective':problem.value,'solver_seconds':problem.solver_stats.solve_time,'blocks':len(blocks),'dimensions':[b[3] for b in blocks]}
        if equality.dual_value is None or inequality.dual_value is None:raise ValueError('No numerical conic dual for enrichment')
        dual=-5*np.r_[equality.dual_value,inequality.dual_value];spectra=[]
        for key,(a,ds,ph,q,ts,hop) in prepared.items():
            price=a-dual[0]*np.eye(len(a))-np.einsum('i,ijk->jk',dual[1:3],ds)-np.diag(dual[3:64]@ts)-np.einsum('i,ijk->jk',dual[64:117],hop)-dual[117]*ph-dual[118]*q
            ev,vec=eigh(price,subset_by_index=[0,min(1,len(a)-1)])
            spectra.append((float(ev[0]),key,ev,vec))
        entry['minimum_full_reduced_eigenvalue']=min(item[0] for item in spectra)
        entry['numerical_full_dual_lower']=(float(rhs@dual)+min(0,entry['minimum_full_reduced_eigenvalue']))/5
        entry['new_basis_directions']=0;entry['limited_blocks']=0
        if iteration<3:
            for lowest,key,ev,vec in sorted(spectra,key=lambda z:z[0]):
                if lowest>=-1e-8:continue
                if key not in bases and len(bases)>=32:entry['limited_blocks']+=1;continue
                old=bases.get(key,np.zeros((len(vec),0)));joined=np.column_stack([old,vec[:,ev<-1e-8]])
                u,singular,_=np.linalg.svd(joined,full_matrices=False);rank=int(sum(singular>1e-8))
                if rank>48:entry['limited_blocks']+=1;continue
                bases[key]=u[:,:rank];entry['new_basis_directions']+=rank-old.shape[1]
        history.append(entry);print(json.dumps(entry),flush=True)
    atoms = []
    minimum_psd = 0.
    for key,v,density,rank in blocks:
        if density.value is None:
            continue
        ev,vec = np.linalg.eigh(density.value)
        minimum_psd = min(minimum_psd,float(ev.min()))
        for i,value in enumerate(ev):
            if value > 1e-10:
                atoms.append({'sector':list(key),'dual_weight':float(value),'coordinates':(v@vec[:,i]).tolist()})
    (args.output / 'dual_atoms.json').write_text(json.dumps({'accepted':False,'atoms':atoms,'solver_status':problem.status},indent=2)+'\n')
    residual = None if moments.value is None else np.asarray(moments.value)-rhs
    files = [Path(__file__), Path(__file__).with_name('two_spectator_family_trust.py'),
             Path(__file__).with_name('signed_charge_numeric.py'), cert_path, proposal_path, ledger_path, candidates_path]
    diagnostic = {'accepted':False,'status':problem.status,'numerical_objective':problem.value,
                  'seconds':time.monotonic()-started,'solver_seconds':problem.solver_stats.solve_time,
                  'blocks':len(blocks),'dimensions':[b[3] for b in blocks],'atoms':len(atoms),
                  'minimum_numerical_psd_eigenvalue':minimum_psd,
                  'maximum_equality_residual':None if residual is None else float(np.max(np.abs(residual[:117]))),
                  'maximum_inequality_residual':None if residual is None else float(np.max(residual[117:])),
                  'history':history,'source_sha256':{str(p.resolve().relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
                  'scope':'Four numerical primal SDP rounds with full-block spectral enrichment, at most32 spans of dimension48 plus determinant anchors.30 solver seconds and80iterations per round. Bounds enforced and limited enrichments recorded. Solver atoms and residuals are untrusted; exact physical source reconstruction and family replay required.'}
    (args.output / 'subspace_diagnostic.json').write_text(json.dumps(diagnostic,indent=2)+'\n')
    print(json.dumps({k:v for k,v in diagnostic.items() if k!='source_sha256'}),flush=True)


if __name__ == '__main__':
    main()
