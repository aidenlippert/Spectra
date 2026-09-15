"""Positivity discovery for a fixed charge metric, with independent LP blocks.

Metric provenance is explicit: this is not fresh metric discovery. All atom
columns are generated/priced algebraically and exact polynomial replay decides
acceptance. A numerical failure never proves cone infeasibility.
"""
from fractions import Fraction as F
from pathlib import Path
import json
import time
from experiments.marginal_joint_coefficient_constructor import input_digest
from experiments.marginal_moment_pricing import MomentDictionary
from experiments.marginal_joint_spinflip import flip_label
from experiments.marginal_polynomial_metric import JointPolynomial, check_positive, replay
from experiments.marginal_number_quotient import complete_bounded_number_ideals


def export_part(ring, polynomial, scale, labels, weights, bound, denominator=10**12):
    proof = {'denominator': denominator, 'bound': round(F(bound)*denominator),
             'positive_indicators': [], 'charge_indicators': [], 'number_multipliers': [[], []]}
    represented = {0: F(proof['bound'], denominator)}
    localizer = {0: -1, **{3 << (2*i): 1 for i in range(ring.sites)}}
    for label, value in zip(labels, weights):
        orbit = sorted({tuple(label), tuple(flip_label(label, ring.sites))})
        integer = max(0, round(float(value)*denominator/len(orbit)))
        if not integer:
            continue
        for family, required, occupied in orbit:
            if family not in ('positive', 'charge'):
                raise ValueError('Unsupported atom family')
            key = 'positive_indicators' if family == 'positive' else 'charge_indicators'
            proof[key].append({'required': required, 'occupied': occupied, 'weight': integer})
            term = {occupied: F(integer, denominator)}
            for i in range(ring.modes):
                if (required ^ occupied) & (1 << i):
                    term = ring.multiply(term, {0: 1, 1 << i: -1})
            if family == 'charge':
                term = ring.multiply(term, localizer)
            represented = ring.add(represented, term)
    residual = ring.add({m: F(v, scale) for m, v in polynomial.items()}, {m: -v for m, v in represented.items()})
    remainder, ideals, stats = complete_bounded_number_ideals(ring, residual)
    proof['number_multipliers'] = [[{'mask': m, 'coefficient': round(v*denominator)} for m, v in sorted(ideal.items()) if round(v*denominator)] for ideal in ideals]
    receipt = check_positive(ring, polynomial, scale, proof)
    receipt.update(quotient_remainder_l1=str(sum(map(abs, remainder.values()), F(0))), quotients=stats)
    if F(receipt['lower']) <= 0:
        raise ValueError('Exact positive part bound was not established')
    return proof, receipt


class IncrementalMaster:
    """Keep the simplex basis while appending zero-cost positivity columns."""
    def __init__(self, rhs):
        import highspy
        import numpy as np
        from scipy.sparse import hstack, eye
        self.highspy = highspy; self.rows = len(rhs); self.columns = 0
        self.solver = highspy.Highs()
        for key, value in [('threads', 1), ('output_flag', False), ('solver', 'simplex'),
                           ('simplex_strategy', 4), ('primal_feasibility_tolerance', 1e-8),
                           ('dual_feasibility_tolerance', 1e-8)]:
            if self.solver.setOptionValue(key, value) != highspy.HighsStatus.kOk:
                raise ValueError('Native option rejected: '+key)
        matrix = hstack((eye(self.rows), -eye(self.rows)), format='csc')
        lp = highspy.HighsLp(); lp.num_row_ = self.rows; lp.num_col_ = 2*self.rows
        lp.col_cost_ = np.ones(2*self.rows); lp.col_lower_ = np.zeros(2*self.rows)
        lp.col_upper_ = np.full(2*self.rows, highspy.kHighsInf)
        lp.row_lower_ = rhs; lp.row_upper_ = rhs
        lp.a_matrix_.format_ = highspy.MatrixFormat.kColwise
        lp.a_matrix_.start_, lp.a_matrix_.index_, lp.a_matrix_.value_ = matrix.indptr, matrix.indices, matrix.data
        if self.solver.passModel(lp) == highspy.HighsStatus.kError:
            raise ValueError('Native residual master rejected')

    def run(self, fresh, remaining):
        import numpy as np
        from types import SimpleNamespace
        h = self.highspy
        fresh = fresh.tocsc(); count = fresh.shape[1]
        if fresh.shape[0] != self.rows:
            raise ValueError('New master columns have wrong row count')
        if count:
            status = self.solver.addCols(count, np.zeros(count), np.zeros(count),
                np.full(count, h.kHighsInf), fresh.nnz, fresh.indptr, fresh.indices, fresh.data)
            if status == h.HighsStatus.kError:
                raise ValueError('Native master rejected appended columns')
            self.columns += count
        self.solver.setOptionValue('time_limit', self.solver.getRunTime()+remaining)
        self.solver.run(); status = self.solver.getModelStatus(); solution = self.solver.getSolution()
        values = np.asarray(solution.col_value)
        return SimpleNamespace(success=status == h.HighsModelStatus.kOptimal,
            message=self.solver.modelStatusToString(status),
            fun=float(self.solver.getInfo().objective_function_value),
            x=np.r_[values[2*self.rows:], values[:2*self.rows]],
            eqlin=SimpleNamespace(marginals=np.asarray(solution.row_dual)))


def solve_part(candidate, name, out, batch=256, max_rounds=80, time_limit=180, bound=None, native=False):
    import numpy as np
    from scipy.sparse import csc_matrix, hstack, eye
    from scipy.optimize import linprog
    if set(candidate) != {'modes','particles','hamiltonian','target_lower','polynomial_metric'}:
        raise ValueError('Fixed bare-H metric candidate required')
    if name not in ('weight','numerator') or batch < 1 or max_rounds < 1 or time_limit <= 0:
        raise ValueError('Valid part and positive search budgets required')
    out=Path(out);out.mkdir(parents=True,exist_ok=True);started=time.monotonic()
    bound=F(1,10) if name=='weight' and bound is None else F(1) if bound is None else F(bound)
    if bound <= 0:
        raise ValueError('Positive export margin required')
    ring=JointPolynomial(candidate);w,k,cost=ring.compile(F(candidate['target_lower']))
    polynomial,scale=(w,cost['metric_scale']) if name=='weight' else (k,cost['numerator_scale'])
    degree=max(m.bit_count() for m in polynomial)
    data={key:candidate[key] for key in ('modes','particles','hamiltonian')}
    dictionary=MomentDictionary(data,F(candidate['target_lower']),feature_degree=0,proof_degree=degree)
    normal=dictionary.quotient.normal({m:F(v,scale) for m,v in polynomial.items()})
    rhs=dictionary.projection@np.array([float(normal.get(m,0)) for m in dictionary.quotient.basis])-float(bound)*dictionary.one
    rows=dictionary.qrows;labels=dictionary.labels(2);active={(0,label) for label in labels};matrix=csc_matrix((rows,0));built=0;identity=eye(rows,format='csc')
    native_master=IncrementalMaster(rhs) if native else None
    prepared=time.monotonic();history=[];accepted=False;reason='Round budget reached'
    for iteration in range(max_rounds):
        remaining=time_limit-(time.monotonic()-prepared)
        if remaining<=0:
            reason='Time budget reached';break
        fresh=[dictionary.column(label) for label in labels[built:]]
        fresh_matrix=csc_matrix(np.array(fresh).T) if fresh else csc_matrix((rows,0))
        if fresh:matrix=hstack((matrix,fresh_matrix),format='csc')
        built=len(labels);master=hstack((matrix,identity,-identity),format='csc')
        if native_master:
            result=native_master.run(fresh_matrix,remaining)
        else:
            result=linprog(np.r_[np.zeros(built),np.ones(2*rows)],A_eq=master,b_eq=rhs,bounds=(0,None),method='highs-ds',options={'time_limit':remaining,'dual_feasibility_tolerance':1e-8,'primal_feasibility_tolerance':1e-8})
        record={'round':iteration,'active_atoms':built,'master_columns':master.shape[1],'master_nonzeros':master.nnz,'status':result.message,'search_seconds':time.monotonic()-prepared}
        if not result.success:
            history.append(record);reason='Restricted solver did not finish';break
        record['phase_one_l1']=float(result.fun)
        if result.fun<1e-6:
            proposal={'candidate_sha256':input_digest(candidate),'labels':labels,'values':result.x[:built].tolist(),'bound':str(bound)}
            (out/'proposal.json').write_text(json.dumps(proposal,indent=2)+'\n')
            try:
                proof,exact=export_part(ring,polynomial,scale,labels,result.x[:built],bound)
                (out/'proof.json').write_text(json.dumps(proof,indent=2)+'\n');record['exact_receipt']=exact;accepted=True;reason='Exact part accepted'
            except ValueError as error:
                record['export_rejection']=str(error);reason='Exact export rejected'
            history.append(record);print(json.dumps(record),flush=True);break
        # Dictionary uses -atom columns; negate our equality dual for pricing.
        fresh,maximum,checked=dictionary.price(np.r_[-result.eqlin.marginals,np.zeros(rows+1)],active,batch)
        selected=[label for block,label in fresh if block==0]
        record.update(maximum_dual_violation=maximum,candidates_priced=checked,added_atoms=len(selected))
        history.append(record);print(json.dumps(record),flush=True)
        (out/'history.json').write_text(json.dumps(history,indent=2)+'\n')
        if not selected:
            reason='No numerical improving column; no exact infeasibility claim';break
        labels.extend(selected);active.update((0,label) for label in selected)
    receipt={'candidate_sha256':input_digest(candidate),'part':name,'bound':str(bound),'proof_degree':degree,'exact_accepted':accepted,'reason':reason,'rows':rows,'active_atoms':len(labels),'preparation_seconds':prepared-started,'total_seconds':time.monotonic()-started,'rounds':len(history),'native_incremental':native,'solver_version':native_master.solver.version() if native_master else 'SciPy bundled HiGHS','dictionary':dictionary.stats,'scope':'Fixed imported metric; generated geometry-only seeds and direct moment pricing. No configuration enumeration or full atom matrix. Exact part replay required, then full joint replay for combined acceptance.'}
    (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');(out/'history.json').write_text(json.dumps(history,indent=2)+'\n');(out/'active_labels.json').write_text(json.dumps(labels)+'\n');print(json.dumps(receipt),flush=True)
    return receipt


def combine(candidate, out):
    out=Path(out)
    certificate=dict(candidate,kind='joint_polynomial_metric_gap_v1')
    for name in ('weight','numerator'):
        receipt=json.loads((out/name/'receipt.json').read_text())
        if receipt['candidate_sha256']!=input_digest(candidate) or not receipt['exact_accepted']:
            raise ValueError('Both exact parts must bind this candidate')
        certificate[name+'_proof']=json.loads((out/name/'proof.json').read_text())
    result=replay(certificate)
    (out/'certificate.json').write_text(json.dumps(certificate,indent=2)+'\n');(out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


if __name__=='__main__':
    import argparse
    from unittest.mock import patch
    parser=argparse.ArgumentParser();parser.add_argument('--candidate',type=Path,required=True);parser.add_argument('--out',type=Path,required=True);parser.add_argument('--part',choices=('weight','numerator'));parser.add_argument('--time-limit',type=float,default=180);parser.add_argument('--native',action='store_true');args=parser.parse_args()
    candidate=json.loads(args.candidate.read_text())
    with patch('experiments.marginal_determinant_tree.DeterminantOracle.action',side_effect=AssertionError('No states')),patch('experiments.marginal_spin_constructor.spin_states',side_effect=AssertionError('No states')),patch('experiments.marginal_polynomial_metric.complete_number_ideals',side_effect=AssertionError('No full lift')):
        if args.part:solve_part(candidate,args.part,args.out/args.part,time_limit=args.time_limit,native=args.native)
        else:
            results=[solve_part(candidate,name,args.out/name,time_limit=args.time_limit,native=args.native) for name in ('weight','numerator')]
            if all(r['exact_accepted'] for r in results):print(json.dumps(combine(candidate,args.out)),flush=True)
