"""Bounded higher-accuracy rerun of the reduced quartic stabilizer SDP."""
import json
from pathlib import Path
from fractions import Fraction as F
from experiments.marginal_stabilizer import split_blocks
from experiments.marginal_adaptive import assemble_and_solve
from experiments.marginal_coefficient import export, symmetric_upper
from experiments.marginal_collective import hopping_polynomial, verify_interval

def run(eps=1e-11, denominator=10**12):
    modes, particles = 10, 5
    blocks, decomposition = split_blocks(modes, 4)
    h = hopping_polynomial(modes, F(1, 5))
    proposal, _ = assemble_and_solve(h, modes, particles, blocks, degree=4,
        groups=[[i, i+5] for i in range(5)], eps=eps, residual_penalty=True)
    cert, receipt = export(h, modes, particles, blocks, proposal,
        denominator=denominator, operator_degree=4)
    cert['variational_upper'] = symmetric_upper(modes, F(1, 5))
    receipt.update(verify_interval(cert))
    receipt.update({'numerical_objective': proposal['numerical_objective'],
                    'residual_penalty': True, 'solver_eps': eps,
                    'factor_denominator': denominator})
    out = Path(__file__).resolve().parents[1] / 'results/marginal_stabilizer_refine'
    out.mkdir(exist_ok=True)
    (out/'certificate.json').write_text(json.dumps(cert)+'\n')
    (out/'receipt.json').write_text(json.dumps(receipt, indent=2)+'\n')
    (out/'structure.json').write_text(json.dumps({'blocks':len(blocks),
        'largest':max(b.get('basis_transform',[]).shape[1] for b in blocks),
        'solver_eps':eps, 'factor_denominator':denominator}, indent=2)+'\n')
    print(json.dumps({k:receipt.get(k) for k in ('lower_float','width_float','residual_l1_float','numerical_objective','status')}), flush=True)
    return receipt

if __name__ == '__main__': run()
