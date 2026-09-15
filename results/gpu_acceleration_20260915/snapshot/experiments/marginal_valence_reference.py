"""Valence-projector reference, weighted DD charge gap, and fresh Temple proof."""
from fractions import Fraction as F
import json
from pathlib import Path
import time

from experiments.marginal_spin_reduction import SpinZeroOracle, symmetrize
from experiments.marginal_symbolic import encode, add, scale
from experiments.marginal_valence_states import valence_states, dimer_singlet_witness, spin_raise_residual
from experiments.marginal_localized_search import construct as seed_gap
from experiments.marginal_sparse_upper import refine as refine_upper, witness_statistics
from experiments.marginal_spin_temple import discover, replay


def build(source, output, *, denominator=10**14, upper_steps=24, target=F(1,10**8), offset=F(1,100)):
    source,output=Path(source),Path(output);target,offset=F(target),F(offset)
    if output.exists(): raise ValueError('Preserve previous valence-reference export')
    if type(denominator) is not int or not 1<=denominator<=10**18:
        raise ValueError('Bounded positive coefficient denominator required')
    if type(upper_steps) is not int or not 1<=upper_steps<=32 or target<=0 or offset<=0:
        raise ValueError('Positive target, offset and bounded upper budget required')
    input_data=json.loads(source.read_text());original={k:input_data[k] for k in ('modes','particles','hamiltonian')}
    original_oracle=SpinZeroOracle(original)
    p=valence_states(original_oracle.modes,original_oracle.particles)
    initial=dimer_singlet_witness(original_oracle.modes,original_oracle.particles)
    if spin_raise_residual(initial,original_oracle.modes): raise ValueError('Initial dimer state is not an exact singlet')
    rounded={w:F(round(x*denominator),denominator) for w,x in original_oracle.h.items() if round(x*denominator)}
    h=symmetrize(rounded,original_oracle.modes)
    nearby=dict(original,hamiltonian=encode(h));oracle=SpinZeroOracle(nearby)
    error=sum((abs(x) for x in add(original_oracle.h,scale(h,-1)).values()),F(0))
    if target<=2*error: raise ValueError('Target below certified coefficient transfer budget')
    output.mkdir(parents=True);started=time.monotonic()
    def write(name,data): (output/name).write_text(json.dumps(data,indent=2)+'\n')
    write('hamiltonian.json',nearby)
    write('selection.json',{'retained_states':p,'independent_upper':initial,
        'selection':'All singly occupied spatial sites at half filling; no Hamiltonian-sector diagonalization selects P.',
        'upper_seed':'Product of adjacent-site spin singlets, checked by exact spin raising.'})
    write('quantization.json',{'source':str(source),'denominator':denominator,'coefficient_norm_error':str(error),
        'coefficient_norm_error_float':float(error),'target':str(target),
        'scope':'Nearby rational spin-symmetric Hamiltonian proposed for arithmetic; final replay recomputes the exact coefficient-norm transfer to the input Hamiltonian.'})
    try:
        print(json.dumps({'phase':'valence_upper','retained_dimension':len(p),'transfer_error':float(error)}),flush=True)
        refine_upper(output/'hamiltonian.json',output/'selection.json',output/'upper',upper_steps)
        upper_cert=json.loads((output/'upper/certificate.json').read_text())
        mu,variance,_=witness_statistics(oracle,upper_cert['independent_upper'])
        print(json.dumps({'phase':'valence_charge_gap','upper':float(mu),'variance':float(variance)}),flush=True)
        seed_gap(output/'hamiltonian.json',output/'charge_gap',p,mu+offset)
        gap=json.loads((output/'charge_gap/certificate.json').read_text());gamma=F(gap['target_lower'])
        if not gamma>mu+offset: raise ValueError('Weighted DD charge gap misses the excitation proposal')
        reduced=dict(nearby,retained_states=p,independent_upper=upper_cert['independent_upper'],complement_lower=str(gamma),
            complement_atoms={'kind':gap['kind'],'blocks':gap['blocks']})
        print(json.dumps({'phase':'valence_temple','complement_lower':float(gamma),'upper':float(mu)}),flush=True)
        certificate,receipt,history=discover(original,reduced,oracle,offset,target)
        receipt.update(source=str(source),selection_dimension=len(p),gap_has_dense_factor=False,
                       gap_explicit_atom_count=sum(len(b['atoms']) for b in gap['blocks']),
                       construction_seconds=time.monotonic()-started,
                       construction_scope='Exact finite input-Hamiltonian energy interval using an occupation-selected valence reference, fresh singlet Krylov witness, complete weighted DD charge gap, fresh directional Schur inertia and Temple variance. Rounding error is recomputed in replay. Full spin-sector actions and comparison weighting remain explicit; no general scaling claim.')
        write('certificate.json',certificate);write('receipt.json',receipt);write('response_history.json',history)
        print(json.dumps({k:receipt[k] for k in ('width_float','retained_dimension','response_dimension','construction_seconds')}),flush=True)
        return receipt
    except Exception as exc:
        write('failure.json',{'status':'not_accepted','error':str(exc),'elapsed_seconds':time.monotonic()-started})
        raise


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--source',required=True);p.add_argument('--output',required=True)
    p.add_argument('--denominator',type=int,default=10**14);p.add_argument('--upper-steps',type=int,default=24)
    p.add_argument('--target',default='1/100000000');p.add_argument('--offset',default='1/100')
    args=p.parse_args();build(args.source,args.output,denominator=args.denominator,upper_steps=args.upper_steps,target=F(args.target),offset=F(args.offset))
