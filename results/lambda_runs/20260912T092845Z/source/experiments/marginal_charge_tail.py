"""Exact brackets for the independently maximized weighted charge-tail envelope.

This brackets a conservative row polynomial on D>=minimum_doublons. It is
not a ground-energy upper bound or a certificate for the complete Q space.
"""
from fractions import Fraction as F
import json
from pathlib import Path
from experiments.marginal_coherent_tree import CoherentCharge
from experiments.marginal_charge_polynomial import add, multiply, scale, indicator, _atoms, _mask
from experiments.marginal_charge_ratio_dp import ratio_upper


def envelope(data):
    o=CoherentCharge(data);minimum=data.get('minimum_doublons')
    if type(minimum) is not int or not 1<=minimum<=o.sites//2:
        raise ValueError('Valid positive ionic-tail threshold required')
    poly=dict(o.oracle.diagonal)
    totals={'ratio_dp_states':0,'ratio_dp_transitions':0,'ratio_metric_powers':0,
            'maximum_memory':0,'empty_source_events':0,'transition_groups':0}
    for (c,a),amp in o.groups.items():
        closed=o.close(c|a,a)
        if closed is None:continue
        lo,hi=o.amplitude_range(amp,*closed)
        if lo>=0:sign=1
        elif hi<=0:sign=-1
        else:raise ValueError('Tail envelope requires fixed-sign grouped amplitudes')
        delta=tuple(((c>>(2*i))&3).bit_count()-((a>>(2*i))&3).bit_count() for i in range(o.sites))
        ratio,cost=ratio_upper(o,delta,*closed,minimum)
        for target,source in [('ratio_dp_states','visited_states'),('ratio_dp_transitions','transitions'),('ratio_metric_powers','metric_scalar_powers')]:totals[target]+=cost[source]
        totals['maximum_memory']=max(totals['maximum_memory'],cost['memory_sites'])
        if ratio is None:totals['empty_source_events']+=1;continue
        poly=add(poly,scale(multiply(indicator(c|a,a),amp),-sign*ratio))
        totals['transition_groups']+=1
    if any(mask.bit_count()>4 for mask in poly):raise ValueError('At most degree-four row envelope required')
    shifts=[{0:F(-o.target),**{1<<i:F(1) for i in range(spin,o.modes,2)}} for spin in (0,1)]
    charge={0:F(-minimum),**{3<<(2*i):F(1) for i in range(o.sites)}}
    totals.update(polynomial_terms=len(poly),determinant_actions=len(o.oracle.cache),occupation_endpoint_evaluations=0)
    return o,poly,shifts,charge,totals


def replay(certificate):
    if certificate.get('kind')!='weighted_charge_tail_bracket_v1':raise ValueError('Unsupported tail-envelope bracket')
    o,f,shifts,charge,cost=envelope(certificate);modes=o.modes
    atoms=_atoms(certificate['positive_indicators'],modes,4)
    localizers=_atoms(certificate['charge_indicators'],modes,2)
    ideals=certificate['number_multipliers']
    if type(ideals) is not list or len(ideals)!=2:raise ValueError('Two number multipliers required')
    ideal={}
    for terms,shift in zip(ideals,shifts):
        if type(terms) is not list or len(terms)>30000:raise ValueError('Bounded number multiplier required')
        p={}
        for item in terms:
            mask=_mask(item['mask'],modes)
            if mask.bit_count()>3 or type(item['coefficient']) is not str:raise ValueError('Invalid number multiplier')
            p=add(p,{mask:F(item['coefficient'])})
        ideal=add(ideal,multiply(shift,p))
    if type(certificate['b']) is not str:raise ValueError('Exact polynomial bound required')
    b=F(certificate['b']);residual=add(f,{0:-b},scale(atoms,-1),scale(multiply(charge,localizers),-1),scale(ideal,-1))
    error=sum(map(abs,residual.values()),F(0));lower=b-error
    state=certificate['witness_state']
    if not o.oracle.valid_state(state) or sum(((state>>(2*i))&3)==3 for i in range(o.sites))<certificate['minimum_doublons']:
        raise ValueError('Physical fixed-spin tail witness required')
    upper=sum((value for mask,value in f.items() if state&mask==mask),F(0))
    if lower>upper:raise ValueError('Inconsistent polynomial bracket')
    return dict(cost,lower=str(lower),upper=str(upper),width=str(upper-lower),lower_float=float(lower),
                upper_float=float(upper),width_float=float(upper-lower),residual_l1=str(error),
                positive_indicators=len(certificate['positive_indicators']),charge_indicators=len(certificate['charge_indicators']),
                number_multiplier_terms=sum(map(len,ideals)),explicit_upper_witnesses=1,
                scope='Exact bracket on the minimum of a conservative weighted-row envelope over the fixed-spin ionic tail. Lower bound uses nonnegative degree4 indicators, (D-minimum) localizers, number identities and residual l1. Upper bound evaluates one physical witness for this polynomial only, not an energy or true weighted-row upper bound. Ratio maxima are exact finite-memory DP; distinct transitions still use independent maxima.')


def export_proposal(source, proposal, output, witness, minimum=2, denominator=10**12):
    source,proposal,output=Path(source),Path(proposal),Path(output)
    if output.exists():raise ValueError('Preserve previous tail bracket')
    if type(denominator) is not int or denominator<1:raise ValueError('Positive rationalization denominator required')
    data=json.loads(source.read_text());solution=json.loads(proposal.read_text())
    cert={k:data[k] for k in ['modes','particles','hamiltonian','metric_rule']}
    cert.update(kind='weighted_charge_tail_bracket_v1',minimum_doublons=minimum,witness_state=witness,
                b='0',positive_indicators=[],charge_indicators=[],number_multipliers=[[],[]])
    if len(solution['labels'])!=len(solution['values']):raise ValueError('Mismatched proposal labels')
    for label,value in zip(solution['labels'],solution['values']):
        rational=F(round(value*denominator),denominator)
        if label==['b']:cert['b']=str(rational);continue
        if not rational:continue
        if label[0]=='ideal':cert['number_multipliers'][label[1]].append({'mask':label[2],'coefficient':str(rational)})
        elif label[0] in ('positive','charge'):
            if rational<0:raise ValueError('Negative proposed positivity weight')
            key='positive_indicators' if label[0]=='positive' else 'charge_indicators'
            cert[key].append({'required':label[1],'occupied':label[2],'weight':str(rational)})
        else:raise ValueError('Unsupported LP proposal label')
    receipt=replay(cert);output.mkdir(parents=True)
    for name,value in [('certificate',cert),('receipt',receipt)]: (output/(name+'.json')).write_text(json.dumps(value,indent=2)+'\n')
    return receipt


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--verify',required=True);args=p.parse_args()
    print(json.dumps(replay(json.loads(Path(args.verify).read_text())),indent=2))
