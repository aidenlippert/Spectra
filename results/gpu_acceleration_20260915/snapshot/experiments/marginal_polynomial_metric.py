"""Joint denominator-free occupation polynomials for a rational charge metric."""
from fractions import Fraction as F
from math import lcm, comb
from itertools import combinations
from experiments.marginal_coherent_tree import CoherentCharge
from experiments.marginal_charge_product import feature_orbits


class JointPolynomial:
    def __init__(self,data):
        modes=data['modes'];sites=modes//2
        # CoherentCharge supplies checked CAR groups and signs; its unit
        # auxiliary metric is not used in the polynomial metric inequality.
        unit={'family':'pair','max_pair_distance':1,'factors':['1']*len(feature_orbits(sites,'pair',1))}
        self.o=CoherentCharge(dict(data,metric_rule=unit));self.modes=modes;self.sites=sites
        self.spin_masks=self.o.spin_masks;self.target=self.o.target
        metric=data['polynomial_metric'];self.denominator=metric['denominator']
        if type(self.denominator) is not int or not 1<=self.denominator<=10**14:raise ValueError('Bounded positive metric denominator required')
        if type(metric['terms']) is not list or not 1<=len(metric['terms'])<=1000:raise ValueError('Bounded charge polynomial required')
        self.terms=[]
        for term in metric['terms']:
            powers=term['powers'];coefficient=term['coefficient']
            if type(powers) is not list or len(powers)!=sites or any(type(p) is not int or not 0<=p<=2 for p in powers) or sum(powers)>6 or type(coefficient) is not int:
                raise ValueError('Integer charge metric with degree at most six required')
            self.terms.append((tuple(powers),coefficient))
        self.zero_on_valence=sum(c for powers,c in self.terms if not any(powers))==0
        self.stats={'polynomial_products':0,'coefficient_products':0,'peak_terms':0,'occupation_endpoint_evaluations':0}

    def add(self,*polys):
        out={}
        for poly in polys:
            for mask,value in poly.items():out[mask]=out.get(mask,0)+value
        return {m:v for m,v in out.items() if v}

    def scale(self,poly,value):return {m:v*value for m,v in poly.items() if v*value}

    def multiply(self,a,b):
        out={};self.stats['polynomial_products']+=1
        for x,u in a.items():
            for y,v in b.items():
                mask=x|y;self.stats['coefficient_products']+=1
                if any((mask&s).bit_count()>self.target for s in self.spin_masks):continue
                out[mask]=out.get(mask,0)+u*v
        out={m:v for m,v in out.items() if v};self.stats['peak_terms']=max(self.stats['peak_terms'],len(out));return out

    def weight(self,charges):
        total=0
        for powers,coefficient in self.terms:
            value=coefficient
            for q,power in zip(charges,powers):value*=q**power
            total+=value
        return F(total,self.denominator)

    def metric_polynomial(self,q):
        powers=[[{0:1},x,self.multiply(x,x)] for x in q];result={}
        for exponents,coefficient in self.terms:
            term={0:coefficient}
            for i,power in enumerate(exponents):
                if power:term=self.multiply(term,powers[i][power])
            result=self.add(result,term)
        return result

    def compile(self,gamma):
        gamma=F(gamma);o=self.o
        hden=lcm(gamma.denominator,*(x.denominator for x in o.oracle.h.values()))
        n=[{1<<i:1} for i in range(self.modes)]
        q=[self.add({0:-1},n[2*i],n[2*i+1]) for i in range(self.sites)]
        weight=self.metric_polynomial(q)
        diagonal={m:int(x*hden) for m,x in o.oracle.diagonal.items()}
        result=self.multiply(self.add(diagonal,{0:-int(gamma*hden)}),weight)
        used=0
        for (c,a),amplitude in o.groups.items():
            closed=o.close(c|a,a)
            if closed is None:continue
            lo,hi=o.amplitude_range(amplitude,*closed)
            if lo>=0:sign=1
            elif hi<=0:sign=-1
            else:raise ValueError('Joint polynomial requires fixed-sign grouped amplitudes')
            source_n=[({0:1} if a&(1<<i) else {}) if (c|a)&(1<<i) else n[i] for i in range(self.modes)]
            amp={}
            for support,x in amplitude.items():amp=self.add(amp,self.scale(source_n[support.bit_length()-1] if support else {0:1},int(x*hden)))
            if not amp:continue
            delta=tuple(((c>>(2*i))&3).bit_count()-((a>>(2*i))&3).bit_count() for i in range(self.sites))
            changed=[self.add({0:delta[i]-1},source_n[2*i],source_n[2*i+1]) for i in range(self.sites)]
            target_weight=self.metric_polynomial(changed)
            # On a valence target q'=0. If v(0)=0, v(q')Q(q')=v(q')
            # throughout the physical sector, so no global projector is needed.
            if not self.zero_on_valence:
                target_p={0:1}
                for charge in changed:target_p=self.multiply(target_p,self.add({0:1},self.scale(self.multiply(charge,charge),-1)))
                target_q=self.add({0:1},self.scale(target_p,-1))
            event={0:1}
            for i in range(self.modes):
                if c&(1<<i):event=self.multiply(event,{0:1,1<<i:-1})
                elif a&(1<<i):event=self.multiply(event,n[i])
            penalty=self.multiply(self.multiply(event,amp),target_weight)
            if not self.zero_on_valence:penalty=self.multiply(penalty,target_q)
            result=self.add(result,self.scale(penalty,-sign));used+=1
        receipt=dict(self.stats,metric_scale=self.denominator,numerator_scale=self.denominator*hden,
            weight_terms=len(weight),numerator_terms=len(result),numerator_degree=max((m.bit_count() for m in result),default=0),
            transition_groups=used,determinant_actions=len(o.oracle.cache),target_valence_exclusion=True)
        if self.zero_on_valence:receipt['projector_free_zero_metric']=True
        return weight,result,receipt


def check_positive(ring,poly,scale,proof):
    """Check a bounded Boolean/localizer decomposition with exact residual."""
    denominator=proof['denominator'];bound=proof['bound'];degree=6
    if type(denominator) is not int or not 1<=denominator<=10**16 or type(bound) is not int:raise ValueError('Exact bounded positivity scale and bound required')
    def mask(value):
        if type(value) is not int or not 0<=value<1<<ring.modes:raise ValueError('Invalid occupation mask')
        return value
    def atoms(raw,maximum_degree):
        if type(raw) is not list or len(raw)>100000:raise ValueError('Bounded positive atom list required')
        result={}
        for item in raw:
            required,occupied=mask(item['required']),mask(item['occupied']);weight=item['weight']
            if occupied&~required or required.bit_count()>maximum_degree or type(weight) is not int or weight<0:raise ValueError('Nonnegative bounded indicator required')
            term={occupied:weight}
            if any((occupied&s).bit_count()>ring.target for s in ring.spin_masks):continue
            for i in range(ring.modes):
                if (required^occupied)&(1<<i):term=ring.multiply(term,{0:1,1<<i:-1})
            result=ring.add(result,term)
        return result
    positive=atoms(proof['positive_indicators'],degree)
    localizer=atoms(proof['charge_indicators'],degree-2)
    charge={0:-1,**{3<<(2*i):1 for i in range(ring.sites)}}
    represented=ring.add({0:bound},positive,ring.multiply(charge,localizer))
    ideals=proof['number_multipliers']
    if type(ideals) is not list or len(ideals)!=2:raise ValueError('Two number identities required')
    for spin,terms in enumerate(ideals):
        if type(terms) is not list or len(terms)>100000:raise ValueError('Bounded number multiplier required')
        multiplier={}
        for item in terms:
            support=mask(item['mask']);value=item['coefficient']
            if support.bit_count()>degree-1 or type(value) is not int:raise ValueError('Integer bounded number multiplier required')
            if any((support&s).bit_count()>ring.target for s in ring.spin_masks):continue
            multiplier=ring.add(multiplier,{support:value})
        shift={0:-ring.target,**{1<<i:1 for i in range(spin,ring.modes,2)}}
        represented=ring.add(represented,ring.multiply(shift,multiplier))
    residual=ring.add(ring.scale(poly,denominator),ring.scale(represented,-scale))
    error=F(sum(abs(x) for x in residual.values()),scale*denominator)
    return {'lower':str(F(bound,denominator)-error),'residual_l1':str(error),
            'positive_indicators':len(proof['positive_indicators']),'charge_indicators':len(proof['charge_indicators']),
            'number_multiplier_terms':sum(map(len,ideals))}


def complete_number_ideals(ring,poly):
    """Proposer-only exact quotient lift; replay still checks the full identity.

    x_S = (sum_{i outside S} x_{S+i} - (N_sigma-target)x_S)/(target-|S|).
    The remainder consists of full-population monomials, so this construction
    may have sector-sized cost. It is not a compressed discovery algorithm.
    """
    current={mask:F(value) for mask,value in poly.items() if value};ideals=[{},{}]
    for index,spin in enumerate(ring.spin_masks):
        for degree in range(ring.target):
            selected=[mask for mask in current if (mask&spin).bit_count()==degree]
            for mask in selected:
                coefficient=current.pop(mask)/F(ring.target-degree)
                ideals[index][mask]=ideals[index].get(mask,F(0))-coefficient
                for i in range(ring.modes):
                    bit=1<<i
                    if spin&bit and not mask&bit:
                        target=mask|bit;current[target]=current.get(target,F(0))+coefficient
        current={mask:value for mask,value in current.items() if value}
    return current,[{mask:value for mask,value in ideal.items() if value} for ideal in ideals]


def replay(certificate):
    if certificate.get('kind')!='joint_polynomial_metric_gap_v1':raise ValueError('Unsupported joint polynomial gap')
    if type(certificate.get('target_lower')) is not str:raise ValueError('Exact gap target required')
    ring=JointPolynomial(certificate);u,k,compilation=ring.compile(F(certificate['target_lower']))
    weight=check_positive(ring,u,compilation['metric_scale'],certificate['weight_proof'])
    numerator=check_positive(ring,k,compilation['numerator_scale'],certificate['numerator_proof'])
    if F(weight['lower'])<=0 or F(numerator['lower'])<0:raise ValueError('Strict metric positivity and nonnegative joint numerator required')
    return {'complement_lower':str(F(certificate['target_lower'])),'weight_positivity':weight,
            'numerator_positivity':numerator,'compilation':compilation,
            'verification_coefficient_products':ring.stats['coefficient_products'],
            'scope':'Exact joint weighted-Q-row certificate: positive rational charge polynomial metric, exact target-valence exclusion, coherent signed CAR amplitudes, and global Boolean/localizer positivity decompositions on fixed spin populations and D>=1. No charge-pattern or occupation endpoint iteration in replay. Polynomial/atom counts may still scale combinatorially; no general representability claim.'}


def check_nonsingleton_separator(ring,poly,scale,proof):
    """Finite dual obstruction, including the positivity checker's L1 allowance."""
    if ring.modes>12:raise ValueError('Finite separator supports at most twelve modes')
    states,weights=proof['states'],proof['weights']
    if (type(states) is not list or type(weights) is not list or len(states)!=len(weights)
            or not 1<=len(states)<=400):raise ValueError('Bounded signed occupation functional required')
    seen=set();fractions=[]
    for state,weight in zip(states,weights):
        if (type(state) is not int or state in seen or not 0<=state<1<<ring.modes
                or any((state&s).bit_count()!=ring.target for s in ring.spin_masks)
                or type(weight) is not str):raise ValueError('Unique physical fixed-spin states and exact weights required')
        seen.add(state);fractions.append(F(weight))
    denominator=lcm(*(v.denominator for v in fractions));moments={};charge_moments={}
    for state,value in zip(states,fractions):
        weight=int(value*denominator);charge=sum(((state>>(2*i))&3)==3 for i in range(ring.sites))-1
        mask=state
        while True:
            moments[mask]=moments.get(mask,0)+weight
            charge_moments[mask]=charge_moments.get(mask,0)+weight*charge
            if not mask:break
            mask=(mask-1)&state
    if moments.get(0)!=denominator:raise ValueError('Functional must map one to one')
    maximum=max(map(abs,moments.values()))
    if maximum>denominator:raise ValueError('Monomial moments must cover the L1 residual allowance')
    def indicator(table,required,occupied):
        empty=required^occupied;subset=empty;value=0
        while True:
            value+=(-1 if subset.bit_count()%2 else 1)*table.get(occupied|subset,0)
            if not subset:return value
            subset=(subset-1)&empty
    positive=localizer=0;minimum=0
    for degree in range(min(6,ring.modes)+1):
        for indices in combinations(range(ring.modes),degree):
            required=sum(1<<i for i in indices);occupied=required
            while True:
                completions=1
                for spin in ring.spin_masks:
                    remaining=ring.sites-(required&spin).bit_count();needed=ring.target-(occupied&spin).bit_count()
                    completions*=comb(remaining,needed) if 0<=needed<=remaining else 0
                if completions>1:
                    value=indicator(moments,required,occupied);positive+=1;minimum=min(minimum,value)
                if degree<=4 and completions:
                    value=indicator(charge_moments,required,occupied);localizer+=1;minimum=min(minimum,value)
                if not occupied:break
                occupied=(occupied-1)&required
    if minimum<0:raise ValueError('Functional is negative on an admitted atom')
    value=F(sum(v*moments.get(m,0) for m,v in poly.items()),scale*denominator)
    if value>=0:raise ValueError('Strictly negative separating value required')
    return {'functional_value':str(value),'functional_value_float':float(value),'normalization':'1',
            'minimum_atom':str(F(minimum,denominator)),'max_abs_monomial_moment':str(F(maximum,denominator)),
            'positive_indicators_checked':positive,'charge_indicators_checked':localizer,
            'functional_states':len(states),'residual_l1_allowance_covered':True}


def replay_separator(certificate):
    if certificate.get('kind')!='joint_polynomial_nonsingleton_separator_v1':raise ValueError('Unsupported polynomial separator')
    if type(certificate.get('target_lower')) is not str:raise ValueError('Exact gap target required')
    ring=JointPolynomial(certificate)
    if ring.modes>12:raise ValueError('Finite separator supports at most twelve modes')
    _,poly,cost=ring.compile(F(certificate['target_lower']))
    result=check_nonsingleton_separator(ring,poly,cost['numerator_scale'],certificate['functional'])
    return dict(result,target_lower=certificate['target_lower'],compilation=cost,
        scope='Exact finite obstruction for this metric and gamma: no nonnegative numerator lower bound from nonsingleton occupation indicators of degree<=6, (D-1) indicators of degree<=4 and fixed-spin number identities, including the coefficient-L1 residual allowance. Signed functional is not a physical probability or energy witness. Other metrics, localizers and SOS/operator cones remain unexcluded.')


def export_proposal(candidate,weight_solution,numerator_solution,output,denominator=10**12):
    import json
    from pathlib import Path
    output=Path(output)
    if output.exists():raise ValueError('Preserve previous polynomial gap')
    if type(denominator) is not int or not 1<=denominator<=10**16:raise ValueError('Bounded positive rationalization denominator required')
    certificate=json.loads(Path(candidate).read_text());certificate['kind']='joint_polynomial_metric_gap_v1'
    for key,path in [('weight_proof',weight_solution),('numerator_proof',numerator_solution)]:
        proposal=json.loads(Path(path).read_text())
        if len(proposal['labels'])!=len(proposal['values']):raise ValueError('Proposal labels and values must match')
        proof={'denominator':denominator,'bound':0,'positive_indicators':[],'charge_indicators':[],'number_multipliers':[[],[]]}
        for label,value in zip(proposal['labels'],proposal['values']):
            integer=round(value*denominator)
            if label==['b']:proof['bound']=integer;continue
            if not integer:continue
            if label[0]=='ideal':proof['number_multipliers'][label[1]].append({'mask':label[2],'coefficient':integer})
            elif label[0] in ('positive','charge'):
                if integer<0:raise ValueError('Negative positivity proposal')
                proof['positive_indicators' if label[0]=='positive' else 'charge_indicators'].append({'required':label[1],'occupied':label[2],'weight':integer})
            else:raise ValueError('Unknown proposal column')
        certificate[key]=proof
    receipt=replay(certificate);output.mkdir(parents=True)
    for name,value in [('certificate',certificate),('receipt',receipt)]: (output/(name+'.json')).write_text(json.dumps(value,indent=2)+'\n')
    return receipt


if __name__=='__main__':
    import argparse,json
    from pathlib import Path
    p=argparse.ArgumentParser();p.add_argument('--verify',required=True);args=p.parse_args()
    certificate=json.loads(Path(args.verify).read_text())
    verify=replay_separator if certificate.get('kind')=='joint_polynomial_nonsingleton_separator_v1' else replay
    print(json.dumps(verify(certificate),indent=2))
