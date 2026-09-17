"""Compile exact action inner products once for many trajectory queries.

The constructor still visits every supplied/reached configuration. A kernel
loaded from disk is not trusted: this API constructs it afresh from original
problem data and the proposal. Later queries use only its small exact matrices.
"""
from fractions import Fraction as F
from math import lcm
import time
from research.intervention_reduction_20260916.exact import digest,exact_moments


def quadratic(gram,z):
    applied=[(sum(g*a for g,(a,b) in zip(row,z)),sum(g*b for g,(a,b) in zip(row,z))) for row in gram]
    return sum(a*x+b*y for (a,b),(x,y) in zip(z,applied))


class CheckedKernel:
    """Only original-input construction; no deserialization accepting shortcut."""
    def __init__(self,data,proposal,orbital):
        started=time.monotonic()
        if type(orbital) is not int or not 0<=2*orbital+1<data['modes']:
            raise ValueError('Observed spatial orbital outside model')
        moments=exact_moments(data,proposal)
        rank=moments['rank'];rows=dict(zip(proposal['configurations'],proposal['vectors']))
        actions=[rows]+moments['integer_actions'];zero=[0]*rank
        labels=sorted(set().union(*(set(a) for a in actions)))
        joint=[sum((a.get(s,zero) for a in actions),[]) for s in labels]
        columns=list(zip(*joint));width=len(columns)
        gram=[[0]*width for _ in range(width)]
        for i,column in enumerate(columns):
            for j in range(i,width):
                gram[i][j]=gram[j][i]=sum(a*b for a,b in zip(column,columns[j]))
        occupation=lambda s: ((s>>(2*orbital))&1)+((s>>(2*orbital+1))&1)
        self.observable=[[sum(occupation(s)*row[i]*row[j] for s,row in rows.items())
                          for j in range(rank)] for i in range(rank)]
        self.metric=[row[:rank] for row in gram[:rank]]
        self.action_gram=gram;self.action_denominators=moments['action_denominators']
        # No configuration vector or action is retained in the reusable object.
        self.moments={k:moments[k] for k in ('metric','rank','amplitudes','counts')}
        self.fixture_sha256=digest(data);self.proposal_sha256=digest(proposal);self.orbital=orbital
        self.modes=data['modes'];self.vector_denominator=proposal['denominator']
        self.query_data={'modes':data['modes']};self.query_proposal={'controls':proposal['controls']}
        self.receipt={'fixture_sha256':self.fixture_sha256,'proposal_sha256':self.proposal_sha256,
            'orbital':orbital,'reduced_dimension':rank,'joint_action_dimension':width,
            'stored_integer_matrix_entries':width**2+rank**2,
            'joint_action_gram_sha256':digest(gram),'observable_gram_sha256':digest(self.observable),
            'construction_seconds':time.monotonic()-started,
            'configuration_records_retained_for_queries':0,**moments['counts']}

    def norm_square(self,z):
        return quadratic(self.metric,z)

    def observable_numerator(self,z):
        return quadratic(self.observable,z)

    def residual_integral(self,controls,shift,coefficients,duration,coefficient_den,initial_norm):
        weights=[F(1)]+controls;dens=self.action_denominators;rank=self.moments['rank']
        common=lcm(*(d*w.denominator for d,w in zip(dens,weights)),shift.denominator)
        scales=[int(common*w/d) for d,w in zip(dens,weights)]
        derivative_scale=common*duration.denominator
        degree=len(coefficients)-1;residual=[]
        for k,c in enumerate(coefficients):
            derivative=coefficients[k+1] if k<degree else [[0,0]]*rank
            vector=[[-derivative_scale*(k+1)*di+duration.numerator*int(common*shift)*a,
                      derivative_scale*(k+1)*dr+duration.numerator*int(common*shift)*b]
                    for (a,b),(dr,di) in zip(c,derivative)]
            for scale in scales:
                vector.extend([[-duration.numerator*scale*a,-duration.numerator*scale*b] for a,b in c])
            residual.append(vector)
        applied=[[(sum(g*a for g,(a,b) in zip(row,z)),sum(g*b for g,(a,b) in zip(row,z)))
                  for row in self.action_gram] for z in residual]
        integral_den=lcm(*range(1,2*degree+2));numerator=0
        for i,z in enumerate(residual):
            for j in range(i,len(residual)):
                dot=sum(a*x+b*y for (a,b),(x,y) in zip(z,applied[j]))
                numerator+=(1 if i==j else 2)*dot*(integral_den//(i+j+1))
        if numerator<0:
            raise AssertionError('Negative exact compiled residual square')
        return F(numerator,integral_den*derivative_scale**2*coefficient_den**2*initial_norm)
