"""Degree-aware density-conditioned factor closure.

This keeps the Hastings-style factor symbolic: F=psi_i+n_j D+theta and
theta=tau+mu n_j+epsilon n_j nu.  Multiplication records every monomial,
including the density-times-cubic terms (degree six after squaring). Dropping
those terms is rejected by the checker.
"""
from fractions import Fraction as F
raise RuntimeError('Rejected exploratory prototype: use density_conditioned.exact_density for the actual paired CAR closure.')
from experiments.marginal_symbolic import product,canonical,adj,mono,add

def factor(psi_i, n_j, D, tau, mu, eps, nu):
    # arguments are tuples representing normal-ordered symbolic words
    return add(mono(psi_i,F(1)), add(*(mono(((n_j,)+w),c) for w,c in D.items())), mono((n_j,),F(mu)), mono((),F(tau)), add(*(mono(((n_j,)+w),eps*c) for w,c in nu.items())))

def closure(f): return canonical(product(adj(f),f))

def check_all_terms(poly, dropped=False):
    if dropped and any(len(w)>=6 for w in poly): raise ValueError('cannot drop density-times-cubic closure terms')
    return {'terms':len(poly),'max_degree':max(map(len,poly),default=0),'degree6_terms':sum(len(w)>=6 for w in poly)}
