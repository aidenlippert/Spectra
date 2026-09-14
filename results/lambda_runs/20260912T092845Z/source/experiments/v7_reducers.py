"""Bounded conventional sparse candidate generators."""
from fractions import Fraction
from collections import deque
from .v7_certificate import Piece
MAX_ORDER,MAX_BASIS=24,1024
def _v(x): return {k:Fraction(a) for k,a in x.items() if Fraction(a)}
def _cost(g,**k):
 for a,b in getattr(g,'cost',{}).items(): k[a]=k.get(a,0)+b
 return k
def _check(n):
 if n<0 or n>MAX_ORDER: raise ValueError('order must be in [0,24]')
def full_taylor(gen,initial,duration,order):
 _check(order); c=_v(initial); a=[c]
 for j in range(order): c={p:x/Fraction(j+1) for p,x in gen.apply(c).items()}; a.append(c)
 return Piece(Fraction(duration),tuple(a)),_cost(gen,trials=0,sorts=0,solve=0,orthogonalization=0,memory=sum(map(len,a)))
def projected_taylor(gen,initial,duration,order,basis):
 _check(order); B=set(basis)
 if len(B)>MAX_BASIS: raise ValueError('basis cap exceeded')
 c=_v(initial); a=[c]
 for j in range(order): c={p:x/Fraction(j+1) for p,x in gen.apply(c).items() if p in B}; a.append(c)
 return Piece(Fraction(duration),tuple(a)),_cost(gen,trials=0,sorts=0,solve=0,orthogonalization=0,memory=sum(map(len,a)))
def bfs_basis(gen,initial,depth,cap=MAX_BASIS):
 if cap<1 or cap>MAX_BASIS: raise ValueError('basis cap exceeded')
 B=set(initial); q=deque((p,0) for p in B); t=0
 while q:
  p,d=q.popleft(); t+=1
  if d>=depth: continue
  for z in gen.apply({p:Fraction(1)}):
   if z not in B:
    if len(B)>=cap: raise RuntimeError('basis cap exceeded')
    B.add(z); q.append((z,d+1))
 return B,_cost(gen,trials=t,sorts=0,solve=0,orthogonalization=0,memory=len(B))
def residual_basis(gen,initial,duration,order,cap=MAX_BASIS,batch=1):
 _check(order); B=set(initial); c=_v(initial); a=[c]; t=0
 for j in range(order):
  raw=gen.apply(c); om={p:x for p,x in raw.items() if p not in B}; t+=len(om)
  for p in sorted(om,key=lambda p:(-abs(om[p]),p))[:max(1,batch)]:
   if len(B)>=cap: raise RuntimeError('basis cap exceeded')
   B.add(p)
  c={p:x/Fraction(j+1) for p,x in raw.items() if p in B}; a.append(c)
 return Piece(Fraction(duration),tuple(a)),_cost(gen,trials=t,sorts=order,solve=0,orthogonalization=0,memory=sum(map(len,a)))
def arnoldi(gen,initial,duration,order,krylov_dim):
 """Two-pass modified Gram-Schmidt and projected exponential polynomial.

 Numeric Krylov calculations propose an exact rational polynomial; the trusted
 checker recomputes its full residual, including numerical inaccuracies.
 """
 _check(order)
 if type(krylov_dim) is not int or not 1<=krylov_dim<=64:
  raise ValueError('krylov_dim must be 1..64')
 s=_v(initial); beta=sum(float(x)**2 for x in s.values())**.5
 if not beta: raise ValueError('nonzero initial observable required')
 V=[{p:float(x)/beta for p,x in s.items()}]; columns=[]; orth=0; solve=0
 for j in range(krylov_dim):
  q={p:Fraction(str(x)) for p,x in V[j].items()}
  w={p:float(x) for p,x in gen.apply(q).items()}
  h=[0.0]*(j+2)
  for repeat in range(2):
   for i,v in enumerate(V):
    dot=sum(w.get(p,0.0)*x for p,x in v.items()); orth+=len(v)
    h[i]+=dot
    for p,x in v.items(): w[p]=w.get(p,0.0)-dot*x; orth+=1
  norm=sum(x*x for x in w.values())**.5; orth+=len(w)
  h[j+1]=norm; columns.append(h)
  if norm<=1e-13 or j+1==krylov_dim: break
  V.append({p:x/norm for p,x in w.items() if x})
 m=len(V); labels=sorted(set().union(*(v.keys() for v in V)))
 if len(labels)>gen.max_terms: raise RuntimeError('Krylov support cap')
 H=[[columns[j][i] if i<len(columns[j]) else 0.0 for j in range(m)] for i in range(m)]
 y=[beta]+[0.0]*(m-1); out=[]
 for k in range(order+1):
  c={}
  for p in labels:
   x=sum(y[i]*V[i].get(p,0.0) for i in range(m)); solve+=m
   if x: c[p]=Fraction(str(x))
  out.append(c)
  if k<order:
   y=[sum(H[i][j]*y[j] for j in range(m))/(k+1) for i in range(m)]; solve+=m*m
 if sum(map(len,out))>50000: raise RuntimeError('total coefficient cap')
 return Piece(Fraction(duration),tuple(out)),_cost(gen,trials=0,sorts=1,solve=solve,orthogonalization=orth,memory=sum(map(len,out))+sum(map(len,V))+m*m)

generate_full_taylor=full_taylor; generate_projected_taylor=projected_taylor; generate_residual_basis=residual_basis; generate_arnoldi=arnoldi
