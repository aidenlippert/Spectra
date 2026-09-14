**A response from a spectral interval**

Use H-b=[[A,B],[B*,D]], Q for a specified doubly occupied spatial orbital,
and P=I-Q. The response is X=p(D)B*, with no inverse in its construction.
For 0<delta<=D<=M, put c=(M+delta)/2, a=(M-delta)/2,
Z=(cI-D)/a, z0=c/a>1. The Chebyshev residual is

 r_k(D)=T_k(Z)/T_k(z0),    ||r_k(D)||<=1/T_k(z0).

Since r_k(0)=1, p_(k-1)(D)=(I-r_k(D))/D is a polynomial. The inverse
notation is removable by the divided-difference recurrence

 q_0=0, q_1=I/a,
 q_(j+1)=2 Z q_j-q_(j-1)+(2/a)T_j(z0)I,
 p_(k-1)=q_k/T_k(z0).

Applying this to B* takes k-1 applications of D and stores three working
vectors or batches. The program is compact; the work done by each action
and the number of right-hand sides remain explicit costs.

For the exact program E=B*-DX=r_k(D)B*. If ||B||<=g, then
||E||^2/delta<=eta=g^2/[delta T_k(z0)^2]. The congruence identity gives

 T*(H-b)T=[[K,E*],[E,D]],
 K=A-BX-X*B*+X*DX,   T=[[I,0],[-X,I]].

Thus K>=eta I suffices for H>=b. The interval and residual certificate
does not certify K. A rounded response matrix can instead be checked by
its actual residual; the same identity is valid for every X.

**A stronger molecular gap from hopping out of the fixed occupancy**

Write one retained real spatial density matrix as L=[[A,v],[v*,l]],
with the last orbital fixed doubly occupied. Define Q(A)=sum A_ij E_ij
on the remaining spatial modes. Exact CAR algebra gives

 Q Q(L)^2 Q = Q[(Q(A)+2l)^2+2||v||^2-Q(vv*)]Q.

The last two terms are the positive squared hopping amplitude out of Q.
They are present even when the within-Q square is dropped. For
R=(1/2)sum_i lambda_i Q(L_i)^2, define
V=sum_i lambda_i v_i v_i* and s0=sum_i lambda_i ||v_i||^2. Then

 QRQ >= Q[s0-(1/2)Q(V)]Q.

Combining this with c+2t_pp+Q(t_remaining) gives a one-body lower bound
with constant c+2t_pp+s0 and matrix t_remaining-V/2. For an even remaining
electron count, twice the sum of the lowest half as many diagonal entries,
minus N_remaining times the absolute off-diagonal row bound, is sufficient.
Add the certified lower tail shift before subtracting b.

For an upper endpoint, use the corresponding largest diagonal sum and row
bound for t_remaining. Each squared within-Q density has norm at most
(N_remaining*max(||A||_1,||A||_infinity)+2|l|)^2. The hopping square is at
most 2||v||^2. Sum the nonnegative weighted upper bounds and add the upper
tail shift. These derivations enumerate no fermionic states. They strengthen
the old rule that used only R>=0. Exact CAR checks test the hopping identity
independently on small rational density matrices.

For the coupling bound, expand only B*=QHP as a canonical operator
polynomial. Every CAR monomial has norm at most one, hence
||B||<=sum_w |coefficient_w(QHP)|. The full shifted-H norm is an independent
upper-bound control, but the compressed-density upper is substantially tighter.

**A charged full-certificate control**

For H6 only, explicitly check K on the retained states in conserved blocks.
N_up and a checked orbital-parity mask give a complete partition of the
N=6 space. Generate the 924 basis labels and count every Hamiltonian action.
This control is an enumerated retained solve, not a compact full solver.

Generate X by the same polynomial response recurrence and round its entries
to rational integers over a fixed denominator. For each conserved block,
compute E=B*-DX exactly and use eta_block=||E||_F^2/delta rounded upward.
The proposal supplies a rational triangular factor L for K-eta_block I.
Replay forms K-L L* -eta_block I exactly and requires its Gershgorin lower
bound to be nonnegative. All sectors, including uncoupled sectors, must pass.
This gives a complete H>=b proof, with explicit response-column and retained
factor costs. It requires no inverse but does expand the response actions.

Use the frozen upper U and b=U-0.001 Ha for this control, below the original
1.6 mHa target and the exact old-family interval floor. Keep the explicit
closure cost separate from the nonenumerating response certificate. Compare
the scalar-response closure on the same partition to isolate what the
polynomial response changes. Existing full-enumeration reference bounds must
also remain visible; crossing the old restricted-family floor alone does not
establish a competitive complete algorithm.
