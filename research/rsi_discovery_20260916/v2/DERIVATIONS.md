# Encoded claims and their limits

These arguments justify the new checking operations. They do not prove that an
arbitrary generated Python program is correct, that the molecular optimization
is universally tractable, or that Spectra improves indefinitely.

## Transporting a checked CAR template

For an input product, replace its distinct mode numbers by their ascending ranks.
The inverse substitution is injective and increasing. A CAR reduction uses only
the creation/annihilation flags, mode equality (the contraction delta), and mode
order (fermionic sorting signs). Those predicates are preserved by this
substitution. Every identity obtained by the CAR rewrite rules on rank labels
therefore transports to the original labels.

The generated program proposes a normal polynomial on rank labels. On each new
pattern, `CertifiedTemplates` checks that polynomial against the original exact
normalizer and separately checks literal ladder action on every basis vector of
its at-most-six-mode local algebra. The result stored in the cache is immutable.
A reuse substitutes labels into that already checked polynomial. It does not
trust another execution of the generated program. This local algebra check is
not construction or enumeration of the molecular fixed-particle-number sector.

`Transport.lean` checks equality-delta transport, order-inversion transport, and
the implication from strict order preservation to injectivity. It also checks
the stated integer cross-term and error-composition identities. Its accepted
axioms are only Lean's `propext`, `Quot.sound`, and `Classical.choice`. The full
Python checker and CAR representation theorem have not been formalized in Lean.

## Symmetric Gram maps

For real symmetric Q, the polynomial associated with a dictionary B is

\[
\sum_i Q_{ii}B_i^\dagger B_i+
\sum_{i<j}Q_{ij}(B_i^\dagger B_j+B_j^\dagger B_i).
\]

Every generated map must reproduce the coefficient of every canonical CAR word
for every upper-triangle coordinate. Both cross terms are mandatory. The real
solver adapter divides off-diagonal coefficients equally between the two full
matrix coordinates, then verifies exact equality to the symmetrized original
map. These half-integer entries are exactly representable in binary floating
point. Spin projection and the numerical optimization subsequently use the
existing Spectra implementation. Final energy acceptance reconstructs rational
certificates from the original Hamiltonian in a separate process without NumPy,
SciPy, CVXPY, PySCF or Quimb imports.

The pilot's generated second method groups the two terms by adjoint orbits at
the template level. A neutral canonical monomial with r creators and r
annihilators acquires reversal sign (-1)^(r(r-1))=1 under adjunction. Its adjoint
therefore exchanges the two ascending index groups. Coincident orbit members
are doubled; opposite coefficients cancel before concrete mode substitution.
The recorded implementation still passes the exact map checker on every new
input. This argument is not substituted for that acceptance gate.

## Physical nulls and optimization equivalence

### Polynomial-operator maps

The follow-up representation task uses P_i = sum_a C_ia B_a with exact integer
coefficients. For a monomial diagonal B_a†B_a, its coefficient in P_i†P_i is
C_ia², and in P_i†P_j+P_j†P_i is 2 C_ia C_ja. For a monomial symmetric cross
term B_a†B_b+B_b†B_a with a<b, the corresponding factors are C_ia C_ib and
C_ia C_jb+C_ib C_ja. These formulas include both cross terms and permit zero or
dependent polynomial operators. All arithmetic must stay exact even when
coefficients exceed binary floating-point precision.

The acceptance implementation expands the polynomial products through the
original exact CAR algebra. Candidate programs may choose their own construction,
including calls to the ordinary or acquired monomial map. All arms have the same
monomial-map interface. The additional common-backend transfer freezes the
candidate programs and replaces the primitive implementation with the ordinary
backend for every arm. A raw speedup that disappears there is insufficient for
this study's discovery-improvement gate. This is a finite scoped algorithm
comparison, not a new theorem about all polynomial maps or scientific discovery.

### Declared ideal witnesses

If K=T(N-n), then K annihilates every input state in the fixed-N sector. Thus
B and B'=B-K have the same action there. That physical identity alone does not
show that a particular truncated SOS optimization permits replacing B by B'.

For a proposed coordinate map B_i -> sum_a C_ia R_a, the new checker constructs
each diagonal and each Hermitian off-diagonal Gram difference. It supplies exact
rational coordinates expressing every difference in the **declared** free
columns. Only when every witness reconstructs exactly is that dictionary map
accepted as equivalent modulo those columns. Positive Gram matrices transport
as C^T Q C, which remains positive semidefinite because
v^T C^T Q C v=(Cv)^T Q(Cv)>=0. The checked free-column corrections account for
the full polynomial difference.

The free scalar energy coordinate is excluded. A nonzero scalar column, an
energy-labelled column, or a collection spanning the identity is refused.
The caller must separately justify that the supplied free columns vanish in the
stated sector; in the executed examples they are explicitly (N-n) times the
current permitted Hermitian number multipliers. The routine does not authorize
an arbitrary newly supplied physical constraint.

The executed examples reduce three declared operators to two at (m,N)=(6,2),
(8,4), and (10,4), with all six upper-triangle witnesses checked. A narrower
declared ideal fails the equivalence gate while the physical null remains
valid. This is an obstruction to this exact representation witness, not a
proof that two entire differently parameterized SOS relaxation optima differ.

## Controlled observable dynamics

Let O be the fixed Hermitian observable basis and write exact coefficient
identities

\[
i[H_r,O_i]=\sum_j (M_r)_{ij}O_j+R_{r,i}.
\]

Each CAR monomial has operator norm at most one. For Gaussian rational
coefficients, the sum of absolute real and imaginary parts therefore bounds
the operator norm of each residual. Under the original controls u,v, this gives
a componentwise defect bound epsilon_i(H)+|u|epsilon_i(D)+|v|epsilon_i(W).
This is a bound on the defining differential equation, not agreement between
two numerical trajectories.

For exact moments x and the reduced solution z with the same initial moments,
the infinity-norm error over a constant-control interval of length t satisfies

\[
\|x(t)-z(t)\|_\infty\le e^{Lt}\|x(0)-z(0)\|_\infty+
\delta(e^{Lt}-1)/L,
\]

where L=||M||_infinity and delta is the maximum residual bound. At L=0 the
second term is delta*t. This follows by variation of constants and
||exp(Mt)||_infinity<=exp(Lt). Interval errors compose across all four phases.

The reduced propagator is a rational Taylor polynomial. For x>=0 and N+2>x,
the exponential tail after order N is bounded above by

\[
\frac{x^{N+1}}{(N+1)!}\frac{1}{1-x/(N+2)}.
\]

Every later term ratio is at most x/(N+2). The same bound with x=||M||t
controls the matrix-exponential truncation. Rational rounding after each phase
is explicitly added to the propagated error; exponential/error bounds round
outward. Inputs outside the stated Taylor envelope are refused.

The original two actuator norms are at most two, and the target population
contrast D has norm at most two and spectrum in [-2,2]. The original uncertainty
budgets add 16*T*0.001 + 4*0.0005 + 4*0.001 to the final observable radius,
using the Duhamel trace-norm bound, the stated initial trace distance, and the
total integrated phase-flip rate. The original four half-unit-duration phases,
amplitude limit 1/2, Hamiltonian and MPS remain fixed.

The four-observable experiment returns [-2,2], so it **does not** prove D<=-3/5.
The conditional branch retains the exact matrices, residual polynomials,
initial moments and remaining obligation. It reads no full-state trajectory or
enumerated reduced embedding. This failure is not a proof of unreachability.

The executable frontier extension allows candidate programs to add up to
fourteen observables to the fixed identity and D. Each is checked for exact
Gaussian-rational coefficients, particle neutrality, Hermiticity, original mode
range, and degree at most four. The entire basis must be independent. Candidate
code only constructs an observable artifact; a separate standard-library-only
process derives the dynamics, initial moments, residual bounds and final
enclosure from the original inputs. The same error argument applies to any
admitted basis. Large commutator closures and exponential envelopes are resource
obstructions, not proofs that no useful reduction exists. Conditional branches
are retained but cannot authorize later acceptance.
