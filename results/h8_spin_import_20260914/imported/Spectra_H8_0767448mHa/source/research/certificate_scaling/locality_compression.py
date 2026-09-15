"""Bounded locality experiment for fermionic cluster certificate candidates.

The certificate is deliberately conservative.  Split a spinful Hubbard chain
into contiguous clusters, solve each cluster exactly in every particle sector,
and add a triangle-inequality penalty for the omitted boundary hoppings.  This
gives a certified lower bound on the full-chain ground energy.  It is not a
claim that the penalty is sharp; the experiment measures whether locality
compresses certificates at all under an explicit, countable budget.
"""
from __future__ import annotations

import argparse, json, math, time
from fractions import Fraction
from pathlib import Path
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh


def basis(norb: int, n: int):
    return [x for x in range(1 << norb) if x.bit_count() == n]


def apply_annihilate(x: int, p: int):
    if not (x >> p) & 1: return None
    return x ^ (1 << p), (-1) ** ((x & ((1 << p) - 1)).bit_count())


def hop(x: int, dst: int, src: int):
    a = apply_annihilate(x, src)
    if a is None: return None
    y, s1 = a
    if (y >> dst) & 1: return None
    s2 = (-1) ** ((y & ((1 << dst) - 1)).bit_count())
    return y | (1 << dst), s1 * s2


def hubbard(nsite: int, n: int, t: float, u: float, periodic: bool = False):
    """Return sparse spinful Hubbard Hamiltonian in a fixed-N sector."""
    states = basis(2 * nsite, n); idx = {x:i for i,x in enumerate(states)}
    rows, cols, vals = [], [], []
    for col, x in enumerate(states):
        diag = u * sum(((x >> (2*i)) & 1) and ((x >> (2*i+1)) & 1)
                       for i in range(nsite))
        rows.append(col); cols.append(col); vals.append(float(diag))
        for i in range(nsite - (0 if periodic else 1)):
            j = (i + 1) % nsite
            for spin in (0, 1):
                p, q = 2*i + spin, 2*j + spin
                for dst, src in ((q,p),(p,q)):
                    out = hop(x, dst, src)
                    if out:
                        y, sign = out
                        rows.append(idx[y]); cols.append(col); vals.append(-t * sign)
    return coo_matrix((vals, (rows, cols)), shape=(len(states), len(states))).tocsr()


def ground(nsite: int, n: int, t: float, u: float, backend="numpy"):
    if n < 0 or n > 2*nsite: return math.inf
    h = hubbard(nsite, n, t, u)
    if h.shape[0] == 1: return float(h[0,0])
    if backend == "cupy":
        try:
            import cupy as cp
            # Cluster sectors are intentionally small; dense GPU batches are
            # useful for hundreds of parameter cases, while global truth stays CPU.
            return float(cp.linalg.eigvalsh(cp.asarray(h.toarray()))[0].get())
        except ImportError:
            raise RuntimeError("--backend cupy requested but CuPy is unavailable")
    return float(eigsh(h, k=1, which="SA", return_eigenvectors=False, tol=1e-11)[0])


def exact_lower_candidate(nsite, n, t, u, cache):
    """Rational LDL replay of H-ell I for small local sectors."""
    if int(t) != t or int(u) != u:
        raise ValueError("--exact requires integer t and u for rational replay")
    key=(nsite,n,int(t),int(u))
    if key in cache: return cache[key]
    h=np.asarray(hubbard(nsite,n,t,u).toarray(), dtype=float)
    if not np.allclose(h, h.T):
        raise ValueError("exact replay refused: local Hamiltonian is not symmetric")
    ev=float(np.linalg.eigvalsh(h)[0])
    # Downward decimal rounding is only a candidate; LDL is the proof gate.
    ell=Fraction(math.floor((ev-1e-7)*10**8),10**8)
    states=basis(2*nsite,n); m=len(states)
    for _ in range(8):
        a=[[Fraction(int(round(h[i,j])))-ell*(i==j) for j in range(m)] for i in range(m)]
        good=True; d=[]; L=[[Fraction(int(i==j)) for j in range(m)] for i in range(m)]
        for i in range(m):
            piv=a[i][i]-sum(L[i][k]*L[i][k]*d[k] for k in range(i))
            if piv <= 0: good=False; break
            d.append(piv)
            for j in range(i+1,m):
                L[j][i]=(a[j][i]-sum(L[j][k]*d[k]*L[i][k] for k in range(i)))/piv
        if good: break
        ell -= Fraction(1,10**8)
    else:
        raise RuntimeError("fractional LDL replay failed after 8 downward candidates")
    cache[key]=(ell,str(ell),m)
    return cache[key]


def cluster_lower(nsite: int, n: int, max_cluster: int, t: float, u: float, backend="numpy", exact=False, exact_cache=None):
    """Decoupled-cluster minimum minus explicit omitted boundary norm."""
    parts = [min(max_cluster, nsite-i) for i in range(0, nsite, max_cluster)]
    exact_mode = exact
    dp = ([Fraction(0)] + [None] * n) if exact_mode else ([0.0] + [math.inf] * n)
    cache = {}
    for size in parts:
        vals = [exact_lower_candidate(size,q,t,u,exact_cache)[0] if exact else ground(size, q, t, u, backend) for q in range(2*size+1)]
        cache[size] = vals
        ndp = ([None] * (n+1)) if exact_mode else ([math.inf] * (n+1))
        for used, old in enumerate(dp):
            if old is None or (not exact_mode and not math.isfinite(old)): continue
            for q, e in enumerate(vals):
                if used + q <= n and (ndp[used+q] is None or old+e < ndp[used+q]): ndp[used+q] = old+e
        dp = ndp
    cuts = len(parts) - 1
    penalty = (Fraction(2) * abs(int(t)) * cuts) if exact_mode else (2.0 * abs(t) * cuts)
    return dp[n] - penalty, dp[n], penalty, parts


def overlap_lower(nsite, max_cluster, t, u, backend="numpy"):
    """Lower candidate from overlapping windows with term ownership.

    Windows start every ``k-1`` sites.  Each onsite and bond term is assigned
    to exactly one window, so the local Hamiltonians sum to the global one.
    Their independent ground energies therefore form a lower bound candidate;
    no boundary norm penalty is added.  Particle numbers are minimized
    independently because overlapping windows do not share a conserved global
    sector.  This is deliberately a stress test, not a claimed optimal rule.
    """
    k = max_cluster; starts = list(range(0, nsite, max(1, k-1)))
    windows = [(s, min(nsite, s+k)) for s in starts]
    onsite_owner = {}
    bond_owner = {}
    for i in range(nsite):
        onsite_owner[i] = next(j for j,(a,b) in enumerate(windows) if a <= i < b)
    for i in range(nsite-1):
        bond_owner[i] = next(j for j,(a,b) in enumerate(windows) if a <= i and i+1 < b)
    total = 0.0; sizes=[]; cache={}
    for wi,(a,b) in enumerate(windows):
        size=b-a; sizes.append(size)
        # Build local Hamiltonian by removing unowned terms from the ordinary
        # window Hamiltonian; dense correction is small for the intended k<=4.
        best=math.inf
        for q in range(2*size+1):
            key=(size,q,tuple(i for i,o in onsite_owner.items() if o==wi),tuple(i for i,o in bond_owner.items() if o==wi),t,u)
            if key not in cache:
                h=hubbard(size,q,t,u).toarray()
                states=basis(2*size,q)
                for col,x in enumerate(states):
                    owned=sum(((x>>(2*i))&1) and ((x>>(2*i+1))&1) for i in range(size) if onsite_owner.get(a+i)==wi)
                    allons=sum(((x>>(2*i))&1) and ((x>>(2*i+1))&1) for i in range(size))
                    h[col,col] += u*(owned-allons)
                # zero hopping entries for unowned local bonds
                owned_bonds={i-a for i,o in bond_owner.items() if o==wi}
                for i in range(size-1):
                    if i not in owned_bonds:
                        for col,x in enumerate(states):
                            for spin in (0,1):
                                for dst,src in ((2*i+spin,2*(i+1)+spin),(2*(i+1)+spin,2*i+spin)):
                                    out=hop(x,dst,src)
                                    if out: h[basis(2*size,q).index(out[0]),col] += t*out[1]
                cache[key]=float(np.linalg.eigvalsh(h)[0])
            best=min(best,cache[key])
        total += best
    return total, sizes


def run(args):
    out = []
    exact_cache = {}
    for L in args.sites:
        N = L  # half filling
        for ratio in args.couplings:
            t, u = 1.0, ratio
            tic = time.perf_counter()
            exact = ground(L, N, t, u) if L <= args.global_max_sites else None
            full_s = time.perf_counter() - tic
            for k in args.clusters:
                if k > L: continue
                if args.exact and k > 3:
                    continue
                tic = time.perf_counter(); lb, dec, pen, parts = cluster_lower(L,N,k,t,u,args.backend,args.exact,exact_cache)
                lb_exact = str(lb) if args.exact else None
                dec_exact = str(dec) if args.exact else None
                pen_exact = str(pen) if args.exact else None
                out.append({"sites":L,"particles":N,"U_over_t":ratio,"cluster":k,
                            "parts":parts,"exact_energy":exact,"lower_bound":float(lb),
                            "decoupled_energy":float(dec),"omitted_coupling_penalty":float(pen),
                            "lower_bound_exact":lb_exact,"decoupled_energy_exact":dec_exact,"omitted_penalty_exact":pen_exact,
                            "product_variational_upper":0.0 if u >= 0 else None,
                            "omitted_penalty_width":float(pen),
                            "gap":float(exact-lb) if exact is not None else None,
                            "relative_gap":float((exact-lb)/max(1,abs(exact))) if exact is not None else None,
                            "global_truth_status":"numerical_exact_reference" if exact is not None else "not_computed",
                            "certificate_status":"numerical_cluster_eigenvalues; operator_penalty_exact",
                            "exact_replay": bool(args.exact),
                            "product_occupation_upper": 0.0 if u >= 0 else None,
                            "full_seconds":full_s,"certificate_seconds":time.perf_counter()-tic,
                            "states_full":math.comb(2*L,N),
                            "cluster_state_max":max(math.comb(2*p,p) for p in parts)})
                if k <= 4:
                    ov, windows = overlap_lower(L,k,t,u,args.backend)
                    out.append({"method":"overlap_owned_terms","sites":L,"particles":N,
                                "U_over_t":ratio,"cluster":k,"overlap_lower_candidate":ov,
                                "overlap_windows":windows,"disjoint_lower_candidate":float(lb),
                                "overlap_minus_disjoint":ov-float(lb),
                                "certificate_status":"numerical local eigensolves; exact term ownership"})
    return out


if __name__ == "__main__":
    ap=argparse.ArgumentParser(); ap.add_argument("--sites", nargs="+", type=int, default=[4,6,8])
    ap.add_argument("--clusters", nargs="+", type=int, default=[2,3,4]); ap.add_argument("--couplings", nargs="+", type=float, default=[0,2,4])
    ap.add_argument("--global-max-sites", type=int, default=8)
    ap.add_argument("--backend", choices=["numpy","cupy"], default="numpy")
    ap.add_argument("--exact", action="store_true", help="fractional LDL replay for k<=3 sectors")
    ap.add_argument("--output", type=Path, default=Path("results/certificate_scaling/locality_compression/summary.json"))
    a=ap.parse_args(); data=run(a); a.output.parent.mkdir(parents=True,exist_ok=True); a.output.write_text(json.dumps(data,indent=2)+"\n")
    print(json.dumps({"records":len(data),"output":str(a.output)}))
