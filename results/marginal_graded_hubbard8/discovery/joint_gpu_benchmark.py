"""Bounded CPU/CuPy benchmark for the actual joint projector search.

This is a timing and numerical-equivalence harness only.  Rational replay of
the Gram bound and exact energy acceptance remain separate gates.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse, hashlib, json, math, os, time
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np
from scipy.linalg import eigvalsh

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / "results/marginal_graded_hubbard8"

from experiments.marginal_local_hubbard_block import sector_matrices
from experiments.marginal_projector_extendibility import _projector_vector
from experiments.marginal_charged_projectors import charged_vectors, joint_overlap_grams


def _profiles(x):
    a, b, p, q, d, e = map(F, x)
    return ([a, b, 10-a-b, 10-a-b, b, a],
            [p, q, 5-2*p-2*q, q, p],
            [d, e, F(5, 2)-2*d-2*e, e, d])


def _sectors(x, local):
    u, t, v = _profiles(x)
    return sector_matrices(6, F(local["U"]), F(local["t"]), u, t,
                           F(local["V"]), v)


def _resident(x, local, half, half_norm, charged_family, charged_norm):
    mats = []
    for key, (_, K, cols) in _sectors(x, local).items():
        scale = np.sqrt(np.asarray([sum(a*a for a in c.values()) for c in cols], dtype=np.float64))
        A = np.asarray(K, dtype=np.float64) / scale[:, None] / scale[None, :]
        def projector(vector, norm):
            w = np.asarray([sum(a*vector.get(s, 0) for s, a in c.items()) for c in cols], dtype=np.float64)
            w /= scale * math.sqrt(float(norm))
            return np.outer(w, w)
        P = projector(half, half_norm)
        Q = sum((projector(v, charged_norm) for v in charged_family), np.zeros_like(A))
        mats.append((key, A, P, Q))
    return mats


def _profiles_resident(origin, local, half, half_norm, charged_family, charged_norm, step):
    base = _resident(origin, local, half, half_norm, charged_family, charged_norm)
    derivatives = []
    for j in range(6):
        point = list(origin); point[j] += F(step)
        pert = _resident(point, local, half, half_norm, charged_family, charged_norm)
        derivatives.append([(b[1]-a[1]) / float(step) for a, b in zip(base, pert)])
    return base, derivatives


def _gram_batch(source_half, source_charged, ratio):
    groups = joint_overlap_grams(source_half, source_charged, 4)
    blocks = []
    for data in groups.values():
        gram = np.asarray(data["gram"], dtype=np.float64)
        norms = np.asarray(data["norms"], dtype=np.float64)
        source = np.asarray(data["sources"])
        weights = norms / np.where(source == 0, 1.0, ratio)
        blocks.append(gram / np.sqrt(weights[:, None] * weights[None, :]))
    return blocks


def _cpu_evals(mats, derivatives, points, alpha, beta, ratio, batch):
    out = []
    origin = np.asarray(points[0], dtype=np.float64)
    for start in range(0, len(points), batch):
        delta = np.asarray(points[start:start+batch], dtype=np.float64) - origin
        mins = np.full(len(delta), np.inf)
        for si, (_, A, P, Q) in enumerate(mats):
            D = np.einsum('bi,ijk->bjk', delta, np.asarray([derivatives[j][si] for j in range(6)]))
            vals = np.linalg.eigvalsh(A + D + (alpha + beta)*P + beta*ratio*Q)
            mins = np.minimum(mins, vals[:, 0])
        out.extend(mins.tolist())
    return np.asarray(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--evaluations", type=int, default=12)
    ap.add_argument("--batch", type=int, default=8)
    ap.add_argument("--ratio", type=float, choices=(0.5,1.0,2.0), default=0.5)
    ap.add_argument("--out", type=Path, default=BASE / "gpu_benchmark")
    args = ap.parse_args()
    if args.evaluations < 1 or args.evaluations > 64 or args.batch < 1 or args.batch > 16 or args.ratio <= 0:
        raise ValueError("positive evaluations, batch and ratio required")
    cert = json.loads((BASE / "six_site_projector/refined_certificate.json").read_text())
    local = cert["local_window"]
    origin = [F(local["onsite_profile"][0]), F(local["onsite_profile"][1]),
              F(local["hopping_profile"][0]), F(local["hopping_profile"][1]),
              F(local["density_profile"][0]), F(local["density_profile"][1])]
    half, half_norm = _projector_vector(cert["vector"], 6)
    charged_source = json.loads((BASE / "charged_projector/source.json").read_text())["physical_state"]
    charged_family, charged_norm = charged_vectors(charged_source)
    prep0 = time.perf_counter()
    mats, derivatives = _profiles_resident(origin, local, half, half_norm, charged_family, charged_norm, F(1, 1000))
    gram_blocks = _gram_batch(cert["vector"], charged_source, args.ratio)
    prep = time.perf_counter() - prep0
    rng = np.random.default_rng(20260912)
    points = [np.asarray(origin, dtype=np.float64)]
    points.extend(np.asarray(origin, dtype=np.float64) + rng.normal(0, 2e-3, 6) for _ in range(args.evaluations-1))
    alpha, beta = float(F(cert["penalty"])), 0.001
    theta_h = float(F(cert["projector_sum_ceiling"]) / cert["windows"])
    theta_joint = {0.5:0.595736407,1.0:0.759879428,2.0:1.23137232075}[args.ratio]
    t0 = time.perf_counter(); cpu = _cpu_evals(mats, derivatives, points, alpha, beta, args.ratio, args.batch); cpu_time = time.perf_counter()-t0
    gram_t0 = time.perf_counter(); cpu_gram = np.asarray([eigvalsh(b, subset_by_index=[len(b)-1, len(b)-1])[0] for b in gram_blocks]); gram_cpu_time = time.perf_counter()-gram_t0
    reference_t0=time.perf_counter();reference=[]
    for point in points:
        delta=point-points[0]
        reference.append(min(float(eigvalsh(A+sum((delta[j]*derivatives[j][si] for j in range(6)),np.zeros_like(A))+(alpha+beta)*P+beta*args.ratio*Q,subset_by_index=[0,0])[0]) for si,(_,A,P,Q) in enumerate(mats)))
    reference_time=time.perf_counter()-reference_t0
    if not np.allclose(cpu,reference,atol=1e-10,rtol=1e-10):raise ValueError('Batched CPU assembly disagrees with scalar reference')
    gpu = {"available": False, "reason": "CuPy not attempted unless installed with a visible CUDA device"}
    try:
        import cupy as cp
        if cp.cuda.runtime.getDeviceCount()<1:raise RuntimeError('CuPy installed but no visible CUDA device')
        device = cp.cuda.Device(); device.use()
        sync = cp.cuda.Stream.null.synchronize
        transfer0 = time.perf_counter()
        dmats = []
        for _, A, P, Q in mats:
            dder = np.asarray([derivatives[j][len(dmats)] for j in range(6)])
            dmats.append((cp.asarray(A), cp.asarray(P), cp.asarray(Q), cp.asarray(dder)))
        dblocks = [cp.asarray(b) for b in gram_blocks];dpoints=cp.asarray(points);sync();transfer=time.perf_counter()-transfer0
        warm0=time.perf_counter();largest=max(dmats,key=lambda m:len(m[0]))[0]
        cp.linalg.eigvalsh(cp.stack([largest,largest]));sync();warmup=time.perf_counter()-warm0
        kt = time.perf_counter(); gpu_rows=[]
        dorigin = dpoints[0]
        for start in range(0, len(points), args.batch):
            delta = dpoints[start:start+args.batch] - dorigin
            gmin = cp.full((len(delta),), cp.inf, dtype=cp.float64)
            for A, P, Q, dder in dmats:
                D = cp.einsum('bi,ijk->bjk', delta, dder)
                vals = cp.linalg.eigvalsh(A + D + (alpha+beta)*P + beta*args.ratio*Q)
                gmin = cp.minimum(gmin, vals[:, 0])
            gpu_rows.append(gmin)
        device_profile=cp.concatenate(gpu_rows);sync();gpu_profile_kernel=time.perf_counter()-kt
        kt=time.perf_counter();gvals=cp.stack([cp.linalg.eigvalsh(b)[-1] for b in dblocks]);sync();gpu_gram_time=time.perf_counter()-kt
        copy0=time.perf_counter();gpu_profile=cp.asnumpy(device_profile);gpu_gram=cp.asnumpy(gvals);sync();copy_time=time.perf_counter()-copy0
        if not np.allclose(gpu_profile,reference,atol=1e-9,rtol=1e-10) or not np.allclose(gpu_gram,cpu_gram,atol=1e-9,rtol=1e-10):
            raise ValueError('GPU eigenvalues fail CPU equivalence tolerance')
        gpu = {"available": True, "device": str(device), "transfer_seconds": transfer, "warmup_seconds": warmup, "profile_kernel_seconds": gpu_profile_kernel, "gram_kernel_seconds": gpu_gram_time,"result_copy_seconds":copy_time,"total_end_to_end_seconds":prep+transfer+warmup+gpu_profile_kernel+gpu_gram_time+copy_time,"profile_max_abs_difference":float(np.max(np.abs(gpu_profile-reference))),"gram_max_abs_difference":float(np.max(np.abs(gpu_gram-cpu_gram))),"validated":True,"resident_bytes":int(sum(A.nbytes+P.nbytes+Q.nbytes+d.nbytes for A,P,Q,d in dmats)+sum(b.nbytes for b in dblocks)+dpoints.nbytes),"matrix_batch_bytes":args.batch*max(len(A) for A,_,_,_ in dmats)**2*8,"solver_workspace":"Not measured; resident_bytes excludes vendor eigensolver workspace."}
    except (ImportError, ModuleNotFoundError) as exc:
        gpu = {"available": False, "reason": f"{type(exc).__name__}: {exc}"}
    args.out.mkdir(parents=True, exist_ok=True)
    receipt = {"source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), "evaluations": args.evaluations, "batch": args.batch,"ratio":args.ratio, "sectors": len(mats), "gram_blocks": len(gram_blocks), "gram_columns": int(sum(len(b) for b in gram_blocks)), "prep_seconds": prep, "cpu": {"profile_kernel_seconds": cpu_time,"scalar_subset_reference_seconds":reference_time,"batched_max_abs_difference":float(np.max(np.abs(cpu-reference))),"total_end_to_end_seconds":prep+reference_time+gram_cpu_time, "gram_kernel_seconds": gram_cpu_time, "profile_values": cpu.tolist(), "gram_maxima": cpu_gram.tolist(), "density_objective_values": ((cpu-alpha*theta_h-beta*theta_joint)/5).tolist(), "environment": {"OPENBLAS_NUM_THREADS": os.environ.get("OPENBLAS_NUM_THREADS")}}, "gpu": gpu, "objective_direction": "larger lower-bound density is better; numerical only", "normalization": "Reflection blocks use inverse square-root orbit norms; weighted Gram is symmetrically normalized on both axes.", "penalties": {"alpha": alpha, "beta": beta, "theta_half": theta_h, "theta_joint": theta_joint}, "scope": "bounded float64 resident-matrix benchmark; no exact certificate or speedup claim"}
    (args.out / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()
