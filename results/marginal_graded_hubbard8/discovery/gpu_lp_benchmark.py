"""Benchmark a frozen Reynolds restricted LP on CPU and optional cuOpt.

The frozen matrix is rebuilt from the accepted crossover proposal, so both
backends (when available) see the same row/column ordering.  This is a
benchmark harness, not a certificate producer: exact polynomial replay stays
the acceptance gate.
"""
from pathlib import Path
import argparse, json, time

import numpy as np
from scipy.sparse import csc_matrix, hstack, eye

from experiments.marginal_reynolds_pricing import ReynoldsMomentDictionary
from experiments.marginal_fixed_metric_pricing import IncrementalMaster
from experiments.marginal_joint_coefficient_constructor import input_digest


ROOT = Path(__file__).resolve().parents[1]


def build_matrix(data, proposal, checkpoint):
    gamma = proposal["gamma"]
    dictionary = ReynoldsMomentDictionary(data, gamma, feature_degree=2, proof_degree=6)
    active = [(int(block), tuple(label)) for block, label in checkpoint["active"]]
    qrows = dictionary.qrows
    rows = 2 * qrows + 1
    metric = dictionary.metric
    columns = []
    for block, label in active:
            c = np.zeros(rows)
            c[block * qrows:(block + 1) * qrows] = -dictionary.column(label)
            columns.append(c)
    atom = csc_matrix(np.asarray(columns, dtype=float).T)
    # The actual joint master has free metric coefficients first, then atom
    # columns, then +/- slacks.  Keep this logical order when computing
    # residuals; HiGHS internally stores the initial slacks first.
    full = hstack((dictionary.metric, -dictionary.metric, atom,
                   eye(rows, format="csc"), -eye(rows, format="csc")), format="csc")
    rhs = np.r_[0.0009 * dictionary.one, 0.0001 * dictionary.one, 1.0]
    return dictionary, full, rhs, len(columns), input_digest(data), active


def cpu_run(matrix, rhs, seconds, mps_path):
    started = time.monotonic()
    master = IncrementalMaster(rhs)
    rows = matrix.shape[0]
    # The final normalization equality must not be satisfied by a slack.
    for col in (rows - 1, 2 * rows - 1):
        master.solver.changeColBounds(col, 0.0, 0.0)
    for key, value in [('solver', 'ipm'), ('run_crossover', 'on'),
                       ('primal_feasibility_tolerance', 1e-9),
                       ('dual_feasibility_tolerance', 1e-9),
                       ('ipm_optimality_tolerance', 1e-9)]:
        master.solver.setOptionValue(key, value)
    # IncrementalMaster stores base slacks first, while its returned vector
    # places appended columns first.  `fresh` is metric+, metric-, atoms and
    # therefore aligns with the logical matrix used below.
    metric_atom_count = matrix.shape[1] - 2 * rows
    appended = matrix[:, :metric_atom_count]
    result = master.run(appended, seconds)
    master.solver.writeModel(str(mps_path))
    x = np.asarray(result.x)
    residual = np.asarray(matrix @ x - rhs)
    return {
        "backend": "cpu_highs",
        "status": result.message,
        "success": bool(result.success),
        "objective": float(result.fun),
        "max_abs_residual": float(np.max(np.abs(residual))) if residual.size else 0.0,
        "elapsed_seconds": time.monotonic() - started,
        "mps": str(mps_path),
    }


def optional_cuopt(mps_path, seconds, matrix, rhs, method='barrier'):
    """Run cuOpt only when installed; return a structured unavailable result."""
    try:
        from cuopt.linear_programming.problem import Problem
        from cuopt.linear_programming.solver_settings import SolverSettings
    except Exception as error:
        return {"backend": "cuopt", "available": False,
                "reason": f"cuOpt unavailable: {type(error).__name__}: {error}"}
    started = time.monotonic()
    try:
        # NVIDIA's public modeling API: readMPS is a classmethod; solve takes
        # SolverSettings. This path still requires a real GPU validation run.
        problem = Problem.readMPS(str(mps_path))
        settings=SolverSettings()
        settings.set_parameter('time_limit',float(seconds))
        settings.set_parameter('method',method)
        settings.set_parameter('crossover',True)
        settings.set_optimality_tolerance(1e-9)
        problem.solve(settings)
        variables={v.VariableName:float(v.Value) for v in problem.getVariables()}
        constraints={c.ConstraintName:float(c.DualValue) for c in problem.getConstraints()}
        # HiGHS MPS writer names native columns c0..cN and rows r0..rM.
        # Recover by name, so MPS parser reorderings cannot corrupt pricing.
        n,m=matrix.shape[1],matrix.shape[0]
        if set(variables)!={f'c{i}' for i in range(n)} or set(constraints)!={f'r{i}' for i in range(m)}:
            raise ValueError('MPS row/column names differ from frozen native model')
        native=np.array([variables[f'c{i}'] for i in range(n)])
        x=np.r_[native[2*m:],native[:2*m]]
        dual=np.array([constraints[f'r{i}'] for i in range(m)])
        if not np.all(np.isfinite(x)) or not np.all(np.isfinite(dual)):
            raise ValueError('Backend returned nonfinite primal or dual values')
        solution=mps_path.parent/f'cuopt_{method}_solution.json'
        solution.write_text(json.dumps({'logical_primal':x.tolist(),'row_dual':dual.tolist()},indent=2)+'\n')
        return {"backend": "cuopt", "available": True,
                'method':method,'objective':float(problem.ObjValue),
                'max_abs_residual':float(np.max(np.abs(matrix@x-rhs))),
                'max_bound_violation':float(max(0.,-min(x),abs(native[m-1]),abs(native[2*m-1]))),
                'solution':str(solution),'exact_certificate_accepted':False,
                "elapsed_seconds": time.monotonic() - started,
                "status": str(problem.Status)}
    except Exception as error:
        return {"backend": "cuopt", "available": True,
                "error": f"{type(error).__name__}: {error}",
                "elapsed_seconds": time.monotonic() - started}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--proposal", type=Path,
                        default=ROOT / "joint_reynolds_crossover/proposal.json")
    parser.add_argument("--data", type=Path, default=ROOT / "hamiltonian.json")
    parser.add_argument("--out", type=Path, default=ROOT / "gpu_benchmark")
    parser.add_argument("--seconds", type=float, default=120)
    parser.add_argument("--cuopt", action="store_true")
    parser.add_argument('--cuopt-method',choices=('barrier','pdlp'),default='barrier')
    args = parser.parse_args()
    if args.seconds <= 0 or not args.proposal.exists() or not args.data.exists():
        raise ValueError("Existing proposal/data and positive time required")
    proposal = json.loads(args.proposal.read_text())
    checkpoint = json.loads((ROOT / "joint_reynolds_crossover/checkpoint.json").read_text())
    data = json.loads(args.data.read_text())
    args.out.mkdir(parents=True, exist_ok=True)
    dictionary, matrix, rhs, atom_count, digest, active = build_matrix(data, proposal, checkpoint)
    mps = args.out / "frozen_restricted.mps"
    results = [cpu_run(matrix, rhs, args.seconds, mps)]
    if args.cuopt:
        results.append(optional_cuopt(mps, args.seconds,matrix,rhs,args.cuopt_method))
    receipt = {"source_sha256": digest, "proposal": str(args.proposal),
               "rows": int(matrix.shape[0]), "columns": int(matrix.shape[1]),
               "atom_columns": atom_count, "nonzeros": int(matrix.nnz),
               "ordering": "checkpoint active metric+, metric-, atom columns; returned CPU vector reorders appended columns before initial slacks; MPS uses native initial-slack order",
               "results": results,
               "scope": "Frozen numerical restricted LP comparison. CPU result is diagnostic; exact export/replay remains the certificate gate. GPU is reported only when cuOpt is already installed."}
    (args.out / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()
