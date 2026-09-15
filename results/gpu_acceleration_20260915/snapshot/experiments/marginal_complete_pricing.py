"""Complete, bounded numerical pricing over principal supports of size 2--4.

This is a diagnostic for a residual matrix.  It enumerates every support and
uses a symmetric eigensolver on each small principal block; it deliberately
does not claim an exact PSD certificate (the input and arithmetic are float).
"""
from itertools import combinations, islice
import heapq
import math


def scan_supports(matrix, *, max_support=4, batch_size=65536, top_k=256,
                  tolerance=1e-8):
    """Exhaustively scan principal blocks of orders 2..``max_support``.

    Returns a JSON-friendly dictionary. ``negative_supports`` contains at most
    ``top_k`` most negative eigenpairs, while ``checked_supports`` is exact.
    The scan has bounded memory and never materializes all quadruples.
    """
    import numpy as np
    a = np.asarray(matrix, dtype=float)
    if a.ndim != 2 or a.shape[0] != a.shape[1] or not np.all(np.isfinite(a)):
        raise ValueError("matrix must be finite and square")
    if not np.allclose(a, a.T, rtol=1e-10, atol=1e-12):
        raise ValueError("matrix must be symmetric")
    n = int(a.shape[0])
    if type(max_support) is not int or not 2 <= max_support <= 4:
        raise ValueError("max_support must be between 2 and min(4,n)")
    max_support = min(int(max_support), n)
    if max_support < 2:
        raise ValueError("matrix must have at least two rows")
    if type(batch_size) is not int or batch_size < 1 or type(top_k) is not int or top_k < 0:
        raise ValueError("batch_size and top_k must be nonnegative bounded integers")
    if isinstance(tolerance, bool) or not math.isfinite(float(tolerance)) or tolerance < 0:
        raise ValueError("tolerance must be finite and nonnegative")

    counts = {}
    minimum = math.inf
    near_boundary = 0
    heap = []
    negative_total = 0
    actual_counts = {}
    per_order = {}
    per_order_heaps = {}
    for order in range(2, int(max_support) + 1):
        total = math.comb(n, order)
        counts[str(order)] = total
        actual_counts[str(order)] = 0
        per_order[str(order)] = {"checked": 0, "negative": 0,
                                 "near_boundary": 0, "minimum_eigenvalue": math.inf}
        per_order_heaps[str(order)] = []
        iterator = combinations(range(n), order)
        while True:
            supports = list(islice(iterator, batch_size))
            if not supports:
                break
            actual_counts[str(order)] += len(supports)
            per_order[str(order)]["checked"] += len(supports)
            idx = np.asarray(supports, dtype=int)
            blocks = a[idx[:, :, None], idx[:, None, :]]
            vals, vecs = np.linalg.eigh(blocks)
            lows = vals[:, 0]
            minimum = min(minimum, float(np.min(lows)))
            local_min = float(np.min(lows))
            near = int(np.count_nonzero(np.abs(lows) <= tolerance))
            neg = int(np.count_nonzero(lows < -tolerance))
            per_order[str(order)]["minimum_eigenvalue"] = min(per_order[str(order)]["minimum_eigenvalue"], local_min)
            per_order[str(order)]["near_boundary"] += near
            per_order[str(order)]["negative"] += neg
            near_boundary += near
            negative_total += neg
            if top_k:
                candidates=np.flatnonzero(lows < -tolerance)
                if len(candidates)>top_k:
                    cutoff=np.partition(lows[candidates],top_k-1)[top_k-1]
                    better=candidates[lows[candidates]<cutoff]
                    # Match heap's lexicographic support tie break. Batch
                    # supports arrive in lexicographic order.
                    ties=candidates[lows[candidates]==cutoff]
                    candidates=np.concatenate((better,ties[-(top_k-len(better)):]))
                # A discarded candidate has at least top_k better candidates
                # in its own batch, so cannot enter either final top_k list.
                for row in candidates:
                    item = (float(lows[row]), tuple(int(x) for x in supports[row]),
                            vecs[row, :, 0].astype(float).tolist())
                    key = (-item[0], item[1])
                    if len(heap) < top_k:
                        heapq.heappush(heap, (key, item))
                    # The heap root is the least negative retained value
                    # (smallest ``-eigenvalue``); evict it for a more
                    # negative direction.
                    elif key > heap[0][0]:
                        heapq.heapreplace(heap, (key, item))
                    oh = per_order_heaps[str(order)]
                    if len(oh) < top_k:
                        heapq.heappush(oh, (key, item))
                    elif key > oh[0][0]:
                        heapq.heapreplace(oh, (key, item))
    items = [x[1] for x in heap]
    items.sort(key=lambda x: (x[0], x[1]))
    return {
        "dimension": n, "max_support": int(max_support),
        "checked_supports": actual_counts, "possible_supports": counts,
        "checked_total": sum(actual_counts.values()), "negative_support_count": negative_total,
        "per_order_stats": per_order,
        "minimum_eigenvalue": float(minimum),
        "near_boundary_count": near_boundary, "tolerance": float(tolerance),
        "negative_supports": [
            {"eigenvalue": x[0], "support": list(x[1]), "direction": x[2]}
            for x in items],
        "negative_supports_by_order": {
            order: [{"eigenvalue": x[0], "support": list(x[1]), "direction": x[2]}
                    for x in sorted((entry[1] for entry in heap), key=lambda x: (x[0], x[1]))]
            for order, heap in per_order_heaps.items()},
        "scope": "Complete finite support coverage, numerical only; no exact PSD or representability claim.",
    }


complete_price = scan_supports
