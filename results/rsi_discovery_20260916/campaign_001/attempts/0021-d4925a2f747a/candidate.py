def propose(payload):
    rows = []
    for block in payload["blocks"]:
        energy, vector = eigenpair(block["H"], payload["den"], False)
        lower = rational(round(energy * 10000000000), 10000000000) - rational("8/10000")
        factor = cholesky(block["H"], payload["den"], lower + rational("1/10000000"), 3)
        rows.append({"alpha_count": block["alpha"], "lower_Ha": str(lower),
                     "factor_denominator": 10 ** 3, "factor": factor,
                     "upper_vector": [int(round(x * 10000000000)) for x in vector]})
    return {"blocks": rows}
