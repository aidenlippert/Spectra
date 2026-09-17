def propose(payload):
    rows = []
    for block in payload["blocks"]:
        energy, vector = eigenpair(block["H"], payload["den"], False)
        lower = rational(round(energy * 10000000000), 10000000000) - rational("12/10000")
        factor = cholesky(block["H"], payload["den"], lower + rational("1/100000000"), 11)
        upper = [int(round(x * 100000000)) for x in vector]
        rows.append({"alpha_count": block["alpha"], "lower_Ha": str(lower),
                     "factor_denominator": 10 ** 11, "factor": factor,
                     "upper_vector": upper})
    return {"blocks": rows}