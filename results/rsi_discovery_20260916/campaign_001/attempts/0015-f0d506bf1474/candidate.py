def propose(payload):
    rows = []
    for block in payload["blocks"]:
        energy, vector = eigenpair(block["H"], payload["den"], True)
        lower = rational(round(energy * 10000000000), 10000000000) - rational("8/10000")
        for digits in [7, 9, 12]:
            factor = cholesky(block["H"], payload["den"], lower + rational("1/1000000"), digits)
            if acquired_margin(block["H"], payload["den"], lower, factor, 10 ** digits) >= 0:
                break
        rows.append({"alpha_count": block["alpha"], "lower_Ha": str(lower),
                     "factor_denominator": 10 ** digits, "factor": factor,
                     "upper_vector": [int(round(x * 10000000000)) for x in vector]})
    return {"blocks": rows}
