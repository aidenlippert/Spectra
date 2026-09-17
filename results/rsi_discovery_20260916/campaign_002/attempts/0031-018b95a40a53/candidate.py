def propose(payload):
    rows = []
    for block in payload["blocks"]:
        energy, vector = eigenpair(block["H"], payload["den"], True)
        lower = rational(round(energy * 10000000000), 10000000000) - rational("8/10000")
        for digits, fd in [(7, 10000000), (9, 1000000000), (12, 1000000000000)]:
            factor = cholesky(block["H"], payload["den"], lower + rational("1/1000000"), digits)
            if acquired_margin(block["H"], payload["den"], lower, factor, fd) >= 0:
                break
        rows.append({"alpha_count": block["alpha"], "lower_Ha": str(lower),
                     "factor_denominator": fd, "factor": factor,
                     "upper_vector": [int(round(x * 10000000000)) for x in vector]})
    return {"blocks": rows}
