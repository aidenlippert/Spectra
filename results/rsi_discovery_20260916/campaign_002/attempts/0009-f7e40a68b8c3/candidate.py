def propose(payload):
    factor = payload["factor"]
    return [[sum(x * y for x, y in zip(a, b)) for b in factor] for a in factor]
