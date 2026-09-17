def propose(payload):
    factor = payload["factor"]
    return [[0] * len(factor) for row in factor]
