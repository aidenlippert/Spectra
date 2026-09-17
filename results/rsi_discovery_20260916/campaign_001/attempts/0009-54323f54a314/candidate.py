def propose(payload):
    factor = payload["factor"]
    n = len(factor)
    result = [[0] * n for i in range(n)]
    for i in range(n):
        for j in range(i + 1):
            result[i][j] = sum(x * y for x, y in zip(factor[i], factor[j]))
            result[j][i] = result[i][j]
    return result
