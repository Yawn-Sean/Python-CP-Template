def getMu(n):
    mu = [0] * (n + 1)
    flg = [0] * (n + 1)
    p = [0] * (n + 1)
    tot = 0
    mu[1] = 1
    for i in range(2, n + 1):
        if flg[i] == 0:
            tot = tot + 1; p[tot] = i; mu[i] = -1
        j = 1
        while j <= tot and i * p[j] <= n:
            flg[i * p[j]] = 1
            if i % p[j] == 0:
                mu[i * p[j]] = 0
                break
            mu[i * p[j]] =  - mu[i]
            j = j + 1
    return mu
