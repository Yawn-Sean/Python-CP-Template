# sum(i = 0 to n - 1) [(a * i + b) / m]
def floor_sum_unsigned(n, m, a, b):
    # mod = ...
    ans = 0
    while True:
        if a >= m:
            ans += n * (n - 1) * (a // m) // 2
            a %= m
        if b >= m:
            ans += n * (b // m)
            b %= m
        # if ans >= mod: ans %= mod
        y_max = a * n + b
        if y_max < m: break
        n, b, m, a, = y_max // m, y_max % m, a, m
    return ans