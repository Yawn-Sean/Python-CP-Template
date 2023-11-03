def matrix_mul(A, B, mod=10 ** 9 + 7):
    n, m = len(A), len(A[0])
    p = len(B[0])
    ans = [[0] * p for _ in range(n)]
    for i in range(n):
        for j in range(m):
            for k in range(p):
                ans[i][k] += A[i][j] * B[j][k]
                ans[i][k] %= mod
    return ans

def matrix_pow(A, n):
    length = len(A)
    tmp = A
    ans = [[0] * length for _ in range(length)]
    for i in range(length):
        ans[i][i] = 1
    for i in range(60):
        if n % 2:
            ans = matrix_mul(ans, tmp)
        tmp = matrix_mul(tmp, tmp)
        n //= 2
    return ans