n, mod = ...

def ext_gcd(a, b):
    if b == 0: return a, 1, 0
    gcd, x, y = ext_gcd(b, a % b)
    return gcd, y, x - (a // b) * y

def mod_inv(a, p):
    gcd, x, y = ext_gcd(a, p)
    return x % p

def count_power(x, p):
    ans = 0
    while x:
        x //= p
        ans += x
    return ans

def c_mod(x, y):
    if y > x or y < 0: return 0
    
    curr, curr_mod = 0, 1
    for p, c, v in prime_factors:
        cnt = min(c, count_power(x, p) - count_power(y, p) - count_power(x - y, p))
        others = facts[v].comb(x, y)
        
        new = pow(p, cnt, v) * others % v
        new_mod = v
        
        curr = (curr - new) * mod_inv(new_mod, curr_mod) * new_mod + new
        curr_mod *= v
        curr %= curr_mod
    return curr

v = mod
prime_factors = []
vals = []
for i in count(2):
    if i * i > v: break
    if v % i == 0:
        cnt = 0
        x = 1
        while v % i == 0:
            v //= i
            cnt += 1
            x *= i
        prime_factors.append((i, cnt, x))
if v > 1: prime_factors.append((v, 1, v))

facts = {p: Factorial(n, p, h) for h, _, p in prime_factors}