poly = 1

if poly:

    # square root of a single value
    def Tonelli_Shanks(N, p):
        if pow(N, p >> 1, p) == p - 1:
            return None
        if p%4==3:
            return pow(N, (p + 1) // 4, p)
        for nonresidue in range(1, p):
            if pow(nonresidue, p>>1, p)==p-1:
                break
        pp = p - 1
        cnt = 0
        while pp % 2 == 0:
            pp //= 2
            cnt += 1
        s = pow(N, pp, p)
        retu = pow(N, (pp + 1) // 2,p)
        for i in range(cnt - 2, -1, -1):
            if pow(s, 1 <<i, p) == p - 1:
                s *= pow(nonresidue, p >> 1 + i, p)
                s %= p
                retu *= pow(nonresidue, p >> 2 + i, p)
                retu %= p
        return retu

    mod = 998244353
    imag = 911660635
    iimag = 86583718
    rate2 = (911660635, 509520358, 369330050, 332049552, 983190778, 123842337, 238493703, 975955924, 603855026, 856644456, 131300601,
                842657263, 730768835, 942482514, 806263778, 151565301, 510815449, 503497456, 743006876, 741047443, 56250497, 867605899)
    irate2 = (86583718, 372528824, 373294451, 645684063, 112220581, 692852209, 155456985, 797128860, 90816748, 860285882, 927414960,
                354738543, 109331171, 293255632, 535113200, 308540755, 121186627, 608385704, 438932459, 359477183, 824071951, 103369235)
    rate3 = (372528824, 337190230, 454590761, 816400692, 578227951, 180142363, 83780245, 6597683, 70046822, 623238099,
                183021267, 402682409, 631680428, 344509872, 689220186, 365017329, 774342554, 729444058, 102986190, 128751033, 395565204)
    irate3 = (509520358, 929031873, 170256584, 839780419, 282974284, 395914482, 444904435, 72135471, 638914820, 66769500,
                771127074, 985925487, 262319669, 262341272, 625870173, 768022760, 859816005, 914661783, 430819711, 272774365, 530924681)

    def butterfly(a):
        n = len(a)
        h = (n - 1).bit_length()
        len_ = 0
        while len_ < h:
            if h - len_ == 1:
                p = 1 << (h - len_ - 1)
                rot = 1
                for s in range(1 << len_):
                    offset = s << (h - len_)
                    for i in range(p):
                        l = a[i + offset]
                        r = a[i + offset + p] * rot % mod
                        a[i + offset] = (l + r) % mod
                        a[i + offset + p] = (l - r) % mod
                    if s + 1 != 1 << len_:
                        rot *= rate2[(~s & -~s).bit_length() - 1]
                        rot %= mod
                len_ += 1
            else:
                p = 1 << (h - len_ - 2)
                rot = 1
                for s in range(1 << len_):
                    rot2 = rot * rot % mod
                    rot3 = rot2 * rot % mod
                    offset = s << (h - len_)
                    for i in range(p):
                        a0 = a[i + offset]
                        a1 = a[i + offset + p] * rot
                        a2 = a[i + offset + p * 2] * rot2
                        a3 = a[i + offset + p * 3] * rot3
                        a1na3imag = (a1 - a3) % mod * imag
                        a[i + offset] = (a0 + a2 + a1 + a3) % mod
                        a[i + offset + p] = (a0 + a2 - a1 - a3) % mod
                        a[i + offset + p * 2] = (a0 - a2 + a1na3imag) % mod
                        a[i + offset + p * 3] = (a0 - a2 - a1na3imag) % mod
                    if s + 1 != 1 << len_:
                        rot *= rate3[(~s & -~s).bit_length() - 1]
                        rot %= mod
                len_ += 2

    def butterfly_inv(a):
        n = len(a)
        h = (n - 1).bit_length()
        len_ = h
        while len_:
            if len_ == 1:
                p = 1 << (h - len_)
                irot = 1
                for s in range(1 << (len_ - 1)):
                    offset = s << (h - len_ + 1)
                    for i in range(p):
                        l = a[i + offset]
                        r = a[i + offset + p]
                        a[i + offset] = (l + r) % mod
                        a[i + offset + p] = (l - r) * irot % mod
                    if s + 1 != (1 << (len_ - 1)):
                        irot *= irate2[(~s & -~s).bit_length() - 1]
                        irot %= mod
                len_ -= 1
            else:
                p = 1 << (h - len_)
                irot = 1
                for s in range(1 << (len_ - 2)):
                    irot2 = irot * irot % mod
                    irot3 = irot2 * irot % mod
                    offset = s << (h - len_ + 2)
                    for i in range(p):
                        a0 = a[i + offset]
                        a1 = a[i + offset + p]
                        a2 = a[i + offset + p * 2]
                        a3 = a[i + offset + p * 3]
                        a2na3iimag = (a2 - a3) * iimag % mod
                        a[i + offset] = (a0 + a1 + a2 + a3) % mod
                        a[i + offset + p] = (a0 - a1 + a2na3iimag) * irot % mod
                        a[i + offset + p * 2] = (a0 + a1 - a2 - a3) * irot2 % mod
                        a[i + offset + p * 3] = (a0 - a1 - a2na3iimag) * irot3 % mod
                    if s + 1 != (1 << (len_ - 2)):
                        irot *= irate3[(~s & -~s).bit_length() - 1]
                        irot %= mod
                len_ -= 2

    # integration
    def integrate(a):
        a = a.copy()
        n = len(a)
        assert n > 0
        a.pop()
        a.insert(0, 0)
        inv = [1, 1]
        for i in range(2, n):
            inv.append(-inv[mod % i] * (mod // i) % mod)
            a[i] = a[i] * inv[i] % mod
        return a

    # differentiation
    def differentiate(a):
        n = len(a)
        assert n > 0
        for i in range(2, n):
            a[i] = a[i] * i % mod
        a.pop(0)
        a.append(0)
        return a

    def convolution_naive(a, b):
        n = len(a)
        m = len(b)
        ans = [0] * (n + m - 1)
        if n < m:
            for j in range(m):
                for i in range(n):
                    ans[i + j] = (ans[i + j] + a[i] * b[j]) % mod
        else:
            for i in range(n):
                for j in range(m):
                    ans[i + j] = (ans[i + j] + a[i] * b[j]) % mod
        return ans

    def convolution_ntt(a, b):
        a = a.copy()
        b = b.copy()
        n = len(a)
        m = len(b)
        z = 1 << (n + m - 2).bit_length()
        a += [0] * (z - n)
        butterfly(a)
        b += [0] * (z - m)
        butterfly(b)
        for i in range(z):
            a[i] = a[i] * b[i] % mod
        butterfly_inv(a)
        a = a[:n + m - 1]
        iz = pow(z, mod - 2, mod)
        for i in range(n + m - 1):
            a[i] = a[i] * iz % mod
        return a

    # convolution with itself
    def convolution_square(a):
        a = a.copy()
        n = len(a)
        z = 1 << (2 * n - 2).bit_length()
        a += [0] * (z - n)
        butterfly(a)
        for i in range(z):
            a[i] = a[i] * a[i] % mod
        butterfly_inv(a)
        a = a[:2 * n - 1]
        iz = pow(z, mod - 2, mod)
        for i in range(2 * n - 1):
            a[i] = a[i] * iz % mod
        return a

    # convolution
    def convolution(a, b):
        """It calculates (+, x) convolution in mod 998244353. 
        Given two arrays a[0], a[1], ..., a[n - 1] and b[0], b[1], ..., b[m - 1], 
        it calculates the array c of length n + m - 1, defined by

        >   c[i] = sum(a[j] * b[i - j] for j in range(i + 1)) % 998244353.

        It returns an empty list if at least one of a and b are empty.

        Complexity
        ----------

        >   O(n log n), where n = len(a) + len(b).
        """
        n = len(a)
        m = len(b)
        if n == 0 or m == 0:
            return []
        if min(n, m) <= 60:
            return convolution_naive(a, b)
        if a is b:
            return convolution_square(a)
        return convolution_ntt(a, b)

    # inversion
    def inverse(a):
        n = len(a)
        assert n > 0 and a[0] != 0
        res = [pow(a[0], mod - 2, mod)]
        m = 1
        while m < n:
            f = a[:min(n, 2 * m)] + [0] * (2 * m - min(n, 2*m))
            g = res + [0] * m
            butterfly(f)
            butterfly(g)
            for i in range(2 * m):
                f[i] = f[i] * g[i] % mod
            butterfly_inv(f)
            f = f[m:] + [0] * m
            butterfly(f)
            for i in range(2 * m):
                f[i] = f[i] * g[i] % mod
            butterfly_inv(f)
            iz = pow(2 * m, mod - 2, mod)
            iz = (-iz * iz) % mod
            for i in range(m):
                f[i] = f[i] * iz % mod
            res += f[:m]
            m <<= 1
        return res[:n]

    # logarithm
    def log(a):
        a = a.copy()
        n = len(a)
        assert n > 0 and a[0] == 1
        a_inv = inverse(a)
        a = differentiate(a)
        a = convolution(a, a_inv)[:n]
        a = integrate(a)
        return a

    # exponentiation
    def exp(a):
        a = a.copy()
        n = len(a)
        assert n > 0 and a[0] == 0
        g = [1]
        a[0] = 1
        h_drv = a.copy()
        h_drv = differentiate(h_drv)
        m = 1
        while m < n:
            f_fft = a[:m] + [0] * m
            butterfly(f_fft)

            if m > 1:
                _f = [f_fft[i] * g_fft[i] % mod for i in range(m)]
                butterfly_inv(_f)
                _f = _f[m // 2:] + [0] * (m // 2)
                butterfly(_f)
                for i in range(m):
                    _f[i] = _f[i] * g_fft[i] % mod
                butterfly_inv(_f)
                _f = _f[:m//2]
                iz = pow(m, mod - 2, mod)
                iz *= -iz
                iz %= mod
                for i in range(m//2):
                    _f[i] = _f[i] * iz % mod
                g.extend(_f)

            t = a[:m]
            t = differentiate(t)
            r = h_drv[:m - 1]
            r.append(0)
            butterfly(r)
            for i in range(m):
                r[i] = r[i] * f_fft[i] % mod
            butterfly_inv(r)
            im = pow(-m, mod - 2, mod)
            for i in range(m):
                r[i] = r[i] * im % mod
            for i in range(m):
                t[i] = (t[i] + r[i]) % mod
            t = [t[-1]] + t[:-1]

            t += [0] * m
            butterfly(t)
            g_fft = g + [0] * (2 * m - len(g))
            butterfly(g_fft)
            for i in range(2 * m):
                t[i] = t[i] * g_fft[i] % mod
            butterfly_inv(t)
            t = t[:m]
            i2m = pow(2 * m, mod - 2, mod)
            for i in range(m):
                t[i] = t[i] * i2m % mod
        
            v = a[m:min(n, 2 * m)]
            v += [0] * (m - len(v))
            t = [0] * (m - 1) + t + [0]
            t=integrate(t)
            for i in range(m):
                v[i] = (v[i] - t[m + i]) % mod

            v += [0] * m
            butterfly(v)
            for i in range(2 * m):
                v[i] = v[i] * f_fft[i] % mod
            butterfly_inv(v)
            v = v[:m]
            i2m = pow(2 * m, mod - 2, mod)
            for i in range(m):
                v[i] = v[i] * i2m % mod
            
            for i in range(min(n - m, m)):
                a[m + i] = v[i]
            
            m *= 2
        return a

    # power
    def power(a, k):
        n = len(a)
        assert n > 0
        if k == 0:
            return [1] + [0] * (n - 1)
        l = 0
        while l < len(a) and not a[l]:
            l += 1
        if l * k >= n:
            return [0] * n
        ic = pow(a[l], mod - 2, mod)
        pc = pow(a[l], k, mod)
        a = log([a[i] * ic % mod for i in range(l, len(a))])
        for i in range(len(a)):
            a[i] = a[i] * k % mod
        a = exp(a)
        for i in range(len(a)):
            a[i] = a[i] * pc % mod
        a = [0] * (l * k) + a[:n - l * k]
        return a

    # square root
    def sqrt(a):
        if len(a) == 0:
            return []
        if a[0] == 0:
            for d in range(1, len(a)):
                if a[d]:
                    if d & 1:
                        return None
                    if len(a) - 1 < d // 2:
                        break
                    res = sqrt(a[d:] + [0] * (d//2))
                    if res == None:
                        return None
                    res = [0] * (d // 2) + res
                    return res
            return [0] * len(a)
        
        # 先求第一项能否找到对应的数字，接下来再凑
        sqr = Tonelli_Shanks(a[0], mod)
        if sqr == None:
            return None
        T = [0] * (len(a))
        T[0] = sqr
        res = T.copy()
        T[0] = pow(sqr, mod - 2, mod) #T:res^{-1}
        m = 1
        two_inv = (mod + 1) // 2
        F = [sqr]
        while m <= len(a) - 1:
            for i in range(m):
                F[i] *= F[i]
                F[i] %= mod
            butterfly_inv(F)
            iz = pow(m, mod-2, mod)
            for i in range(m):
                F[i] = F[i] * iz % mod
            delta = [0] * (2 * m)
            for i in range(m):
                delta[i + m] = F[i] - a[i] - (a[i + m] if i + m < len(a) else 0)
            butterfly(delta)
            G = [0] * (2 * m)
            for i in range(m):
                G[i] = T[i]
            butterfly(G)
            for i in range(2 * m):
                delta[i] *= G[i]
                delta[i] %= mod
            butterfly_inv(delta)
            iz = pow(2 * m, mod - 2, mod)
            for i in range(2 * m):
                delta[i] = delta[i] * iz % mod
            for i in range(m, min(2 * m, len(a))):
                res[i] = -delta[i] * two_inv%mod
                res[i] %= mod
            if 2 * m > len(a) - 1:
                break
            F = res[:2 * m]
            butterfly(F)
            eps = [F[i] * G[i] % mod for i in range(2 * m)]
            butterfly_inv(eps)
            for i in range(m):
                eps[i] = 0
            iz = pow(2 * m, mod - 2, mod)
            for i in range(m,2*m):
                eps[i] = eps[i] * iz % mod
            butterfly(eps)
            for i in range(2 * m):
                eps[i] *= G[i]
                eps[i] %= mod
            butterfly_inv(eps)
            for i in range(m, 2 * m):
                T[i] = -eps[i]*iz
                T[i] %= mod
            iz = iz * iz % mod

            m <<= 1
        return res
