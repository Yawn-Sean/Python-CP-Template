class OperationalFenwickTree:
    def __init__(self, n, op):
        self.n = n
        self.bit = [0] * n
        self.op = op

    def result(self, r):
        res = 0
        while r >= 0:
            res = self.op(res, self.bit[r])
            r = (r & (r + 1)) - 1
        return res

    def update(self, idx, val):
        while idx < self.n:
            self.bit[idx] = self.op(self.bit[idx], val)
            idx = idx | (idx + 1)
