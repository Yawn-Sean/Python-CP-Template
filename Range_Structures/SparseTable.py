class SparseTable:
    def __init__(self, data, merge_method):
        self.note = [0] * (len(data) + 1)
        self.merge_method = merge_method
        l, r, v = 1, 2, 0
        while True:
            for i in range(l, r):
                if i >= len(self.note):
                    break
                self.note[i] = v
            else:
                l *= 2
                r *= 2
                v += 1
                continue
            break
        self.ST = [[0] * len(data) for _ in range(self.note[-1]+1)]
        self.ST[0] = data
        for i in range(1, len(self.ST)):
            for j in range(len(data) - (1 << i) + 1):
                self.ST[i][j] = merge_method(self.ST[i-1][j], self.ST[i-1][j + (1 << (i-1))])

    def query(self, l, r):
        pos = self.note[r-l+1]
        return self.merge_method(self.ST[pos][l], self.ST[pos][r - (1 << pos) + 1])