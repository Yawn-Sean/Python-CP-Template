block_size = ...
qs = [i for i in range(q)]
ans = [None] * q
qs.sort(key=lambda x: (right[x] // block_size, left[x] * (-1) ** (right[x] // block_size)))

l = 0
r = 0
cur = 0

def add(i, x=1):
    global cur

def delete(i):
    add(i, -1)

for i in qs:
    nl, nr = left[i], right[i]
    while r < nr:
        add(r)
        r += 1
    while nl < l:
        l -= 1
        add(l)
    while nr < r:
        r -= 1
        delete(r)
    while l < nl:
        delete(l)
        l += 1
    ans[i] = cur
