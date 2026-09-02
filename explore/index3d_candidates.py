#!/usr/bin/env python3
"""Exploration for T16: the 3d index of an ideal triangulation from the tetrahedron index
   I_D(m,e)(q) = sum_{n >= max(0,-e)} (-1)^n q^{n(n+1)/2 - (n+e/2) m} / ((q)_n (q)_{n+e}),
glued with a candidate rule (m_i, e_i) = (k.(G_p1 - G_p3), k.(G_p2 - G_p3))_i over integer edge weights k (one edge fixed
to 0), for the six permutations p of the quad columns (z, z', z'').  Acceptance: a candidate whose total index is the
same q-series for several triangulations of m004 and differs for m003.  Series in x = q^{1/2}, truncated at x-degree D."""
import sys, itertools, functools, numpy as np, snappy
D = int(sys.argv[1]) if len(sys.argv) > 1 else 12          # x-degree (q^{D/2})
W = 40                                                       # allowed negative x-exponent window
L = 2 * W + D + 1                                            # exponents -W .. D+W (headroom for shifts)
def zero(): return np.zeros(L, dtype=object)
def mono(exp, c=1):
    s = zero()
    if -W <= exp <= D + W: s[exp + W] = c
    return s
def mul(a, b):
    out = zero()
    ia = np.nonzero(a)[0]; ib = np.nonzero(b)[0]
    for i in ia:
        for j in ib:
            e = i + j - W                             # index of the product's exponent (exp = e - W)
            if 0 <= e <= D + 2 * W: out[e] += a[i] * b[j]
    return out
def qpoch_inv(n):
    """1/(q)_n as a series in x (q = x^2), nonnegative exponents only."""
    s = mono(0)
    for i in range(1, n + 1):
        # 1/(1 - q^i) = sum_j q^{ij}
        g = zero()
        j = 0
        while 2 * i * j <= D + W: g[2 * i * j + W] = 1; j += 1
        s = mul(s, g)
    return s
@functools.lru_cache(maxsize=None)
def tet_index(m, e):
    s = zero(); n = max(0, -e)
    while True:
        lead = n * (n + 1) - 2 * n * m - e * m          # x-exponent of the leading monomial
        if lead > D + 2 * W + 60 and n > abs(m) + abs(e) + 5: break   # crude: beyond window for this and larger n
        if n > 60: break
        if lead <= D:
            term = mul(mul(qpoch_inv(n), qpoch_inv(n + e)), mono(lead, (-1) ** n))
            s = s + term
        n += 1
    s[D + W + 1:] = 0                                          # report only exponents <= D
    return s
def index_of(M, perm, kmax):
    G = M.gluing_equations(); N = M.num_tetrahedra()
    rows = [list(G[i]) for i in range(N)]                    # the N edge rows (drop cusp rows)
    def col(i, t): return [rows[E][3 * i + t] for E in range(N)]
    total = zero()
    p1, p2, p3 = perm
    for k in itertools.product(range(-kmax, kmax + 1), repeat=N - 1):
        kk = list(k) + [0]
        prod = mono(0)
        ok = True
        for i in range(N):
            m = sum(kk[E] * (rows[E][3*i + p1] - rows[E][3*i + p3]) for E in range(N))
            e = sum(kk[E] * (rows[E][3*i + p2] - rows[E][3*i + p3]) for E in range(N))
            t = tet_index(int(m), int(e))
            if not t.any(): ok = False; break
            prod = mul(prod, t)
            if not prod.any(): ok = False; break
        if ok: total = total + prod
    return total
def show(s):
    nz = [(i - W, int(s[i])) for i in np.nonzero(s)[0]]
    return ' '.join(f"{c}x^{e}" for e, c in nz) or '0'
if __name__ == '__main__':
    kmax = 5
    tris = {'m004': snappy.Manifold('m004')}
    base = snappy.Manifold('m004')
    seen = set()
    for t in range(40):
        R = base.copy(); R.randomize()
        key = (R.num_tetrahedra(), str(R.gluing_equations()))
        if R.num_tetrahedra() <= 4 and key not in seen and R.num_tetrahedra() > 2:
            seen.add(key); tris[f'm004_r{len(seen)}'] = R
        if len(seen) >= 3: break
    tris['m003'] = snappy.Manifold('m003')
    print({k: v.num_tetrahedra() for k, v in tris.items()})
    for perm in itertools.permutations(range(3)):
        res = {}
        for name, M in tris.items():
            res[name] = show(index_of(M, perm, kmax))
        same = len({v for k, v in res.items() if k.startswith('m004')}) == 1
        diff = res['m003'] != res['m004']
        print(f"perm {perm}: m004-invariant={same} m003-differs={diff}")
        for k, v in res.items(): print(f"    {k}: {v[:120]}")
        sys.stdout.flush()
