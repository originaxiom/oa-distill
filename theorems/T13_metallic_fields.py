#!/usr/bin/env python3
"""T13 -- THE METALLIC FIBRE FIELDS, EXACTLY.  Statement: for the metallic rules sigma_m: a -> a^m b, b -> a, the
fixed points of the squared Fricke action on the Markoff surface (kappa = 0), up to the SL(2) lift signs, have
x-coordinates whose minimal polynomials over Q are found by integer relations at 60 digits.  m = 1 gives x^2 - 3x + 3
(discriminant -3); the other m give the fields of the metallic bundles' fibre characters.  Comparison class: m = 1..5,
four sign twists.  This is the rule-side statement of which member of the class carries conductor 3."""
import mpmath as mp
from _common import say, write
mp.mp.dps = 60
def cheb(n, x):
    S = [mp.mpf(1), x]
    for k in range(2, n + 2): S.append(x * S[-1] - S[-2])
    return S
def Fm(p, m):
    x, y, z = p; S = cheb(m + 1, x); Sm2 = S[m - 2] if m >= 2 else mp.mpf(0)
    return (S[m - 1] * z - Sm2 * y, x, S[m] * z - S[m - 1] * y)
def G(v, m, e):
    p = tuple(v); q = Fm(Fm(p, m), m)
    return mp.matrix([q[0] - e[0]*p[0], q[1] - e[1]*p[1], q[2] - e[2]*p[2], p[0]**2 + p[1]**2 + p[2]**2 - p[0]*p[1]*p[2]])
def newton(v, m, e, iters=60):
    for it in range(iters):
        g = G(v, m, e)
        if mp.norm(g) < mp.mpf(10) ** (-50): return v, True
        J = mp.matrix(4, 3); h = mp.mpf(10) ** (-25)
        for k in range(3):
            dv = mp.matrix(3, 1); dv[k] = h
            gk = G(v + dv, m, e)
            for i in range(4): J[i, k] = (gk[i] - g[i]) / h
        # least squares step
        try: step = mp.qr_solve(J, g)[0]
        except Exception: return v, False
        v = v - step
    return v, mp.norm(G(v, m, e)) < mp.mpf(10) ** (-40)
from snappy import pari
pari.set_real_precision(60)
def minpoly(val, maxdeg=12):
    """integer relation among 1, val, ..., val^d by PARI lindep on the complex vector at 60 digits."""
    w = pari(f'{mp.nstr(val.real, 58)}+{mp.nstr(val.imag, 58)}*I')
    for d in range(1, maxdeg + 1):
        rel = pari.lindep([w ** i for i in range(d + 1)])
        if len(rel) < d + 1: continue
        rel = [int(rel[i]) for i in range(d + 1)]
        if rel[d] == 0 or max(abs(c) for c in rel) > 10**8: continue
        if rel[d] < 0: rel = [-c for c in rel]
        if abs(sum(c * val ** i for i, c in enumerate(rel))) < mp.mpf(10) ** (-40):
            return rel
    return None
import random
random.seed(1)
out = {}
for m in range(1, 6):
    found = {}
    for e in [(1, 1, 1), (-1, -1, 1), (-1, 1, -1), (1, -1, -1)]:
        for t in range(60):
            v = mp.matrix([mp.mpc(random.gauss(0, 2), random.gauss(0, 2)) for _ in range(3)])
            v, ok = newton(v, m, e)
            if not ok or min(abs(v[i]) for i in range(3)) < mp.mpf('1e-8'): continue
            rel = minpoly(v[0])
            key = str(rel) if rel else 'unrecognised'
            found.setdefault(key, set()).add(str(e))
    out[m] = {k: sorted(s) for k, s in found.items()}
    say(f"m={m}: x-minimal polynomials (coefficients low->high) at nondegenerate fixed Markoff points: {out[m]}")
# discriminants / fields of the quadratic ones, degree summary
summary = {}
for m, d in out.items():
    degs = sorted({len(eval(k)) - 1 for k in d if k != 'unrecognised'})
    summary[m] = degs
    say(f"m={m}: degrees of recognised minimal polynomials: {degs}")
write('T13', dict(minpolys={str(m): d for m, d in out.items()}, degrees={str(m): d for m, d in summary.items()}))
