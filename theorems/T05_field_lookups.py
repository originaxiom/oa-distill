#!/usr/bin/env python3
"""T05 -- THE FIELD LOOKUPS, STATED AS LOOKUPS.  Statement: the route (field discriminant -> conductor N -> SL(2,Z/N)
-> finite subgroup of SU(2)? -> McKay label) emits a label only for N in {1, 3, 5} (checked N <= 24); imaginary quadratic
conductors reach only N = 3, so on hyperbolic manifolds the route says E6 or nothing.  Among the first 1200 orientable
cusped census manifolds, 14 have shape field Q(sqrt-3); their volumes are 12, 24, 30 times the covolume of PSL(2,O_-3).
Comparison class: all conductors; the census."""
import itertools, math
from fractions import Fraction
import snappy
from snappy import pari
from _common import say, write
def rat(x, maxden=64, tol=1e-9):
    f = Fraction(x).limit_denominator(maxden); return f if abs(float(f) - x) < tol else None
def quad_disc(z):
    b, c = rat(-2*z.real), rat(abs(z)**2)
    if b is None or c is None: return None
    d = b*b - 4*c
    if d >= 0: return None
    n = d.numerator * d.denominator
    for s in range(2, 60):
        while n % (s*s) == 0: n //= s*s
    return n
def mul(x, y, N): return ((x[0]*y[0]+x[1]*y[2]) % N, (x[0]*y[1]+x[1]*y[3]) % N, (x[2]*y[0]+x[3]*y[2]) % N, (x[2]*y[1]+x[3]*y[3]) % N)
def su2_label(N):
    if N == 1: return 'A_0'
    G = [(a, b, c, d) for a, b, c, d in itertools.product(range(N), repeat=4) if (a*d - b*c) % N == 1]
    I = (1, 0, 0, 1)
    def order(g):
        h, k = g, 1
        while h != I: h = mul(h, g, N); k += 1
        return k
    ords = [order(g) for g in G]; n = len(G); inv = ords.count(2); mo = max(ords); spec = sorted(set(ords))
    if inv > 1: return None
    if mo == n // 2: return f'D_{n//4+2}'
    if n == 24 and spec == [1, 2, 3, 4, 6]: return 'E6'
    if n == 48 and spec == [1, 2, 3, 4, 6, 8]: return 'E7'
    if n == 120 and spec == [1, 2, 3, 4, 5, 6, 10]: return 'E8'
    return None
labels = {N: su2_label(N) for N in range(1, 25)}
say(f"McKay label by conductor N<=24: { {N: l for N, l in labels.items() if l} } ; all others: none")
imag = [3, 4, 7, 8, 11, 15, 19, 20, 23, 24]
say(f"imaginary quadratic conductors <= 24 that get a label: {[N for N in imag if labels[N]]}")
hits, scanned = [], 0
for M in snappy.OrientableCuspedCensus():
    scanned += 1
    if scanned > 1200: break
    sh = [complex(z) for z in M.tetrahedra_shapes('rect')]
    if {quad_disc(z) for z in sh} == {-3}: hits.append((M.name(), float(M.volume())))
Vtet = float(snappy.Manifold('m004').volume()) / 2
L2 = float(pari.lfun(pari.lfuncreate(-3), 2)); covol = 3**1.5 * (math.pi**2/6) * L2 / (4*math.pi**2)
idx = {n: round(v / covol, 6) for n, v in hits}
say(f"census manifolds sharing Q(sqrt-3) among the first {scanned-1}: {len(hits)}: {[n for n, _ in hits]}")
say(f"covolume PSL(2,O_-3) = {covol:.9f} = V_tet/{Vtet/covol:.6f}; indices: {idx}")
write('T05', dict(labels={N: l for N, l in labels.items() if l}, imaginary_labelled=[N for N in imag if labels[N]],
                  n_sharing=len(hits), sharing=[n for n, _ in hits], bianchi_indices=idx))
