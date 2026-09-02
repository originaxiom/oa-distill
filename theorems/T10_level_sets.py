#!/usr/bin/env python3
"""T10 -- CHAIN AND OBJECT ARE ON DIFFERENT LEVEL SETS OF THE SAME INVARIANT.  Statement: the on-site Fibonacci
chain (T08) has Fricke invariant kappa = 4 + 4V^2 >= 4 for every potential strength V (exact); the object's fibre
character (T09) has kappa = 0 (parabolic puncture).  No V puts the chain on the object's level set.  Conversely, a
Dehn filling (p,q) of m004 moves the fibre boundary (the longitude of 4_1) off parabolicity: kappa(p,q) =
tr(longitude) + 2 = l + 1/l + 2; the chain's level sets need l real and > 0, l != 1.  We compute the longitude
eigenvalue for all hyperbolic fillings |p|,|q| <= 8 and report whether any is real.  Comparison class: all V; all
fillings in the box."""
import cmath, math, sympy as sp, snappy
from _common import say, write
E, V = sp.symbols('E V', real=True)
xx, yy = E - V, E + V; zz = sp.expand(xx * yy - 2)          # tr T_a, tr T_b, tr(T_a T_b) for T = [[E - v, -1],[1, 0]]
kap = sp.simplify(xx**2 + yy**2 + zz**2 - xx*yy*zz)
say(f"chain invariant kappa(E, V) = {sp.factor(kap)}  (E-independent; >= 4; = 4 iff V = 0)")
say("object's fibre character: kappa = 0 (T09).  Intersection of level sets: none for real V.")
M0 = snappy.Manifold('m004')
real_fills = []; table = []
for p in range(-8, 9):
    for q in range(0, 9):
        if q == 0 and p != 1: continue
        if math.gcd(abs(p), q) != 1: continue
        M = snappy.Manifold('m004'); M.dehn_fill((p, q))
        if 'positively' not in M.solution_type(): continue
        H = M.cusp_info()[0]['holonomies']; Hl = complex(H[1])
        l = cmath.exp(Hl); tr = l + 1 / l; kappa_pq = tr + 2
        is_real = abs(tr.imag) < 1e-9
        table.append(((p, q), round(tr.real, 6), round(tr.imag, 6), round(kappa_pq.real, 6)))
        if is_real: real_fills.append(((p, q), round(tr.real, 6)))
say(f"hyperbolic fillings in the box: {len(table)}; fillings with a REAL longitude trace: {real_fills}")
say("first rows ((p,q), Re tr l, Im tr l, Re kappa):"); [say("   ", r) for r in table[:10]]
# sanity: the longitude is the fibre boundary of the fibred knot 4_1: at the complete structure its trace is -2
G = M0.fundamental_group(); lon = G.peripheral_curves()[0][1]
Lm = G.SL2C(lon); trL = complex(Lm[0, 0] + Lm[1, 1])
say(f"complete structure: longitude word {lon}, trace {trL:.9f} (parabolic, fibre boundary: kappa = tr + 2 = {trL.real + 2:.2e})")
write('T10', dict(chain_kappa=str(sp.factor(kap)), n_fillings=len(table), real_longitude_fillings=real_fills,
                  table=table, longitude_trace_complete=[trL.real, trL.imag]))
