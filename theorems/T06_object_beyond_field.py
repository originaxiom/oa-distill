#!/usr/bin/env python3
"""T06 -- THE OBJECT BEYOND ITS FIELD.  Statement: within the 14-member Q(sqrt-3) class, H1, the cusp shape and the
number of covers by degree single out m004 (covers separate all 14); the sister m003's (1,q) closings differ from m004's
at every slope; the two 2-tetrahedron triangulations have different gluing matrices and different deformation curves,
each verified against a filled structure; hence the DGG theories T[m004], T[m003] differ.  Comparison class: the 14."""
import cmath, sympy as sp, snappy
from _common import say, write
FAM = ["m003","m004","m202","m203","m206","m207","m208","m410","m412","s118","s119","s594","s595","s596"]
rows = {}
for n in FAM:
    M = snappy.Manifold(n); sh = M.cusp_info('shape')[0]
    rows[n] = dict(H1=str(M.homology()), shape=(round(float(sh.real), 6), round(float(sh.imag), 6)),
                   covers=[len(M.covers(d)) for d in (2, 3, 4, 5, 6)], amph=bool(M.symmetry_group().is_amphicheiral()))
    say(f"  {n:5s} H1={rows[n]['H1']:12s} shape={rows[n]['shape']} covers(2..6)={rows[n]['covers']} amph={rows[n]['amph']}")
sep = [k for k in ('H1', 'shape', 'covers') if sum(1 for n in FAM if rows[n][k] == rows['m004'][k]) == 1]
say(f"separators of m004 in the class: {sep}; covers-by-degree distinct across all 14: {len({str(r['covers']) for r in rows.values()}) == 14}; chiral members: {sum(1 for r in rows.values() if not r['amph'])}")
fills = {}
for n in ('m004', 'm003'):
    fills[n] = {}
    for q in range(3, 7):
        M = snappy.Manifold(n); M.chern_simons(); M.dehn_fill((1, q)); M.volume()
        fills[n][q] = (round(float(M.volume()), 5), round(float(M.chern_simons()), 5), str(M.homology()))
    say(f"  {n} (1,q) q=3..6: {fills[n]}")
z1, z2, m, l = sp.symbols('z1 z2 m l'); zs = [z1, z2]
curves = {}
for n in ('m004', 'm003'):
    M = snappy.Manifold(n); rect = M.gluing_equations('rect'); edges, mer, lon = rect[:-2], rect[-2], rect[-1]
    prod = lambda a, b: sp.Mul(*[zs[i]**a[i] * (1 - zs[i])**b[i] for i in range(2)])
    E = [sp.numer(sp.together(prod(e[0], e[1]) - e[2])) for e in edges]
    Mq = sp.numer(sp.together(prod(mer[0], mer[1]) - mer[2]*m)); Lq = sp.numer(sp.together(prod(lon[0], lon[1]) - lon[2]*l))
    G = sp.groebner(E + [Mq, Lq], z1, z2, m, l, order='lex')
    A = sp.factor(sp.gcd_list([g for g in G.exprs if not (g.has(z1) or g.has(z2))]))
    N = snappy.Manifold(n); N.dehn_fill((1, 5)); zz = [complex(z) for z in N.tetrahedra_shapes('rect')]
    hol = N.cusp_info()[0]['holonomies']; em, el = cmath.exp(complex(hol[0])), cmath.exp(complex(hol[1]))
    val = abs(complex(A.subs({m: em, l: el})))
    curves[n] = (str(A), val)
    gm = str(M.gluing_equations()).replace('\n', ' ')
    say(f"  {n} gluing matrix {gm}\n     curve {A}\n     |curve(exp H_m, exp H_l)| at (1,5) = {val:.2e}")
write('T06', dict(rows=rows, separators=sep, n_chiral=sum(1 for r in rows.values() if not r['amph']),
                  fillings=fills, curves={n: c[0] for n, c in curves.items()}, curve_residuals={n: c[1] for n, c in curves.items()}))
