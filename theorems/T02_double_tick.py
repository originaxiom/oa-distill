#!/usr/bin/env python3
"""T02 -- THE DOUBLE TICK.  Statement: the incidence matrix F of sigma has det -1, so sigma is orientation-reversing on
the once-punctured torus; the mapping torus of sigma is the Gieseking manifold (non-orientable); the orientation double
cover of the Gieseking manifold is m004, the mapping torus of sigma^2 = A = [[2,1],[1,1]]; F commutes with A, so sigma acts
fibrewise on m004 as an orientation-reversing symmetry.  Comparison class: the one-tick and two-tick mapping tori."""
import itertools, sympy as sp, snappy
from _common import say, write
F = sp.Matrix([[1, 1], [1, 0]]); A = F * F
G = snappy.Manifold('m000'); C = G.orientation_cover(); m004 = snappy.Manifold('m004')
say(f"det F = {F.det()}; F^2 = {A.tolist()}; F A == A F: {F*A == A*F}")
say(f"m000 orientable: {G.is_orientable()}, tetrahedra {G.num_tetrahedra()}, volume {float(G.volume()):.9f}")
say(f"orientation cover of m000 is isometric to m004: {C.is_isometric_to(m004)}; volume ratio {float(m004.volume())/float(G.volume()):.9f}")
bund = {}
for name in ('b--R', 'b-+R'):
    B = snappy.Manifold(name); bund[name] = (bool(B.is_orientable()), bool(B.is_isometric_to(G)))
say(f"one-tick bundles {bund}")
cent = [(a, b, c, d) for a, b, c, d in itertools.product(range(-5, 6), repeat=4)
        if a*d - b*c == -1 and sp.Matrix([[a, b], [c, d]]) * A == A * sp.Matrix([[a, b], [c, d]])]
say(f"det -1 elements commuting with A (|entries|<=5): {cent}")
write('T02', dict(detF=int(F.det()), A=[[int(v) for v in r] for r in A.tolist()], commute=bool(F*A == A*F), m000_orientable=bool(G.is_orientable()),
                  cover_is_m004=bool(C.is_isometric_to(m004)), one_tick_bundles=bund, det_minus1_centralisers=cent))
