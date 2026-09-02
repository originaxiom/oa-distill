#!/usr/bin/env python3
"""T03 -- SYMMETRY AND MIRROR.  Statement: Sym(m004) has order 8 (D4), four of its isometries reverse orientation
(cusp-map determinant -1), and the group is amphichiral; A -> A^-1 is realised by an orientation-preserving fibre map
(the bundle is invertible); the discrete faithful character (tr a, tr b, tr ab) lies in Q(sqrt-3)^3 and its complex
conjugate is realised by explicit endomorphisms a->u, b->v of pi_1 (relator -> +-I): the mirror acts on characters as
complex conjugation.  Comparison class: the isometries of m004; words of length <= 7."""
import itertools, sympy as sp, snappy, numpy as np
from _common import say, write
M = snappy.Manifold('m004'); S = M.symmetry_group()
isos = M.is_isometric_to(M, return_isometries=True)
dets = []
for iso in isos:
    C = iso.cusp_maps()[0]; dets.append(int(C[0, 0]*C[1, 1] - C[0, 1]*C[1, 0]))
say(f"|Sym| = {S.order()}, amphichiral = {S.is_amphicheiral()}, isometries {len(isos)}, orientation-reversing {dets.count(-1)}")
A = sp.Matrix([[2, 1], [1, 1]])
invp = [(a, b, c, d) for a, b, c, d in itertools.product(range(-5, 6), repeat=4)
        if a*d - b*c == 1 and sp.Matrix([[a, b], [c, d]]) * A == A.inv() * sp.Matrix([[a, b], [c, d]])]
say(f"A -> A^-1 by det +1 maps (invertible bundle): {invp[:3]}")
Gp = M.fundamental_group(); rel = Gp.relators()[0]; gens = Gp.generators()
a = np.array([[complex(Gp.SL2C('a')[i, j]) for j in range(2)] for i in range(2)])
b = np.array([[complex(Gp.SL2C('b')[i, j]) for j in range(2)] for i in range(2)])
mats = {'a': a, 'b': b, 'A': np.linalg.inv(a), 'B': np.linalg.inv(b)}
def ev(w):
    m = np.eye(2, dtype=complex)
    for ch in w: m = m @ mats[ch]
    return m
tr = lambda m: m[0, 0] + m[1, 1]
x, y, z = tr(a), tr(b), tr(a @ b)
target = (x.conjugate(), y.conjugate(), z.conjugate())
inK = lambda v: abs(2*v.real - round(2*v.real)) < 1e-9 and abs(2*v.imag/3**.5 - round(2*v.imag/3**.5)) < 1e-9
say(f"character (tr a, tr b, tr ab) = ({x:.6f}, {y:.6f}, {z:.6f}); all in (1/2)Z[sqrt-3]: {all(map(inK, (x, y, z)))}")
words = ['']
for L in range(1, 8):
    for t in itertools.product('abAB', repeat=L):
        w = ''.join(t)
        if any(w[i] == w[i+1].swapcase() for i in range(L-1)): continue
        words.append(w)
W = {w: ev(w) for w in words if w}
cu = [w for w, m in W.items() if abs(tr(m) - target[0]) < 1e-8]
cv = [w for w, m in W.items() if abs(tr(m) - target[1]) < 1e-8]
hits = []
for u in cu:
    for v in cv:
        if abs(tr(W[u] @ W[v]) - target[2]) > 1e-8: continue
        sub = {'a': W[u], 'b': W[v], 'A': np.linalg.inv(W[u]), 'B': np.linalg.inv(W[v])}
        r = np.eye(2, dtype=complex)
        for ch in rel: r = r @ sub[ch]
        if np.allclose(r, np.eye(2), atol=1e-7) or np.allclose(r, -np.eye(2), atol=1e-7):
            hits.append((u, v, '+I' if np.allclose(r, np.eye(2), atol=1e-7) else '-I'))
say(f"relator {rel}; endomorphisms realising the conjugate character (words <= 7): {len(hits)}; first: {hits[:3]}")
write('T03', dict(sym_order=S.order(), amphichiral=bool(S.is_amphicheiral()), n_isometries=len(isos), n_orientation_reversing=dets.count(-1),
                  invertible_examples=invp[:3], character=[[x.real, x.imag], [y.real, y.imag], [z.real, z.imag]],
                  n_conjugate_realisations=len(hits), first=hits[:3]))
