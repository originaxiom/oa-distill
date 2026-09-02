#!/usr/bin/env python3
"""A1 -- build the Lie algebra e6 from its root system (Frenkel-Kac cocycle), as 78 x 78 adjoint matrices over Q,
and verify: 72 roots, Jacobi identity, Killing form nondegenerate, rank 6.  Everything downstream of E6 in the record's
chain is recomputed on this algebra.  Saved to e6_data.npz."""
import numpy as np, itertools, os
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__))
# E6 Cartan matrix (Bourbaki: 1-3-4-5-6 chain, 2 attached to 4)
C = np.array([[2,0,-1,0,0,0],[0,2,0,-1,0,0],[-1,0,2,-1,0,0],[0,-1,-1,2,-1,0],[0,0,0,-1,2,-1],[0,0,0,0,-1,2]])
# simply laced: (alpha_i, alpha_j) = C_ij.  Positive roots by closure in simple-root coordinates.
simple = [tuple(int(i == j) for j in range(6)) for i in range(6)]
def ip(a, b): return int(np.array(a) @ C @ np.array(b))
pos = set(simple); frontier = list(simple)
while frontier:
    new = []
    for r in frontier:
        for s in simple:
            t = tuple(x + y for x, y in zip(r, s))
            # t is a root iff (r, s) = -1 (string length via reflection): for simply laced, r+s root <=> (r,s) = -1
            if ip(r, s) == -1 and t not in pos:
                pos.add(t); new.append(t)
    frontier = new
pos = sorted(pos, key=lambda r: (sum(r), r))
roots = pos + [tuple(-x for x in r) for r in pos]
print("positive roots:", len(pos), " total:", len(roots))
assert len(pos) == 36
idx = {r: i for i, r in enumerate(roots)}
# Frenkel-Kac cocycle: eps(alpha_i, alpha_j) = -1 if i < j and C_ij = -1 ... standard: eps bimultiplicative with
# eps(a_i, a_j) = (-1)^{(a_i,a_j)} for i < j, = (-1)^{(a_i,a_i)/2} = -1 for i = j, = 1 for i > j.
import sys
ORIENT = int(sys.argv[1]) if len(sys.argv) > 1 else 1        # 1: sign on i<j ; -1: sign on i>j
HSIGN = int(sys.argv[2]) if len(sys.argv) > 2 else -1        # [e_r, e_-r] = HSIGN * h_r
def eps(a, b):
    s = 1
    for i in range(6):
        for j in range(6):
            if a[i] and b[j]:
                lower = (i < j) if ORIENT == 1 else (i > j)
                if lower: e = -1 if (C[i, j] * a[i] * b[j]) % 2 else 1
                elif i == j: e = -1 if (a[i] * b[j]) % 2 else 1     # eps(a_i, a_i) = (-1)^{(a_i,a_i)/2} = -1
                else: e = 1
                s *= e
    return s
# basis: h_1..h_6 (coroots), then e_r for r in roots.  Brackets:
# [h_i, e_r] = (a_i, r) e_r ; [e_r, e_-r] = h_r = sum r_i h_i ; [e_r, e_s] = eps(r,s) e_{r+s} if r+s root else 0
n = 6 + 72
def E(i): v = np.zeros(n, dtype=object); v[i] = Fraction(1); return v
def bracket(u, v):
    out = np.zeros(n, dtype=object)
    for i in range(n):
        if u[i] == 0: continue
        for j in range(n):
            if v[j] == 0: continue
            c = u[i] * v[j]
            if i < 6 and j < 6: continue
            if i < 6 and j >= 6:
                r = roots[j - 6]; out[j] += c * ip(simple[i], r)
            elif i >= 6 and j < 6:
                r = roots[i - 6]; out[i] -= c * ip(simple[j], r)
            else:
                r, s = roots[i - 6], roots[j - 6]
                t = tuple(x + y for x, y in zip(r, s))
                if all(x == 0 for x in t):
                    for k in range(6): out[k] += HSIGN * c * r[k]   # [e_r, e_-r] = HSIGN * h_r
                elif t in idx:
                    out[6 + idx[t]] += c * eps(r, s)
    return out
# adjoint matrices
ad = []
for i in range(n):
    Mx = np.zeros((n, n), dtype=object)
    for j in range(n):
        col = bracket(E(i), E(j))
        for k in range(n): Mx[k, j] = col[k]
    ad.append(Mx)
# Jacobi on random triples (exact)
import random
random.seed(0); bad = 0
for _ in range(2000):
    i, j, k = random.sample(range(n), 3)
    lhs = bracket(E(i), bracket(E(j), E(k))) + bracket(E(j), bracket(E(k), E(i))) + bracket(E(k), bracket(E(i), E(j)))
    if any(x != 0 for x in lhs): bad += 1
print("Jacobi failures on 2000 random triples:", bad, " (ORIENT, HSIGN) =", (ORIENT, HSIGN))
# Killing form and rank
adf = np.array([[float(x) for x in row] for Mx in ad for row in Mx]).reshape(n, n, n)
K = np.array([[np.trace(adf[i] @ adf[j]) for j in range(n)] for i in range(n)])
print("Killing form rank:", np.linalg.matrix_rank(K), " ad-representation: [ad_i, ad_j] = ad_[i,j] check on samples:",
      all(np.allclose(adf[i] @ adf[j] - adf[j] @ adf[i], sum(float(c) * adf[k] for k, c in enumerate(bracket(E(i), E(j))) if c != 0) if any(c != 0 for c in bracket(E(i), E(j))) else np.zeros((n, n))) for i, j in [(0, 7), (6, 42), (10, 60), (3, 30)]))
np.savez(os.path.join(HERE, 'e6_data.npz'), ad=adf, K=K, roots=np.array(roots), C=C)
print("e6 built and saved; Jacobi", "OK" if bad == 0 else "FAIL")
