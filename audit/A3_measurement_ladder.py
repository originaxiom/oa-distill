#!/usr/bin/env python3
"""A3 -- the measurement ladder (record's B874): centralizers in e6 of subtori of the charge torus C = <x8,x14,x16,x22>
(the 2T-invariants of A2).  dim Cent(S) = 78 - rank of the stacked ad(x), x in S.  Also the structure of Cent(C):
derived algebra and centre dimensions.  CHOICE INVENTORY: which coordinate subtorus is 'measured'."""
import numpy as np, itertools, os
HERE = os.path.dirname(os.path.abspath(__file__))
d = np.load(os.path.join(HERE, 'e6_data.npz')); ad = d['ad']; n = 78
a2 = np.load(os.path.join(HERE, 'a2_data.npz')); inv = a2['inv']; deg = list(a2['degrees'])
X = {int(m): inv[i] for i, m in enumerate(deg)}
def opmat(v): return sum(v[i] * ad[i] for i in range(n) if abs(v[i]) > 0)
def cent(vectors):
    M = np.vstack([opmat(v) for v in vectors]); s = np.linalg.svd(M, compute_uv=False)
    r = int(np.sum(s > 1e-8 * s[0])); return n - r, M
def kernel(M):
    u, s, vh = np.linalg.svd(M); r = int(np.sum(s > 1e-8 * s[0])); return vh[r:].conj().T
print("centralizer dimensions of coordinate subtori (record B874: single 30 or 12; C = 12):")
for k in range(1, 5):
    for sub in itertools.combinations(sorted(X), k):
        dim, _ = cent([X[m] for m in sub]); print(f"  {sub}: {dim}")
dimC, M = cent([X[m] for m in X]); K = kernel(M)                      # basis of Cent(C)
# structure of Cent(C): derived algebra = span of brackets; centre = elements commuting with all of Cent(C)
B = np.column_stack([opmat(K[:, i]) @ K[:, j] for i in range(K.shape[1]) for j in range(K.shape[1])])
der = np.linalg.matrix_rank(B, tol=1e-7)
Z = kernel(np.vstack([opmat(K[:, i]) @ K for i in range(K.shape[1])]))   # x in Cent(C) (coords) with [x, Cent(C)] = 0
print(f"Cent(C): dim {dimC}, derived algebra dim {der}, centre dim {Z.shape[1]}  (record: 12 = su(2,1)+C over R, i.e. sl3 (8) + centre 4 over C)")
np.savez(os.path.join(HERE, 'a3_data.npz'), centC=K)
