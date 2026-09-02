#!/usr/bin/env python3
"""A5 -- robustness checks on the A4 discrepancy: (1) can the invariants be made real; (2) structure of the generic
joint centralizer in stage 2; (3) with the UNTUNED first charge x8, do the other pencils (x14,x22), (x16,x22) in Cent(x8)
ever reach a 14-dimensional joint centralizer?"""
import numpy as np, scipy.linalg as la, os, itertools
HERE = os.path.dirname(os.path.abspath(__file__))
d = np.load(os.path.join(HERE, 'e6_data.npz')); ad = d['ad']; n = 78
a2 = np.load(os.path.join(HERE, 'a2_data.npz')); inv = a2['inv']; deg = list(a2['degrees'])
X = {}
for i, m in enumerate(deg):
    v = inv[i]; ph = v[np.argmax(np.abs(v))]; v = v / (ph / abs(ph)); X[int(m)] = v
    print(f"x{m}: max |Im| after phase rotation = {np.abs(v.imag).max():.2e}  (real vector up to phase: {np.abs(v.imag).max() < 1e-9})")
X = {m: v.real if np.abs(v.imag).max() < 1e-9 else v for m, v in X.items()}
def opmat(v): return sum(v[i] * ad[i] for i in range(n) if abs(v[i]) > 0)
def kernel(M, rel=1e-8):
    u, s, vh = np.linalg.svd(M); r = int(np.sum(s > rel * s[0])); return vh[r:].conj().T
def nullity(M, rel=1e-8):
    s = np.linalg.svd(M, compute_uv=False); return int(np.sum(s < rel * s[0]))
def structure(J):
    """dim, derived-algebra dim, centre dim of the subalgebra spanned by the columns of J (orthonormalised first)."""
    Q, _ = np.linalg.qr(J); k = Q.shape[1]
    B = np.column_stack([opmat(Q[:, i]) @ Q[:, j] for i in range(k) for j in range(k)])
    closure = np.linalg.norm(B - Q @ (Q.conj().T @ B)) / max(1e-12, np.linalg.norm(B))   # must be ~0 for a subalgebra
    Bp = Q.conj().T @ B
    sB = np.linalg.svd(Bp, compute_uv=False); der = int(np.sum(sB > 1e-7 * sB[0])) if sB[0] > 1e-12 else 0
    if closure > 1e-6: der = f"{der} (NOT CLOSED: {closure:.1e})"
    Z = kernel(np.vstack([opmat(Q[:, i]) @ Q for i in range(k)]), rel=1e-7)
    return k, der, Z.shape[1]
def rank_drops(A, B, generic_rank, seeds=(11, 22), bound=1e5):
    res = []
    for sd in seeds:
        rg = np.random.default_rng(sd)
        R = rg.normal(size=(generic_rank, A.shape[0])) + 1j * rg.normal(size=(generic_rank, A.shape[0])); S = rg.normal(size=(A.shape[1], generic_rank)) + 1j * rg.normal(size=(A.shape[1], generic_rank))
        w = la.eigvals(R @ A @ S, -R @ B @ S); res.append([t for t in w if np.isfinite(t) and abs(t) < bound])
    if not res[0] or not res[1]: return []
    out = []
    for t0 in res[0]:
        if min(abs(t0 - u) for u in res[1]) < 1e-5 * max(1, abs(t0)) and all(abs(t0 - c) > 1e-4 * max(1, abs(t0)) for c in out): out.append(t0)
    return out
A = {m: opmat(X[m]) for m in X}
# (1)/(2): stage 1 with real invariants; structure of generic stage-2 centralizer
ts = rank_drops(A[8], A[16], n - nullity(A[8] + 0.3137 * A[16]))
print("stage-1 jump points with real invariants:", [np.round(t, 4) for t in sorted(ts, key=abs)])
for t0 in sorted(ts, key=abs)[:3]:
    K = kernel(A[8] + t0 * A[16]); print(f"  t*={t0:.4f}: Cent(x1) dim/derived/centre = {structure(K)}")
    B14, B16 = A[14] @ K, A[16] @ K
    Jg = K @ kernel(B14 + 0.2718 * B16); print(f"     generic joint centralizer (u generic): dim/derived/centre = {structure(Jg)}")
# (3): untuned x8 with other pencils
K8 = kernel(A[8]); print(f"\nuntuned x1 = x8: dim Cent = {K8.shape[1]}")
for (p, q) in [(14, 22), (16, 22), (14, 16)]:
    Bp, Bq = A[p] @ K8, A[q] @ K8
    g = nullity(Bp + 0.2718 * Bq); us = rank_drops(Bp, Bq, K8.shape[1] - g)
    print(f"  pencil x{p} + u x{q} in Cent(x8): generic joint dim {g}; jump points: " +
          ", ".join(f"u={u0:.4f}->{structure(K8 @ kernel(Bp + u0 * Bq))}" for u0 in sorted(us, key=abs)[:6]) if us else f"  pencil x{p} + u x{q}: generic joint dim {g}; no finite jump points")
# also the full torus C and each 3-subtorus structures
for sub in [(8, 14, 16, 22), (8, 14, 22), (14, 16, 22)]:
    K = kernel(np.vstack([A[m] for m in sub])); print(f"Cent{sub}: dim/derived/centre = {structure(K)}")
