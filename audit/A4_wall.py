#!/usr/bin/env python3
"""A4 -- the two 'measurements' of the record (B874 / B892 / B950), recomputed on the audited e6.
Stage 1 (first measurement): the pencil x1(t) = x8 + t x16.  Note x16 is CENTRAL in Cent(x8) (A3: Cent(x8,x16) =
Cent(x8) = 30), so adding t x16 changes the centralizer: find the t where dim Cent(x1(t)) jumps (the record's 'tuned
combination at a cubic root t*').
Stage 2 (second measurement): for each such t*, the pencil y(u) = x14 + u x16 inside Cent(x1(t*)): find u where the
joint centralizer dim jumps (the record: to 14, algebra su(3)+su(2)+u(1)^3 after the B950 correction, at complex u).
Method for rank drops: generalised eigenvalues of two independent random compressions to the generic rank, confirmed
by singular values on the full matrix.  CHOICE INVENTORY: first charge x8; the tuning t*; the pencil line (x14, x16);
the branch of the algebraic t*, u*; complex vs real."""
import numpy as np, scipy.linalg as la, os
HERE = os.path.dirname(os.path.abspath(__file__))
d = np.load(os.path.join(HERE, 'e6_data.npz')); ad = d['ad']; n = 78
a2 = np.load(os.path.join(HERE, 'a2_data.npz')); inv = a2['inv']; deg = list(a2['degrees'])
X = {int(m): inv[i] for i, m in enumerate(deg)}
def opmat(v): return sum(v[i] * ad[i] for i in range(n) if abs(v[i]) > 0)
def kernel(M, rel=1e-8):
    u, s, vh = np.linalg.svd(M); r = int(np.sum(s > rel * s[0])); return vh[r:].conj().T
def nullity(M, rel=1e-8):
    s = np.linalg.svd(M, compute_uv=False); return int(np.sum(s < rel * s[0]))
def rank_drops(A, B, generic_rank, seeds=(11, 22), bound=1e5):
    """finite t with rank(A + t B) < generic_rank, via compressed generalised eigenvalues, common to two compressions."""
    res = []
    for sd in seeds:
        rg = np.random.default_rng(sd)
        R = rg.normal(size=(generic_rank, A.shape[0])) + 1j * rg.normal(size=(generic_rank, A.shape[0]))
        S = rg.normal(size=(A.shape[1], generic_rank)) + 1j * rg.normal(size=(A.shape[1], generic_rank))
        w = la.eigvals(R @ A @ S, -R @ B @ S)
        res.append([t for t in w if np.isfinite(t) and abs(t) < bound])
    if not res[0] or not res[1]: return []
    out = []
    for t0 in res[0]:
        if min(abs(t0 - u) for u in res[1]) < 1e-5 * max(1, abs(t0)) and all(abs(t0 - c) > 1e-6 * max(1, abs(t0)) for c in out):
            out.append(t0)
    return out
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
A8, A16 = opmat(X[8]), opmat(X[16])
gen_null = nullity(A8 + 0.3137 * A16); print(f"stage 1: generic dim Cent(x8 + t x16) = {gen_null}  (t generic)")
ts = rank_drops(A8, A16, n - gen_null)
print("stage 1 jump points t* (dim Cent, derived, centre):")
stage1 = []
for t0 in sorted(ts, key=lambda z: abs(z)):
    K = kernel(A8 + t0 * A16); k, der, z = structure(K); stage1.append((t0, K)); print(f"   t* = {t0:.6f}  ->  {k}, derived {der}, centre {z}")
print("   (record B874: the tuned x8 + t* x16 sits at a cubic root t*; dim 26 attained over the algebraic closure, never over R)")
# stage 2
A14 = opmat(X[14])
for t0, K in stage1:
    B14, B16 = A14 @ K, A16 @ K                       # restricted to Cent(x1)
    gen2 = nullity(B14 + 0.2718 * B16); print(f"\nstage 2 at t* = {t0:.6f}: generic joint dim Cent(x1, x14 + u x16) = {gen2}")
    us = rank_drops(B14, B16, K.shape[1] - gen2)
    for u0 in sorted(us, key=lambda z: abs(z)):
        J = K @ kernel(B14 + u0 * B16); k, der, z = structure(J)
        print(f"   u* = {u0:.6f} (Im {u0.imag:+.3e}) -> joint centralizer dim {k}, derived {der}, centre {z}")
print("\n(record B892/B950: a joint centralizer of dim 14 = su(3)+su(2)+u(1)^3 (derived 11, centre 3) at a COMPLEX wall point; the 12-dim SM algebra is NOT attained)")
