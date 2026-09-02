#!/usr/bin/env python3
"""A2 -- the principal sl2 in e6, the binary tetrahedral group 2T inside SL(2), and Cent_{e6}(2T).
Method: decompose e6 under the principal sl2 into irreducibles V_{2k} (dims 3, 9, 11, 15, 17, 23), build the monomial
basis of each by f-strings from a highest-weight vector, and let g in SL(2) act by Sym^{2k}(g) exactly.  Checks: the
action is a group homomorphism and preserves the bracket.  Recomputes the record's B854: Cent(2T) = u(1)^4.
CHOICES: the sl2 class (principal, as the record); the SL(2) identification; nothing else."""
import numpy as np, os, itertools
from math import comb
HERE = os.path.dirname(os.path.abspath(__file__))
d = np.load(os.path.join(HERE, 'e6_data.npz')); ad = d['ad']; C = d['C']; roots = [tuple(r) for r in d['roots']]; n = 78
def opmat(v): return sum(v[i] * ad[i] for i in range(n) if v[i] != 0)
def br(u, v): return opmat(u) @ v
def vec(i): v = np.zeros(n); v[i] = 1; return v
c = np.linalg.solve(C.T.astype(float), 2 * np.ones(6)); h = np.zeros(n); h[:6] = c
e = sum(vec(6 + roots.index(tuple(int(i == j) for j in range(6)))) for i in range(6))
neg = [6 + roots.index(tuple(-int(i == j) for j in range(6))) for i in range(6)]
A = np.array([[br(e, vec(neg[i]))[k] for i in range(6)] for k in range(n)]); fi, *_ = np.linalg.lstsq(A, h, rcond=None)
f = sum(fi[i] * vec(neg[i]) for i in range(6))
adh, ade, adf = opmat(h), opmat(e), opmat(f)
assert np.allclose(ade @ adf - adf @ ade, adh) and np.allclose(adh @ ade - ade @ adh, 2 * ade)
# isotypic decomposition: highest-weight vectors = ker(ade) within each ad h eigenspace (adh is diagonal in this basis)
wts = np.round(np.diag(adh)).astype(int)
blocks = []          # (m=2k, unit basis vectors w_0..w_m, diagonal rescaling d_j to the monomial basis)
for wt in sorted(set(wts), reverse=True):
    if wt <= 0: continue
    idxs = np.where(wts == wt)[0]
    sub = ade[:, idxs]
    u, s, vh = np.linalg.svd(sub); null = vh[np.sum(s > 1e-9):].conj()
    for row in null:
        v = np.zeros(n, dtype=complex); v[idxs] = row; v /= np.linalg.norm(v)
        m = wt; ws = [v]; a = []
        for j in range(m):
            nxt = br(f, ws[-1]); nrm = np.linalg.norm(nxt); a.append(nrm); ws.append(nxt / nrm)   # f w_j = a_j w_{j+1}
        dcoef = [1.0]
        for j in range(m): dcoef.append(dcoef[-1] * (m - j) / a[j])                                  # w_j = d_j m_j
        blocks.append((m, ws, np.array(dcoef)))
print("irreps under the principal sl2 (dims):", sorted(len(ws) for _, ws, _ in blocks), " total", sum(len(ws) for _, ws, _ in blocks))
W = np.column_stack([w for _, ws, _ in blocks for w in ws]).astype(complex)
Winv = np.linalg.inv(W); print("basis change cond:", np.linalg.cond(W))
# consistency: e w_{j} should equal j * (d_{j-1}/d_j)^{-1}... in monomials e m_j = j m_{j-1}; w_j = d_j m_j =>
# e w_j = d_j j m_{j-1} = j d_j / d_{j-1} w_{j-1}
ok = True
for m, ws, dc in blocks:
    for j in range(1, m + 1):
        ok &= np.allclose(br(e, ws[j]), j * dc[j] / dc[j - 1] * ws[j - 1], atol=1e-6 * max(1, abs(j * dc[j] / dc[j - 1])))
print("e-action consistent with the monomial model in every string:", ok)
def sym_power(g, m):
    Mx = np.zeros((m + 1, m + 1), dtype=complex)
    for j in range(m + 1):
        poly = np.zeros(m + 1, dtype=complex); poly[0] = 1
        for _ in range(m - j): poly = np.convolve(poly, [g[0, 0], g[1, 0]])[:m + 1]
        for _ in range(j): poly = np.convolve(poly, [g[0, 1], g[1, 1]])[:m + 1]
        Mx[:, j] = poly
    return Mx
def Ad(g):
    B = np.zeros((n, n), dtype=complex); pos = 0
    for m, ws, dc in blocks:
        S = sym_power(g, m); Dm = np.diag(dc); B[pos:pos + m + 1, pos:pos + m + 1] = np.linalg.inv(Dm) @ S @ Dm; pos += m + 1
    return W @ B @ Winv
I2 = np.array([[1j, 0], [0, -1j]]); J2 = np.array([[0, 1], [-1, 0]]); K2 = np.array([[0, 1j], [1j, 0]])
g1 = I2; g2 = 0.5 * (np.eye(2) + I2 + J2 + K2)
# the 24 elements of 2T (closure), then the invariants block by block: P_m = (1/24) sum_g Sym^m(g) is the projector
# onto 2T-invariants in Sym^m; this needs no e6 numerics.  The e6 vector of an invariant monomial combination u is
# sum_j u_j m_j = sum_j (u_j / d_j) w_j, assembled in extended precision because d_j spans a wide range.
import mpmath as mp
mp.mp.dps = 40
els = [np.eye(2, dtype=complex)]; frontier = [np.eye(2, dtype=complex)]
def key(m): return tuple(np.round(m.flatten(), 8))
seen = {key(els[0])}
while frontier:
    new = []
    for m in frontier:
        for g in (g1, g2):
            pgm = m @ g
            if key(pgm) not in seen: seen.add(key(pgm)); els.append(pgm); new.append(pgm)
    frontier = new
print("|2T| =", len(els))
def sym_power_mp(g, m):
    Mx = mp.matrix(m + 1, m + 1)
    for j in range(m + 1):
        poly = [mp.mpc(0)] * (m + 1); poly[0] = mp.mpc(1)
        for _ in range(m - j):
            new = [mp.mpc(0)] * (m + 1)
            for i in range(m + 1):
                if poly[i] != 0:
                    new[i] += poly[i] * g[0, 0]
                    if i + 1 <= m: new[i + 1] += poly[i] * g[1, 0]
            poly = new
        for _ in range(j):
            new = [mp.mpc(0)] * (m + 1)
            for i in range(m + 1):
                if poly[i] != 0:
                    new[i] += poly[i] * g[0, 1]
                    if i + 1 <= m: new[i + 1] += poly[i] * g[1, 1]
            poly = new
        for i in range(m + 1): Mx[i, j] = poly[i]
    return Mx
invariants = []
for m, ws, dc in blocks:
    P = mp.matrix(m + 1, m + 1)
    for g in els:
        gm = mp.matrix([[mp.mpc(g[0, 0]), mp.mpc(g[0, 1])], [mp.mpc(g[1, 0]), mp.mpc(g[1, 1])]])
        P += sym_power_mp(gm, m)
    P = P / len(els)
    # rank of P = number of invariants; column space basis via SVD in float (P entries O(1))
    Pf = np.array([[complex(P[i, j]) for j in range(m + 1)] for i in range(m + 1)])
    u, sv, vh = np.linalg.svd(Pf); r = int(np.sum(sv > 1e-9))
    for k in range(r):
        uvec = u[:, k]                                    # invariant combination of monomials m_j
        coeffs = [mp.mpc(uvec[j]) / mp.mpf(dc[j]) for j in range(m + 1)]
        # normalise coefficients (so the e6 vector is O(1))
        mx = max(abs(cf) for cf in coeffs); coeffs = [cf / mx for cf in coeffs]
        x = sum(complex(coeffs[j]) * ws[j] for j in range(m + 1))
        x /= np.linalg.norm(x)
        invariants.append((m, x))
    print(f"  Sym^{m}: {r} invariant(s)")
print("dim Cent_e6(2T) =", len(invariants), " (record B854: 4); degrees:", [m for m, _ in invariants])
brs = max(np.linalg.norm(br(x, y)) for _, x in invariants for _, y in invariants)
print("abelian (max |[x_i,x_j]| = %.2e):" % brs, brs < 1e-8)
# invariance check with the (float) Ad of the generators, tolerance loose because Ad is float
A1, A2 = Ad(g1), Ad(g2)
print("Ad(g) x = x residuals:", [f"{np.linalg.norm(A2 @ x - x):.1e}" for _, x in invariants])
inv = np.array([x for _, x in invariants])
np.savez(os.path.join(HERE, 'a2_data.npz'), h=h, e=e, f=f, inv=inv, degrees=np.array([m for m, _ in invariants]), W=W)
