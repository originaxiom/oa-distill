"""A8 -- the record's non-abelian hatch (B1098), checked on the audit's own e6.

Claim (B1098): composing the object's Zariski-dense SL(2,C) holonomy with an sl(2) embedding into e6 leaves as
unbroken algebra the centralizer of that sl(2).  For the principal sl(2) of one su(3) factor of su(3)^3 in e6
(the "A2 class") the centralizer is su(3)+su(3): dimension 16, rank 4, semisimple.  For the sl(2) of a single
root (the minimal "A1 class") the centralizer has dimension 35 (su(6)), rank 5.

The density lemma itself is standard (an irreducible non-elementary subgroup of SL(2,C) is Zariski dense, so its
centralizer under any algebraic map equals the centralizer of the image of SL(2)); this script checks the two
centralizers it is applied to, on the e6 built in A1 (no data from the record).
"""
import numpy as np, os, itertools
HERE = os.path.dirname(os.path.abspath(__file__))
d = np.load(os.path.join(HERE, 'e6_data.npz')); ad = d['ad']; C = d['C']; roots = [tuple(int(x) for x in r) for r in d['roots']]
n = 78
ridx = {r: 6 + i for i, r in enumerate(roots)}
def ip(a, b): return int(np.array(a) @ C @ np.array(b))
def vec(i): v = np.zeros(n); v[i] = 1.0; return v
def coroot(r): return sum(r[i] * vec(i) for i in range(6))
def opmat(v): return np.tensordot(v, ad, axes=(0, 0))
def nullspace(M, tol=1e-9):
    u, s, vt = np.linalg.svd(M)
    return vt[(s > tol).sum():].T if M.shape[0] >= M.shape[1] else vt[np.sum(s > tol):].T
def centralizer(elts):
    M = np.vstack([opmat(x) for x in elts])
    u, s, vt = np.linalg.svd(M, full_matrices=True)
    rank = int((s > 1e-9).sum())
    return vt[rank:].T          # columns span the joint kernel

# --- three mutually orthogonal A2 subsystems (the trinification su(3)^3)
pos = [r for r in roots if sum(r) > 0]
def a2_from(r1):
    for r2 in pos:
        if ip(r1, r2) == -1:
            return [r1, r2]
    return None
def orth(r, S): return all(ip(r, s) == 0 for s in S)
A2s = []
used = []
for r in pos:
    if not orth(r, used): continue
    pair = None
    for r2 in pos:
        if ip(r, r2) == -1 and orth(r2, used):
            pair = [r, r2]; break
    if pair is None: continue
    A2s.append(pair); used += pair + [tuple(x + y for x, y in zip(*pair))]
    if len(A2s) == 3: break
assert len(A2s) == 3, A2s
for i, j in itertools.combinations(range(3), 2):
    assert all(ip(a, b) == 0 for a in A2s[i] for b in A2s[j])
print("three mutually orthogonal A2 subsystems found:", A2s)

def sl2_from_roots(rs, weights):
    """principal-type sl2: h = sum weights_i * coroot(r_i), e = sum e_{r_i}; f solved in span of e_{-r_i}."""
    h = sum(w * coroot(r) for w, r in zip(weights, rs))
    e = sum(vec(ridx[r]) for r in rs)
    negs = [vec(ridx[tuple(-x for x in r)]) for r in rs]
    # [e, f] = h with f = sum c_k neg_k  ->  linear solve
    A = np.stack([opmat(e) @ g for g in negs], axis=1)
    c, *_ = np.linalg.lstsq(A, h, rcond=None)
    f = sum(ck * g for ck, g in zip(c, negs))
    assert np.allclose(opmat(e) @ f, h, atol=1e-9)
    assert np.allclose(opmat(h) @ e, 2 * e, atol=1e-9) and np.allclose(opmat(h) @ f, -2 * f, atol=1e-9)
    return h, e, f

def analyse(name, h, e, f, expect_dim):
    Z = centralizer([h, e, f]); dim = Z.shape[1]
    # derived algebra inside Z
    brs = [opmat(Z[:, i]) @ Z[:, j] for i in range(dim) for j in range(i + 1, dim)]
    D = np.linalg.matrix_rank(np.array(brs), tol=1e-8) if brs else 0
    # rank = dimension of the centralizer, within Z, of a generic element of Z
    rng = np.random.default_rng(0)
    x = Z @ rng.standard_normal(dim)
    inner = np.linalg.matrix_rank(opmat(x) @ Z, tol=1e-8)
    rank = dim - inner
    print(f"{name}: centralizer dim {dim} (expected {expect_dim}), derived dim {D}, rank {rank}")
    return Z, dim, D, rank

# the A2 class: principal sl2 of the third su(3) (h = 2(coroot1 + coroot2), e = e1 + e2)
h, e, f = sl2_from_roots(A2s[2], [2, 2])
Z, dim, D, rank = analyse("A2 class (principal sl2 of one trinification su(3))", h, e, f, 16)
assert dim == 16 and D == 16 and rank == 4
# it is exactly the other two su(3)'s: their root vectors and coroots lie in Z
S = []
for k in (0, 1):
    a, b = A2s[k]; c = tuple(x + y for x, y in zip(a, b))
    for r in (a, b, c):
        S += [vec(ridx[r]), vec(ridx[tuple(-x for x in r)])]
    S += [coroot(a), coroot(b)]
S = np.stack(S, axis=1)
assert np.linalg.matrix_rank(np.hstack([Z, S]), tol=1e-8) == 16
print("  the centralizer is spanned by the other two su(3)'s (16 = 8 + 8): su(3)+su(3), rank 4  -- B1098 A2 row reproduced")

# the minimal A1 class: sl2 of a single root
h1, e1, f1 = sl2_from_roots([A2s[2][0]], [1])
Z1, dim1, D1, rank1 = analyse("A1 class (single-root sl2)", h1, e1, f1, 35)
assert dim1 == 35 and D1 == 35 and rank1 == 5
print("  su(6) (dim 35, rank 5)  -- B1098 A1 row reproduced")

# the principal sl2 of e6 itself (the chain's entrance): centralizer 0
from numpy.linalg import matrix_rank
h6 = sum(coroot(r) for r in pos)            # 2 rho^vee
e6 = sum(vec(ridx[r]) for r in [tuple(int(i == j) for j in range(6)) for i in range(6)])
negs = [vec(ridx[tuple(-int(i == j) for j in range(6))]) for i in range(6)]
A = np.stack([opmat(e6) @ g for g in negs], axis=1); c, *_ = np.linalg.lstsq(A, h6, rcond=None); f6 = sum(ck * g for ck, g in zip(c, negs))
assert np.allclose(opmat(e6) @ f6, h6, atol=1e-8)
Zp = centralizer([h6, e6, f6]); print("principal sl2 of e6: centralizer dim", Zp.shape[1], "(expected 0)")
assert Zp.shape[1] == 0
print("A8 PASS")
