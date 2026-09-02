#!/usr/bin/env python3
"""T09 -- THE OBJECT'S POINT ON THE RULE'S CHARACTER VARIETY.  Statement: the rule acts on Fricke coordinates by
F: (x,y,z) -> (z, x, xz - y) (T03/T08); the fixed curve of F^2 (the monodromy of the object) is y = z = x/(x-1);
imposing the puncture condition of the once-punctured torus (Markoff: kappa = x^2 + y^2 + z^2 - xyz = 0, i.e.
tr[A,B] = -2) gives x^2 - 3x + 3 = 0, discriminant -3.  So the field Q(sqrt-3), the two conjugate points (the mirror
pair), and the fibre character of the geometric holonomy of m004 follow from the RULE and the PUNCTURE alone, before
any tetrahedron is drawn.  Verified: SL(2,C) matrices with these traces admit a parabolic T conjugating (A,B) to
sigma^2(A,B) = (ABA, AB), with T commuting with the puncture [A,B]: a representation of the bundle group, complete at
the cusp.  Comparison class: the metallic rules sigma_m: a -> a^m b, b -> a, m = 1..4 (fixed Markoff points exist for each; only m = 1 is\nidentified exactly here)."""
import sympy as sp, numpy as np, itertools
from _common import say, write
x, y, z = sp.symbols('x y z')
F = lambda p: (p[2], p[0], p[0]*p[2] - p[1])
kappa = lambda p: sp.expand(p[0]**2 + p[1]**2 + p[2]**2 - p[0]*p[1]*p[2])
say(f"F preserves kappa: {sp.simplify(kappa(F((x, y, z))) - kappa((x, y, z))) == 0}")
F2 = F(F((x, y, z)))
sol = sp.solve([sp.Eq(F2[0], x), sp.Eq(F2[1], y), sp.Eq(F2[2], z)], [y, z], dict=True)
say(f"fixed curve of F^2 (the monodromy sigma^2): {sol}")
curve = sol[0]
poly = sp.factor(sp.numer(sp.together(kappa((x, y, z)).subs(curve))))
say(f"kappa = 0 (Markoff / parabolic puncture) on the fixed curve: {poly} = 0")
minpoly = [f for f, _ in sp.factor_list(poly)[1] if sp.degree(f, x) >= 2][0]
roots = sp.solve(minpoly, x)
say(f"minimal polynomial {minpoly}, discriminant {sp.discriminant(minpoly, x)}, roots {roots}")
X = roots[0]; Y = sp.simplify(curve[y].subs(x, X)); Z = sp.simplify(curve[z].subs(x, X))
say(f"the object's fibre character: (x, y, z) = ({X}, {Y}, {Z}); y = z = conj(x): {sp.simplify(Y - sp.conjugate(X)) == 0}")
# numerical verification: build A, B with these traces, apply sigma^2, find the conjugator T
xv, yv, zv = complex(X), complex(Y), complex(Z)
A = np.array([[xv, -1], [1, 0]], dtype=complex)                       # tr A = x, det 1
# B = [[p, q],[r, s]] with p+s = y, ps - qr = 1, tr(AB) = z: choose q = 1, solve
p_, s_ = sp.symbols('p s')
# general solve: B = [[p, 1],[r, s]]: ps - r = 1; tr(AB) = x p - r + s... compute symbolically
r_ = sp.symbols('r')
Bsym = sp.Matrix([[p_, 1], [r_, s_]]); Asym = sp.Matrix([[xv, -1], [1, 0]])
eqs = [sp.Eq(p_ + s_, yv), sp.Eq(p_*s_ - r_, 1), sp.Eq((Asym*Bsym).trace(), zv)]
solB = sp.solve(eqs, [p_, r_, s_], dict=True)[0]
B = np.array([[complex(solB[p_]), 1], [complex(solB[r_]), complex(solB[s_])]])
comm = A @ B @ np.linalg.inv(A) @ np.linalg.inv(B)
say(f"tr A = {np.trace(A):.6f}, tr B = {np.trace(B):.6f}, tr AB = {np.trace(A@B):.6f}, tr[A,B] = {np.trace(comm):.6f} (parabolic puncture: -2)")
A2, B2 = A @ B @ A, A @ B                                             # sigma^2: a -> aba, b -> ab
# find T with T A T^-1 = A2, T B T^-1 = B2: linear system T A - A2 T = 0, T B - B2 T = 0
def lin(M, M2):
    rows = []
    for i in range(2):
        for j in range(2):
            row = np.zeros(4, dtype=complex)
            for k in range(2):
                row[i*2 + k] += M[k, j]; row[k*2 + j] -= M2[i, k]
            rows.append(row)
    return np.array(rows)
S = np.vstack([lin(A, A2), lin(B, B2)])
u, sv, vh = np.linalg.svd(S)
T = vh[-1].conj().reshape(2, 2); T = T / np.sqrt(np.linalg.det(T))
say(f"conjugator T: smallest singular value {sv[-1]:.2e}; ||T A T^-1 - ABA|| = {np.linalg.norm(T@A@np.linalg.inv(T) - A2):.2e}, ||T B T^-1 - AB|| = {np.linalg.norm(T@B@np.linalg.inv(T) - B2):.2e}")
# sigma^2 sends the puncture [A,B] to [ABA, AB] = (ABA)[A,B](ABA)^-1, so the peripheral section is T' = (ABA)^-1 T,
# which must commute with [A,B] and be parabolic (complete cusp)
Tp = np.linalg.inv(A @ B @ A) @ T
say(f"tr T = {np.trace(T):.6f}; corrected section T' = (ABA)^-1 T: tr T' = {np.trace(Tp):.6f} (parabolic: +-2), "
    f"||[T', [A,B]]|| = {np.linalg.norm(Tp@comm - comm@Tp):.2e} (abelian cusp group)")
conj_res = float(max(np.linalg.norm(T@A@np.linalg.inv(T) - A2), np.linalg.norm(T@B@np.linalg.inv(T) - B2)))
# comparison class: metallic rules sigma_m: a -> a^m b, b -> a.  With S_n the Chebyshev polynomials of the second kind
# in x = tr A, A^m = S_{m-1} A - S_{m-2} I, so the Fricke action is F_m(x,y,z) = (S_{m-1} z - S_{m-2} y, x, S_m z - S_{m-1} y).
# Fixed points of F_m^2 on the Markoff surface, up to the SL(2) lift signs (x,y,z) -> (e1 x, e2 y, e1 e2 z), are found
# numerically (Newton from random starts; symbolic elimination is unreliable at these degrees) and their minimal
# polynomials over Q recovered with PARI algdep.  m = 1 must reproduce x^2 - 3x + 3.
from snappy import pari
pari.set_real_precision(40)
def cheb_num(n, xx):
    S = [1.0 + 0j, xx]
    for k in range(2, n + 2): S.append(xx * S[-1] - S[-2])
    return S
def Fnum(p, m):
    xx, yy, zz = p; S = cheb_num(m + 1, xx); Sm2 = S[m - 2] if m >= 2 else 0
    return (S[m - 1] * zz - Sm2 * yy, xx, S[m] * zz - S[m - 1] * yy)
def Gfun(v, m, e):
    p = tuple(v); q = Fnum(Fnum(p, m), m)
    return np.array([q[0] - e[0]*p[0], q[1] - e[1]*p[1], q[2] - e[2]*p[2], p[0]**2 + p[1]**2 + p[2]**2 - p[0]*p[1]*p[2]])
def minpoly_of(val, maxdeg=8):
    """integer relation among 1, w, ..., w^k (PARI lindep) at double precision; accepted if small coefficients
    and residual < 1e-8; the algdep entry point rejects complex input in this cypari build."""
    w = pari(f'{val.real!r}+{val.imag!r}*I')
    for k in range(1, maxdeg + 1):
        v = pari.lindep([w ** i for i in range(k + 1)])
        if len(v) < k + 1: continue
        coeffs = [int(v[i]) for i in range(k + 1)]
        if max(abs(c) for c in coeffs) > 10**4 or coeffs[k] == 0: continue
        val_at = sum(c * complex(val) ** i for i, c in enumerate(coeffs))
        if abs(val_at) < 1e-8 * max(1.0, abs(val) ** k):
            if coeffs[k] < 0: coeffs = [-c for c in coeffs]
            return ' + '.join(f'{c}*x^{i}' if i else f'{c}' for i, c in enumerate(coeffs) if c) .replace('+ -', '- ')
    return None
rng = np.random.default_rng(0)
fields = {}
for mm in range(1, 5):
    found = {}
    for e in [(1, 1, 1), (-1, -1, 1), (-1, 1, -1), (1, -1, -1)]:
        for t in range(250):
            v = 2 * (rng.normal(size=3) + 1j * rng.normal(size=3))
            for it in range(80):
                g = Gfun(v, mm, e)
                if np.linalg.norm(g) < 1e-13: break
                J = np.zeros((4, 3), dtype=complex); h = 1e-7
                for k in range(3):
                    dv = np.zeros(3, dtype=complex); dv[k] = h; J[:, k] = (Gfun(v + dv, mm, e) - g) / h
                v = v - np.linalg.lstsq(J, g, rcond=None)[0]
            if np.linalg.norm(Gfun(v, mm, e)) < 1e-11 and min(abs(v)) > 1e-6:
                # polish to high precision by a few more exact-Jacobian steps is unnecessary for algdep at 1e-11? refine:
                for it in range(20):
                    g = Gfun(v, mm, e); J = np.zeros((4, 3), dtype=complex); h = 1e-9
                    for k in range(3):
                        dv = np.zeros(3, dtype=complex); dv[k] = h; J[:, k] = (Gfun(v + dv, mm, e) - g) / h
                    v = v - np.linalg.lstsq(J, g, rcond=None)[0]
                mp = minpoly_of(complex(v[0])) or 'unrecognised'
                found.setdefault(mp, set()).add(str(e))
    fields[mm] = {mp: sorted(sg) for mp, sg in found.items()}
    say(f"  metallic m={mm}: fixed Markoff points exist: {bool(found)}; minimal polynomials of x recognised at double precision (by sign twist): {fields[mm]}")
say("  (m = 1 is exact above; for m >= 2 the points exist but their fields are not identified at this precision:")
say("   a comparison-class statement only -- the metallic fibre fields need exact arithmetic to name.)")
write('T09', dict(fixed_curve={str(k): str(v) for k, v in curve.items()}, minpoly=str(minpoly), discriminant=int(sp.discriminant(minpoly, x)),
                  fibre_character=[str(X), str(Y), str(Z)], tr_comm=[np.trace(comm).real, np.trace(comm).imag],
                  conj_residual=conj_res,
                  trT=[np.trace(Tp).real, np.trace(Tp).imag], cusp_abelian_residual=float(np.linalg.norm(Tp@comm - comm@Tp)),
                  metallic_fixed_points={str(k): v for k, v in fields.items()}))
