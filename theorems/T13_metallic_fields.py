#!/usr/bin/env python3
"""T13 -- THE METALLIC FIBRE FIELDS, EXACTLY.  Statement: for the metallic rules sigma_m: a -> a^m b, b -> a, the
squared Fricke action F_m^2 (Chebyshev form, T09) has a fixed CURVE on the character variety; intersecting it with the
Markoff surface kappa = 0 (parabolic puncture), up to the SL(2) lift signs, gives the fibre character's x-coordinate as
a root of an exact polynomial over Q:
    m = 1: x^2 - 3x + 3 (disc -3)      m = 2: x^4 - 4x^2 + 8 (invariant field Q(x^2) = Q(i))
    m = 3: an octic                     m = 4: x^8 - 7x^6 + 13x^4 + 8x^2 - 32 (invariant field a quartic)
each of which is EXACTLY the minimal polynomial of the corresponding bundle's fibre generator trace computed from
SnapPy's holonomy (T13b).  From the rule side, a quadratic field with conductor 3 occurs for m = 1 alone.
Comparison class: m = 1..4, four lift-sign twists; other sign twists give real (Fuchsian) or elliptic points."""
import sympy as sp
from _common import say, write
x, y, z = sp.symbols('x y z')
def cheb(n, xx):
    S = [sp.Integer(1), xx]
    for k in range(2, n + 2): S.append(sp.expand(xx * S[-1] - S[-2]))
    return S
def Fm(p, m):
    S = cheb(m + 1, p[0]); Sm2 = S[m - 2] if m >= 2 else sp.Integer(0)
    return (sp.expand(S[m - 1] * p[2] - Sm2 * p[1]), p[0], sp.expand(S[m] * p[2] - S[m - 1] * p[1]))
SIGNS = [(1, 1, 1), (-1, -1, 1), (-1, 1, -1), (1, -1, -1)]
out = {}
for m in range(1, 5):
    out[m] = {}
    for e in SIGNS:
        q = Fm(Fm((x, y, z), m), m)
        eqs = [sp.expand(q[0] - e[0]*x), sp.expand(q[1] - e[1]*y), sp.expand(q[2] - e[2]*z)]
        zsol = sp.solve(eqs[1], z)[0]                                     # the middle component is linear in z
        e1 = sp.numer(sp.together(eqs[0].subs(z, zsol))); e3 = sp.numer(sp.together(eqs[2].subs(z, zsol)))
        curve = sp.gcd(sp.Poly(e1, x, y), sp.Poly(e3, x, y)).as_expr()   # the common component = the fixed curve
        assert sp.Poly(curve, x, y).total_degree() > 0, "fixed locus is not a curve"
        kap = sp.numer(sp.together((x**2 + y**2 + z**2 - x*y*z).subs(z, zsol)))
        R = sp.resultant(curve, kap, y)
        fac = [(int(sp.degree(f, x)), str(f)) for f, _ in sp.factor_list(sp.Poly(R, x).as_expr())[1] if sp.degree(f, x) >= 1]
        out[m][str(e)] = fac
        say(f"m={m} sign {e}: fixed curve degree {sp.Poly(curve, x, y).total_degree()}; Markoff intersection: {fac}")
geometric = {1: 'x**2 - 3*x + 3', 2: 'x**4 - 4*x**2 + 8', 3: 'x**8 - 4*x**7 + 4*x**6 - x**5 + 8*x**4 - 11*x**3 - 4*x**2 + 3*x + 6',
             4: 'x**8 - 7*x**6 + 13*x**4 + 8*x**2 - 32'}
found = {m: any(f == geometric[m] for fl in out[m].values() for _, f in fl) for m in geometric}
say(f"geometric fibre polynomial (matching T13b's SnapPy holonomy) present on the rule side: {found}")
quad_disc = {m: [int(sp.discriminant(sp.sympify(f), x)) for fl in out[m].values() for d, f in fl if d == 2 and sp.sympify(f).subs(x, sp.I).is_real is False] for m in out}
say(f"quadratic factors' discriminants per m: {quad_disc}  (conductor 3 <=> disc -3: m = 1 only)")
write('T13', dict(factors={str(m): v for m, v in out.items()}, geometric=geometric, geometric_found={str(m): v for m, v in found.items()},
                  quadratic_discriminants={str(m): v for m, v in quad_disc.items()}))
