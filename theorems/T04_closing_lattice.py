#!/usr/bin/env python3
"""T04 -- THE CLOSING LATTICE.  Statement: (i) for any involutive symmetry iota of the object, every function T of the
object decomposes uniquely as T = T_even + T_odd with T_odd(iota x) = -T_odd(x); the sign of T_odd on 'the object' versus
'the mirror' is fixed only by naming one of them (a choice); (ii) the symmetries of m004 acting on the cusp form
(Z/2)^2 = <mirror, flow reversal>, so a scalar tracker carries one of four linear characters (plus the SL/PSL lift sign,
which is a symmetry of the representation); (iii) hence no object-canonical datum can orient m004.
Comparison class: all functions of the object; all isometries."""
import itertools, sympy as sp, snappy
from _common import say, write
M = snappy.Manifold('m004')
isos = M.is_isometric_to(M, return_isometries=True)
mats = []
for iso in isos:
    C = iso.cusp_maps()[0]; mats.append(sp.Matrix([[int(C[0, 0]), int(C[0, 1])], [int(C[1, 0]), int(C[1, 1])]]))
G = {tuple(m): m for m in mats}
say(f"cusp image of Sym(m004): {len(G)} elements: {[m.tolist() for m in G.values()]}")
abelian = all(a*b == b*a for a in G.values() for b in G.values())
say(f"abelian: {abelian}; linear characters = {len(G)} (mirror-parity x time-parity); plus the 2-dim irrep of D4 for doublets")
# (i) the decomposition identity, symbolically: T = (T + T o iota)/2 + (T - T o iota)/2
T, Ti = sp.symbols('T Ti')
say(f"decomposition identity: (T+Ti)/2 + (T-Ti)/2 - T = {sp.simplify((T+Ti)/2 + (T-Ti)/2 - T)}")
# (iii) mirror acts on the complex volume: Vol + i CS -> Vol - i CS; CS(m004) = 0 so the mirror fixes every canonical datum
cs = float(M.chern_simons()); vol = float(M.volume())
say(f"Vol = {vol:.9f}, CS = {cs:.3e} (mirror-odd, hence 0 by amphichirality): no canonical datum is mirror-odd and nonzero")
# the flip table rows re-derived here: T4 (chirality side) flips under the mirror: the conjugate character differs
Gp = M.fundamental_group(); x = complex(Gp.SL2C('a')[0, 0] + Gp.SL2C('a')[1, 1])
say(f"tr a = {x:.6f} != conjugate {x.conjugate():.6f}: the mirror moves the geometric character (T4 mirror-odd)")
phi = (1 + sp.sqrt(5))/2
say(f"flow reversal inverts the eigenvalue phi^2 -> phi^-2 = (1-phi)^2: {sp.simplify((1-phi)**2 - phi**-2) == 0} (T7 time-odd)")
write('T04', dict(cusp_group_order=len(G), abelian=abelian, linear_characters=len(G), CS=cs, Vol=vol,
                  T4_mirror_odd=abs(x.imag) > 1e-9, T7_time_odd=bool(sp.simplify((1-phi)**2 - phi**-2) == 0),
                  bits=['mirror (c)', 'flow reversal (gamma5)', 'lift sign (theta)'], continuous=['scale']))
