#!/usr/bin/env python3
"""T17 -- THE OBJECT'S ONE beta-ODD OUTPUT IS 2-TORSION, AND UNDER THE DECLARED DICTIONARY IT IS THE STRONG CP BIT.
Statement.  (i) The Chern-Simons invariant of m004 is odd under every orientation-reversing self-isometry; the closing
lattice's two axes c (mirror) and gamma5 (flow reversal) are both orientation-reversing self-isometries, so CS is
(c-odd, gamma5-odd) -- the lattice's ABSENT axis -- and hence 2-torsion: CS in {0, 1/4} mod 1/2 (B1224/B1227), m004 at 0.
(ii) Under the dictionary c = P, gamma5 = T (T07, declared input A8) the (c-odd, gamma5-odd) type is (P-odd, C-even,
T-odd) = the E type; the Standard Model's only dimensionless E-type parameter is the strong phase theta-bar, whose
CP-conserving values are the 2-torsion of R/2piZ, {0, pi}.  The unique group isomorphism of the two Z/2's maps m004's
datum to theta-bar = 0.  The K-type phases (delta_CKM, delta_PMNS: C-odd, T-odd) map to the gamma5 row, a free bit.
So the dictionary predicts: strong CP conserved, weak CP free.  (iii) Failability (MB12): a chiral manifold has generic
CS (593 of 594 chiral census manifolds are not 2-torsion), so the prediction is not vacuous; the class-sibling m003
sits at 1/4 -> theta-bar = pi; the closed worlds m004(1,n) have definite-sign CS with |CS| -> 0 only as n -> infinity.
Not claimed: the value of theta-bar beyond the bit; that the object's mirror IS the strong sector's CP (declared, A8);
B813's refutation of the VALUE identification CS = theta stands (this is a Z/2 -> Z/2 map, no coefficient slot).
Comparison class: the 14-member field class of m004, the one-cusped census, m004's (1,n) closings."""
import snappy, warnings
from _common import say, write
warnings.filterwarnings('ignore')
M = snappy.Manifold('m004'); G = M.symmetry_group()
# (i) the two closing axes are orientation-reversing self-isometries; CS is odd under both
isos = G.isometries(); rows = []
for iso in isos:
    A = iso.cusp_maps()[0]; det = int(round(A.det())); msign = int(round(A[0, 0]))   # action on the meridian = H1 generator
    rows.append((det, msign))
classes = sorted(set(rows))
say(f"m004: |Sym| = {G.order()}, amphichiral = {G.is_amphicheiral()}; (cusp det, H1 sign) classes = {classes}")
assert classes == [(-1, -1), (-1, 1), (1, -1), (1, 1)]
say("  mirror c: det -1, H1 +1; flow reversal gamma5: det -1, H1 -1; both orientation-reversing, product preserving")
say("  CS(-M) = -CS(M), so CS is c-odd and gamma5-odd: the absent (odd,odd) axis -> 2 CS = 0 mod 1/2")
cs0 = float(M.chern_simons()) % 0.5
say(f"  CS(m004) mod 1/2 = {cs0:.12f}")
# the field class (T05) and the census as comparison classes
cls = ['m003','m004','m202','m203','m206','m207','m208','m410','m412','s118','s119','s594','s595','s596']
table = {}
for n in cls:
    N = snappy.Manifold(n); g = N.symmetry_group(); c = float(N.chern_simons()) % 0.5
    table[n] = dict(amph=bool(g.is_amphicheiral()), cs=round(c, 9))
say("  field class: " + ", ".join(f"{n}:{'A' if v['amph'] else 'c'}:{v['cs']}" for n, v in table.items()))
amph_ok = all((min(v['cs'], 0.5 - v['cs']) < 1e-6 or abs(v['cs'] - 0.25) < 1e-6) for v in table.values() if v['amph'])
assert amph_ok and abs(table['m003']['cs'] - 0.25) < 1e-6 and table['m004']['cs'] < 1e-6
tot = tor = 0
for N in snappy.OrientableCuspedCensus(num_cusps=1)[:600]:
    try:
        if N.symmetry_group().is_amphicheiral(): continue
        c = float(N.chern_simons()) % 0.5
    except Exception: continue
    tot += 1; tor += (min(c, 0.5 - c) < 1e-6 or abs(c - 0.25) < 1e-6)
say(f"  chiral one-cusped census (first 600): {tot}, at a 2-torsion value: {tor} (the bite: a chiral object predicts generic CP violation)")
# (ii) the type bookkeeping under c = P, gamma5 = T
TYPES = {'W': (-1, -1, 1), 'K': (1, -1, -1), 'E': (-1, 1, -1)}   # (eP, eC, eT)
obj = {'T4 (chirality row)': (-1, 1), 'T7 = T3 (time row)': (1, -1), 'CS (absent axis)': (-1, -1)}   # (c, gamma5) parities
typ = {ax: [t for t, e in TYPES.items() if (e[0], e[2]) == (ec, eg)][0] for ax, (ec, eg) in obj.items()}
say(f"  under c = P, gamma5 = T: {typ}")
assert typ['CS (absent axis)'] == 'E'
say("  SM dimensionless beta-odd parameters by type: theta-bar_QCD = E (P-odd, C-even, T-odd); delta_CKM, delta_PMNS = K")
say("  prediction of the dictionary: the E-type parameter sits at 2-torsion -> theta-bar in {0, pi}; the K-type phases are a free bit")
say("  m004's datum 0 -> theta-bar = 0 (the identity of Z/2); m003's 1/4 -> theta-bar = pi")
# (iii) the closed worlds
M.chern_simons(); fill = {}
for n in (2, 5, 10, 40):
    M.dehn_fill((1, n)); a = float(M.chern_simons()); M.dehn_fill((1, -n)); b = float(M.chern_simons()); fill[n] = (round(a, 6), round(b, 6))
    assert a < 0 < b and abs(a + b) < 1e-9
M.dehn_fill((0, 0))
say(f"  closings m004(1,+-n): CS = {fill} -- definite sign, mirror pairs opposite, |CS| -> 0 as n -> infinity (B303)")
write('T17', dict(sym_order=G.order(), parity_classes=classes, cs_m004=cs0, field_class=table, chiral_census=[tot, tor],
                  types=typ, fillings=fill, prediction='theta-bar in {0,pi}, m004 -> 0; weak phases free', dictionary='c=P,gamma5=T (A8)'))
