#!/usr/bin/env python3
"""T07 -- THE LOCKING DICTIONARY.  Statement: a world's discrete choices (h: weak handedness, m: matter label, t: arrow)
mod CPT form a rank-2 space; a nonzero measured sign is CPT-even and of one of four types by (eP, eC, eT); the object's
closing lattice (T4, T6 mirror-odd; T7 = T3 flow-odd; no (odd,odd) axis) admits six dictionaries to {P, C, T}, of which
three pass the axes' own semantics and one (c = P, gamma5 = T) is populated by measured signs.
THE LIMIT (added 2026-09-02, the seat testing its own proposal): the product of two signs of the same type is EVEN
type (eP eC eT = (+,+,+)), and EVEN quantities are exactly what no closing bit reaches.  So the dictionary fixes no
relative sign between two handedness signs, two CP/T asymmetries, or two EDMs: its only content is CPT's own type
bookkeeping.  The four "predictions" of the first version are withdrawn; no falsifiable sign statement survives.
Comparison class: all dictionaries."""
import itertools
from _common import say, write
worlds = list(itertools.product((0, 1), repeat=3))
CPT = lambda w: tuple((c + 1) % 2 for c in w)
classes = {min(w, CPT(w)) for w in worlds}
say(f"8 worlds, {len(classes)} physical classes mod CPT (rank 2)")
TYPES = {'EVEN': (1, 1, 1), 'W': (-1, -1, 1), 'K': (1, -1, -1), 'E': (-1, 1, -1)}
idx = {'P': 0, 'C': 1, 'T': 2}
axes = {'T4': (-1, 1), 'T6': (-1, 1), 'T7': (1, -1), 'T3': (1, -1), 'absent(odd,odd)': (-1, -1)}
adm = {}
for cf, gf in itertools.permutations('PCT', 2):
    table = {ax: [t for t, e in TYPES.items() if e[idx[cf]] == a and e[idx[gf]] == b][0] for ax, (a, b) in axes.items()}
    ok = TYPES[table['T4']][0] == -1 and TYPES[table['T7']][2] == -1
    say(f"  c={cf}, gamma5={gf}: {table} {'ADMISSIBLE' if ok else 'rejected'}")
    if ok: adm[f"c={cf},gamma5={gf}"] = table
say(f"admissible: {list(adm)}; testable now (rows populated by measured signs): c=P,gamma5=T")
# the limit: the ratio of two signs of one type is EVEN, and EVEN is outside every dictionary's image
mul = lambda a, b: tuple(x * y for x, y in zip(a, b))
ratios = {t: mul(e, e) for t, e in TYPES.items()}
assert all(r == TYPES['EVEN'] for r in ratios.values())
image = {TYPES[t] for table in adm.values() for t in table.values()}
assert TYPES['EVEN'] not in image
say("the product of two signs of the same type is EVEN for every type; EVEN is in no dictionary's image:")
say("  so no relative sign between two handedness signs, two CP/T asymmetries or two EDMs is fixed by any locking.")
say("  The dictionary's whole content is CPT's type bookkeeping.  No sign prediction survives (first-version P1-P4 withdrawn).")
write('T07', dict(classes=len(classes), admissible=adm, testable='c=P,gamma5=T', same_type_ratio_is_even=True,
                  even_in_image=False, predictions=[]))
