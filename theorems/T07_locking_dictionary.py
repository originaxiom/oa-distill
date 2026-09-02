#!/usr/bin/env python3
"""T07 -- THE LOCKING DICTIONARY.  Statement: a world's discrete choices (h: weak handedness, m: matter label, t: arrow)
mod CPT form a rank-2 space; a nonzero measured sign is CPT-even and of one of four types by (eP, eC, eT); the object's
closing lattice (T4, T6 mirror-odd; T7 = T3 flow-odd; no (odd,odd) axis) admits six dictionaries to {P, C, T}, of which
three pass the axes' own semantics and one (c = P, gamma5 = T) is populated by measured signs.  Predictions follow.
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
preds = ["P1: all K-type signs (CP/T rate asymmetries) are one bit; test: the sign of the leptonic CP phase",
         "P2: all W-type signs (handedness) are one bit up to the lift sign; needs a physical reading of the chord (none yet)",
         "P3: E-type signs (EDMs) are products of the two rows: relative signs of different EDMs are fixed; two measured EDMs decide",
         "P4: no K-type sign in the chirality row, no W-type sign in the time row (decidable now; not violated)"]
for p in preds: say("  " + p)
write('T07', dict(classes=len(classes), admissible=adm, testable='c=P,gamma5=T', predictions=preds))
