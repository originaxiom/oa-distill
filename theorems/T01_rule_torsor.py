#!/usr/bin/env python3
"""T01 -- THE RULE.  Statement: the substitution sigma: a->ab, b->a is one point of a free K4-torsor of equivalent
rules (reversal-conjugation and swap-conjugation); its word-reversal is an INNER automorphism; the fixed-point language
is closed under reversal (no reading arrow); the only intrinsic arrow is that the rule is not surjective on the free
monoid (the word 'bb' has no preimage).  Comparison class: the four Fibonacci-type rules on {a,b}."""
import itertools, sympy as sp
from _common import say, write
R = {'a': 'ab', 'b': 'a'}
def apply(rule, w): return ''.join(rule[c] for c in w)
def inv(s): return ''.join(c.swapcase() for c in reversed(s))
def red(s):
    out = []
    for c in s:
        if out and out[-1] == c.swapcase(): out.pop()
        else: out.append(c)
    return ''.join(out)
# 1. the K4 orbit under reversal-conjugation (reverse every image) and swap-conjugation (a<->b)
def rev(rule): return {k: v[::-1] for k, v in rule.items()}
def swap(rule):
    s = lambda w: w.translate(str.maketrans('ab', 'ba'))
    return {s(k): s(v) for k, v in rule.items()}
orbit = {}
frontier = [('', R)]
while frontier:
    tag, r = frontier.pop()
    key = tuple(sorted(r.items()))
    if key in orbit: continue
    orbit[key] = tag or 'id'
    frontier += [(tag + 'R', rev(r)), (tag + 'S', swap(r))]
say(f"orbit of the rule under <reversal, swap>: size {len(orbit)}: {[dict(k) for k in orbit]}")
# 2. reversal is inner: sigma_rev = w sigma w^-1 for a word w
def conj(w, s): return red(w + s + inv(w))
found = None
for L in range(0, 5):
    for t in itertools.product('abAB', repeat=L):
        w = ''.join(t)
        if red(w) != w: continue
        if all(conj(w, R[x]) == rev(R)[x] for x in R): found = w; break
    if found is not None: break
say(f"reversal = conjugation by w = {found!r}")
# 3. abelianization and its determinant; sigma^2 = the figure-eight monodromy
M = sp.Matrix([[R[y].count(x) for y in 'ab'] for x in 'ab'])
say(f"incidence matrix {M.tolist()}, det {M.det()}, M^2 = {(M*M).tolist()}")
# 4. language closure under reversal (Sturmian), and the arrow: 'bb' has no preimage
w = 'a'
for _ in range(20): w = apply(R, w)
def factors(w, L): return {w[i:i+L] for i in range(len(w)-L+1)}
closed = all(factors(w, L) == {f[::-1] for f in factors(w, L)} for L in range(1, 13))
complexity = [len(factors(w, L)) for L in range(1, 9)]
say(f"fixed point length {len(w)}; language reversal-closed for lengths 1..12: {closed}; complexity p(n), n=1..8: {complexity}")
preimages = {}
for n in range(1, 9):
    for t in itertools.product('ab', repeat=n):
        img = apply(R, ''.join(t)); preimages.setdefault(img, []).append(''.join(t))
say(f"'bb' has a preimage among words of length <= 8: {'bb' in preimages};  'aaa' has: {'aaa' in preimages} ({preimages.get('aaa')})")
write('T01', dict(orbit_size=len(orbit), reversal_conjugator=found, det=int(M.det()), M2=[[int(v) for v in r] for r in (M*M).tolist()],
                  reversal_closed=closed, complexity=complexity, bb_has_preimage='bb' in preimages, aaa_preimage=preimages.get('aaa')))
