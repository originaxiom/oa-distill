"""A7 -- the hypercharge question at the su(3)+su(3) landing, settled at the weight level.

The record (B1098/B1100/B1102) composes the object's SL(2,C) holonomy with the principal sl(2) of one
su(3) factor of su(3)^3 in e6.  The centralizer is the other two factors, su(3)_1 + su(3)_2 (rank 4).
The 27 of e6 restricts to su(3)^3 as (3,3bar,1) + (3bar,1,3) + (1,3,3bar); the eaten factor turns
each of its 3's into a spin-1 triplet, so under su(3)_1 + su(3)_2 the 27 is

        (3, 3bar)  +  3 x (3bar, 1)  +  3 x (1, 3).

This is all that is needed to test B1102's two claims by hand:

  side 1:  there are exactly 18 rational directions Y in the rank-4 Cartan whose 27 eigenvalues are the
           Standard Model hypercharge multiset  6Y = {1 x6, 2 x6, -3 x4, -4 x3, -2 x3, 0 x2, 3 x2, 6 x1};
  side 2:  none of them lies in the Cartan of a single factor (none commutes with a whole su(3)).

Parametrize Y by its values a = (a1,a2,a3) on the weights of the 3 of su(3)_1 (sum 0) and
b = (b1,b2,b3) on the weights of the 3 of su(3)_2 (sum 0).  The 27 eigenvalues are
    {a_i - b_j}  (9)   u   {-a_i} x3 each  (9)   u   {b_j} x3 each  (9).
Side 2 is a one-line count: if a = 0 or b = 0 the multiset has at least nine zeros; the target has two.
Side 1 is a finite search: the six values -a_i, b_j each occur with multiplicity >= 3, so they are drawn
from the target's multiplicity->=3 values {1, 2, -3, -4, -2}; enumerate.
"""
from fractions import Fraction
from itertools import product
from collections import Counter

TARGET = Counter({1: 6, 2: 6, -3: 4, -4: 3, -2: 3, 0: 2, 3: 2, 6: 1})
assert sum(TARGET.values()) == 27 and sum(k * v for k, v in TARGET.items()) == 0

def spectrum(a, b):
    c = Counter()
    for ai in a:
        for bj in b:
            c[ai - bj] += 1
    for ai in a:
        c[-ai] += 3
    for bj in b:
        c[bj] += 3
    return c

# side 2: a pure direction has >= 9 zeros
def zeros_if_pure():
    worst = 27
    vals = [-4, -3, -2, 0, 1, 2, 3, 6]
    for b in product(vals, repeat=3):
        if sum(b) != 0:
            continue
        worst = min(worst, spectrum((0, 0, 0), b)[0])
        worst = min(worst, spectrum(b, (0, 0, 0))[0])
    return worst

# side 1: exhaustive search over the multiplicity->=3 values (the record's own completeness argument)
BIG = [1, 2, -3, -4, -2]
solutions = []
for na in product(BIG, repeat=3):          # na = -a
    if sum(na) != 0:
        continue
    a = tuple(-x for x in na)
    for b in product(BIG, repeat=3):
        if sum(b) != 0:
            continue
        if spectrum(a, b) == TARGET:
            solutions.append((a, b))

# broader control: allow every target value for a and b, confirm nothing new appears
ALL = list(TARGET)
ctrl = 0
for na in product(ALL, repeat=3):
    if sum(na) != 0:
        continue
    a = tuple(-x for x in na)
    for b in product(ALL, repeat=3):
        if sum(b) != 0 or spectrum(a, b) != TARGET:
            continue
        ctrl += 1

pure = [s for s in solutions if all(x == 0 for x in s[0]) or all(x == 0 for x in s[1])]

# orbits under S3 x S3 (Weyl group of su(3)_1 + su(3)_2 acting on the two triples)
def canon(s):
    return (tuple(sorted(s[0])), tuple(sorted(s[1])))
orbits = sorted(set(canon(s) for s in solutions))

# the mirror: swap the two factors and conjugate (3 <-> 3bar): (a, b) -> (-b, -a)
def mirror(s):
    a, b = s
    return (tuple(-x for x in b), tuple(-x for x in a))
mirror_closed = all(canon(mirror(s)) in set(map(canon, solutions)) for s in solutions)
swap_closed = all(canon((s[1], s[0])) in set(map(canon, solutions)) for s in solutions)

print("minimal number of zero eigenvalues for a direction pure on one factor:", zeros_if_pure(), "(target has 2)")
print("number of rational directions carrying the SM hypercharge multiset:", len(solutions), "(control over all target values:", ctrl, ")")
print("directions pure on one factor among them:", len(pure))
print("S3xS3 orbits:", len(orbits), orbits)
print("closed under the mirror (swap + conjugation):", mirror_closed, "; closed under the plain swap:", swap_closed)
for s in solutions[:3]:
    print("  example  a/6 =", [Fraction(x, 6) for x in s[0]], " b/6 =", [Fraction(x, 6) for x in s[1]])
assert len(solutions) == 18 and ctrl == 18 and len(pure) == 0 and len(orbits) == 2 and mirror_closed and not swap_closed
print("A7 PASS: B1102 side 1 (18 directions), side 2 (none color-commuting), B1118 (two S3xS3 orbits fused by the mirror, not by the swap) reproduced at the weight level")
