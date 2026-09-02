#!/usr/bin/env python3
"""T15 -- THE METALLIC CHAINS.  Statement: the chain built from the metallic rule sigma_m: a -> a^m b, b -> a has gap
labels in Z + Z*omega_m with omega_m = 1/(lambda_m + 1) the frequency of the letter b, lambda_m = (m + sqrt(m^2+4))/2
(omega_1 = 1/phi^2, omega_2 = 1 - 1/sqrt2, omega_3 = (5 - sqrt13)/6); the golden omega fails for m = 2, 3.  Experiment class: silver- and bronze-mean
quasicrystal chains (realised in photonic lattices).  Comparison class: m = 1, 2, 3 and the wrong omega."""
import numpy as np
from scipy.linalg import eigh_tridiagonal
from _common import say, write
def fixed_point(m, iters):
    w = 'a'; rule = {'a': 'a' * m + 'b', 'b': 'a'}
    for _ in range(iters): w = ''.join(rule[c] for c in w)
    return w
def labels(ev, omega, n_gaps=10, ncand=8):
    ev = np.sort(ev); N = len(ev); gaps = np.diff(ev); order = np.argsort(gaps)[::-1][:n_gaps]
    res = []
    for i in sorted(order):
        ids = (i + 1) / N
        best = min(((abs(((n * omega) % 1) - ids), n) for n in range(-ncand, ncand + 1)), key=lambda t: t[0])
        res.append((round(ids, 6), best[1], round(best[0], 6)))
    return res
V = 0.5; out = {}
for m in (1, 2, 3):
    lam = (m + (m * m + 4) ** .5) / 2; omega = 1 / (lam + 1)
    iters = {1: 19, 2: 12, 3: 9}[m]
    w = fixed_point(m, iters); N = len(w)
    d = np.array([V if c == 'a' else -V for c in w]); ev = eigh_tridiagonal(d, np.ones(N - 1), eigvals_only=True)
    lab = labels(ev, omega); bad = labels(ev, ((1 + 5 ** .5) / 2) ** -1) if m != 1 else None
    freq_b = w.count('b') / N
    mx = max(r[2] for r in lab); mxbad = max(r[2] for r in bad) if bad else None
    # note: with |n| <= 8 only, a wrong omega cannot fit by equidistribution; that is what makes the test bite
    out[m] = dict(sites=N, omega=omega, freq_b=freq_b, max_residual=mx, max_residual_wrong_omega=mxbad, labels=lab)
    say(f"m={m}: {N} sites, omega_m={omega:.6f} (letter-b frequency {freq_b:.6f}); labels of the 10 largest gaps in Z + Z*omega_m: max residual {mx:.2e} (3/N = {3/N:.1e})"
        + (f"; with the golden omega instead: max residual {mxbad:.3f}" if bad else ""))
write('T15', {str(m): v for m, v in out.items()})
