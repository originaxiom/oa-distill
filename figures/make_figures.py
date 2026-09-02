#!/usr/bin/env python3
"""Figures for PAPER_DRAFT: (1) spectrum of the Fibonacci chain vs V; (2) gap widths vs V (log-log) for the labels of T14."""
import numpy as np, json, os
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.linalg import eigh_tridiagonal
HERE = os.path.dirname(os.path.abspath(__file__))
w = 'a'
for _ in range(15): w = ''.join('ab' if c == 'a' else 'a' for c in w)      # 1597 sites for the picture
N = len(w); sgn = np.array([1.0 if c == 'a' else -1.0 for c in w])
Vs = np.linspace(0, 1.2, 121)
fig, ax = plt.subplots(figsize=(6, 7))
for V in Vs:
    ev = eigh_tridiagonal(V * sgn, np.ones(N - 1), eigvals_only=True)
    ax.plot(np.full(N, V), ev, ',', color='k', alpha=0.6)
ax.set_xlabel('modulation depth V (on-site ±V, unit hopping)'); ax.set_ylabel('energy E')
ax.set_title('Spectrum of the Fibonacci chain (rule a→ab, b→a), 1597 sites')
fig.tight_layout(); fig.savefig(os.path.join(HERE, 'spectrum_vs_V.png'), dpi=160); plt.close(fig)
T14 = json.load(open(os.path.join(HERE, '..', 'theorems', 'out', 'T14.json')))
fig, ax = plt.subplots(figsize=(6, 4.5))
for n, ws in T14['widths'].items():
    ax.loglog(T14['V'], ws, 'o-', ms=3, label=f'n={n}')
ax.set_xlabel('V'); ax.set_ylabel('gap width'); ax.set_title('Gap opening by label (T14): all slopes ≈ 1')
ax.legend(fontsize=7, ncol=2); fig.tight_layout(); fig.savefig(os.path.join(HERE, 'gap_widths_loglog.png'), dpi=160); plt.close(fig)
print('figures written')
