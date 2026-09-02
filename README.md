# oa-distill — the golden mapping torus, distilled to what is proved

This repository holds only what survived independent recomputation in the origin-axiom physics seat (September 2026):
seven theorem modules, each a statement with its **comparison class**, a script that proves it in about a minute, and a
test that fails if its numbers drift. There is no progress log, no verdict file, no arc. Choices are listed once, with
their price, in `AXIOMS.md`. What we do not claim is in `NOT_CLAIMED.md`. The discipline that keeps this tree honest is
in `RULES.md`.

| module | statement (one line) | proof |
|---|---|---|
| T01 | The rule σ: a→ab, b→a is one point of a free K4-torsor of equivalent rules; reversal is inner; the Fibonacci language is reversal-closed; the arrow is monoid non-surjectivity. | `theorems/T01_rule_torsor.py` |
| T02 | σ is orientation-reversing (det −1); its mapping torus is the Gieseking manifold; the object m004 is the orientation double cover, the mapping torus of σ²; σ acts on m004 as a mirror. | `theorems/T02_double_tick.py` |
| T03 | Sym(m004) = D4 with four orientation-reversing isometries; the mirror acts on characters as complex conjugation on ℚ(√−3), realised by explicit automorphism words; the bundle is invertible. | `theorems/T03_symmetry_and_mirror.py` |
| T04 | Every function of the object decomposes under D4; the closing lattice has rank 2 on the manifold (mirror, flow reversal) plus the lift sign; no object-canonical datum can orient m004. | `theorems/T04_closing_lattice.py` |
| T05 | The route field → conductor → SL(2,ℤ/N) → McKay emits only A₀, E6, E8; on hyperbolic manifolds only E6; 14 of the first 1200 census manifolds share m004's field, one commensurability class. | `theorems/T05_field_lookups.py` |
| T06 | Within that class, H1, cusp shape and subgroup growth single out m004; the sister's closings differ at every slope; the two triangulations define different deformation curves and different DGG theories. | `theorems/T06_object_beyond_field.py` |
| T07 | Reality's discrete choices mod CPT are rank 2; three dictionaries to the object's lattice survive; one is testable now and yields four sign predictions. | `theorems/T07_locking_dictionary.py` |

Run everything: `./run_all.sh` (SnapPy ≥ 3.3, sympy ≥ 1.12). Tests: `pytest -q tests`.

Provenance: origin-axiom, branch `claude/physics-seat-evaluation-8dkbrl`, cells R53–R56 and Phases F/G. Nothing here is
banked to origin-axiom `main` by this repository.
