# oa-distill — the golden mapping torus, distilled to what is proved

This repository holds only what survived independent recomputation in the origin-axiom physics seat (September 2026):
eight theorem modules, each a statement with its **comparison class**, a script that proves it in about a minute, and a
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
| T08 | The rule's chain (on-site ±V) has gap labels in ℤ + ℤ/φ; its transfer-matrix trace map is the rule's Fricke action with invariant 4 + 4V² conserved; the K4 rules agree, the 4-letter rule does not. The rule's measured physics; the scale is the experimenter's (A7). | `theorems/T08_rule_spectrum.py` |
| T09 | The fixed curve of the monodromy's Fricke action is y = z = x/(x−1); with the puncture (Markoff κ = 0) it gives x² − 3x + 3 = 0: ℚ(√−3), the mirror pair, and the object's fibre character follow from the rule and the puncture alone; verified as a complete cusped representation of the bundle group. Metallic rules give degree ≥ 3 fields. | `theorems/T09_object_point.py` |
| T10 | The chain sits on κ = 4(1 + V²) ≥ 4, the object on κ = 0: no potential reaches the object; no hyperbolic filling (\|p\|,\|q\| ≤ 8) gives the object a real fibre-boundary trace: the two never share a level set. | `theorems/T10_level_sets.py` |
| T08 | The Fibonacci chain built from the rule has gaps labelled by ℤ + ℤ/φ (gap-labelling theorem; measured in photonic, polaritonic and cold-atom chains); the spectrum is the bounded-orbit set of the rule's own Fricke action, with the Fricke invariant conserved; all four K4 rules give the same labels, the four-letter rule does not. This is the rule's measured physics; its scale is the experimenter's (A7). | `theorems/T08_rule_spectrum.py` |

Run everything: `./run_all.sh` (SnapPy ≥ 3.3, sympy ≥ 1.12). Tests: `pytest -q tests`.

Provenance: origin-axiom, branch `claude/physics-seat-evaluation-8dkbrl`, cells R53–R56 and Phases F/G. Nothing here is
banked to origin-axiom `main` by this repository.
