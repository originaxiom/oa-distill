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
| T07 | Reality's discrete choices mod CPT are rank 2; three dictionaries to the object's lattice survive; but any relative sign within one CPT type is EVEN and outside every dictionary's image, so no sign prediction exists (the first version's four withdrawn). | `theorems/T07_locking_dictionary.py` |
| T08 | The rule's chain (on-site ±V) has gap labels in ℤ + ℤ/φ; its transfer-matrix trace map is the rule's Fricke action with invariant 4 + 4V² conserved; the K4 rules agree, the 4-letter rule does not. The rule's measured physics; the scale is the experimenter's (A7). | `theorems/T08_rule_spectrum.py` |
| T09 | The fixed curve of the monodromy's Fricke action is y = z = x/(x−1); with the puncture (Markoff κ = 0) it gives x² − 3x + 3 = 0: ℚ(√−3), the mirror pair, and the object's fibre character follow from the rule and the puncture alone; verified as a complete cusped representation of the bundle group. Metallic rules give degree ≥ 3 fields. | `theorems/T09_object_point.py` |
| T10 | The chain sits on κ = 4(1 + V²) ≥ 4, the object on κ = 0: no potential reaches the object; no hyperbolic filling (\|p\|,\|q\| ≤ 8) gives the object a real fibre-boundary trace: the two never share a level set. | `theorems/T10_level_sets.py` |
| T13 | Exact: the fixed curve of the metallic monodromy's Fricke action meets the Markoff surface at x² − 3x + 3 (m = 1), x⁴ − 4x² + 8 (m = 2, invariant field ℚ(i)), an octic (m = 3), x⁸ − 7x⁶ + 13x⁴ + 8x² − 32 (m = 4); each equals the minimal polynomial of the bundle's fibre-generator trace from SnapPy's holonomy (T13b). From the rule side, conductor 3 is m = 1 alone. | `theorems/T13_metallic_fields.py`, `T13b_bundle_trace_fields.py` |
| T14 | Every gap of the rule's chain opens linearly in the modulation V (exponents 0.85–1.10 over ten labels), unlike the Harper model; the width ratios at fixed V are the scale-free observables and are tabulated. | `theorems/T14_gap_law.py` |
| T15 | The metallic chains' gap labels lie in ℤ + ℤ·ω_m with ω_m = 1/(λ_m + 1) the letter frequency (residual < 3/N for m = 1, 2, 3); the golden ω fails for m = 2, 3. | `theorems/T15_metallic_chains.py` |
| T17 | The object's Chern-Simons invariant is E-type (c-odd, γ₅-odd) and sits at 2-torsion; under the dictionary c=P, γ₅=T, the E-type SM parameter is θ̄_QCD; 2-torsion gives θ̄ ∈ {0,π}, m004 at 0 (strong CP conserved); weak CP phases are free bits. Bite: 1/594 chiral census manifolds sit at 2-torsion CS. | `theorems/T17_beta_odd_bit.py` |

Run everything: `./run_all.sh` (SnapPy ≥ 3.3, sympy ≥ 1.12). Tests: `pytest -q tests`.

Provenance: origin-axiom, branch `claude/physics-seat-evaluation-8dkbrl`, cells R53–R56 and Phases F/G. Nothing here is
banked to origin-axiom `main` by this repository.

## Audit (not theorems): the record's chain from E6 to "SM structure", rebuilt from scratch

`audit/AUDIT.md` with scripts `audit/A1`–`A8`: e6 built and verified; the principal sl(2) and the binary tetrahedral group;
Cent(2T) = u(1)⁴ (record reproduced); the centralizer ladder (reproduced); the complete stratification of the 2T-torus
by root hyperplanes (new): generic su(3)+4U(1), six hyperplanes of su(3)+su(2)+3U(1), SU(5)+2U(1) and SO(10)+U(1) strata,
and the rank theorem: the 12-dimensional Standard Model algebra can never be a centralizer in e6. Two U(1)s are forced.
§7–§9 (after reading the whole unification record): the record caught this first (B950/B952) and computed two ways past
it — the 27-VEV route (external point) and the holonomy hatch; A8 reproduces the hatch's centralizers, A7 the 18 exact
hypercharge directions at the weight level, and reads them: the rank-4 landing is trinification with colour broken to
SO(3)_C. The Higgs doublets are derived structure (10 ⊂ 27); the Higgs mechanism is external. `tests/test_audit.py`.
