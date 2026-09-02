# Audit: the record's chain from E6 to "Standard Model structure", recomputed from scratch

**Date:** 2026-09-02 (§1–§6); **§7–§9 added the same day after reading the whole unification record.**
**Question (owner):** did the seat go all the way along the chain to the SM structure, before values? — and, after §6 was
written: *"Higgs sector and 12 vs 14 dimensions is solved in the repo, no?"*
Scripts `A1`–`A8` in this directory; every number below is from them. Nothing here is a theorem module of this tree:
the tree does not claim a gauge group (NOT_CLAIMED.md). This is an audit of what the record's construction produces
when rebuilt independently, with every choice named — and, in §7–§9, of what the record already holds downstream of it.

## What was rebuilt

| step | record | audit | result |
|---|---|---|---|
| e6 itself | assumed | A1: built from the E6 root system with a Frenkel–Kac cocycle; 72 roots; Jacobi verified on 2000 random triples; Killing rank 78 | exact |
| the sl(2) carrying 2T | B854: principal | A2: principal sl(2) constructed (h = 2ρ∨, e = Σ e_α, f solved); e6 = V₂ ⊕ V₈ ⊕ V₁₀ ⊕ V₁₄ ⊕ V₁₆ ⊕ V₂₂ (dims 3, 9, 11, 15, 17, 23) | exact |
| 2T ⊂ SL(2) acting on e6 | B854 | A2: 24 elements; action by Sym^m(g) in string bases; homomorphism and bracket-preservation verified | exact |
| Cent(2T) | B854: u(1)⁴ | A2: one invariant in each of Sym⁸, Sym¹⁴, Sym¹⁶, Sym²², none in Sym², Sym¹⁰; abelian to 3·10⁻¹³ | **reproduced** |
| the ladder | B874: Cent(C) = 12, Cent(x8) = 30 | A3: Cent(C) = sl(3) ⊕ 4u(1) (dim 12, derived 8, centre 4); Cent(x8) = 30; other charges 24, 46, 40; pairs 12–30 | **reproduced** (and extended) |
| the "second measurement wall" | B892/B950: joint centralizer of dim 14 = su(3)⊕su(2)⊕u(1)³ at a complex wall point | A4/A5: on the record's pencils the joint centralizer is 12 generically and jumps to 18 (A3 ⊕ 3u(1)); the tuned first charge x8 + t*x16 has centralizer so(10) ⊕ u(1) (46) at three real t*. A6 then shows why: the 14-dimensional algebra is not a point but **six hyperplanes** of C | **reproduced as a stratum, not as a wall point** |
| the complete list of centralizers on C | not in the record | A6: root-hyperplane stratification of C (66 nonzero root functionals on 15 hyperplanes) | new |
| the non-abelian hatch's two landings | B1098: A2 class → su(3)⊕su(3) (16, rank 4); minimal A1 → su(6) (35, rank 5) | A8: both centralizers recomputed on A1's e6; the A2 centralizer is spanned by the other two trinification su(3)'s; the principal sl(2) of e6 has centralizer 0 | **reproduced** |
| the exact hypercharge at the A2 landing | B1102: 18 rational directions carry the SM 6Y multiset; none commutes with a full su(3); B1118: two S₃×S₃ orbits fused by the mirror | A7: reproduced at the weight level with a one-line reason for side 2 (a direction pure on one factor has ≥ 9 zero eigenvalues; the target has 2) | **reproduced, and read** (§8) |

## A6 — every centralizer reachable by "measuring" on the torus

C = ⟨x8, x14, x16, x22⟩ consists of semisimple elements; for x ∈ C, Cent(x) = 𝔥′ ⊕ Σ_{α(x)=0} 𝔤_α, so
dim Cent(x) = 6 + #{roots vanishing at x}. The strata of the arrangement of the 15 hyperplanes:

| codimension in C | (dim, #roots) of the centralizer, with the number of index tuples producing it | reading |
|---|---|---|
| 0 (generic point) | (12, 6) | A2 ⊕ 4u(1) = su(3) + four U(1)s |
| 1 (a hyperplane) | (14, 8) × 6, (18, 12) × 9 | **A2⊕A1 ⊕ 3u(1) = su(3)+su(2)+3U(1) on six hyperplanes**; A3 ⊕ 3u(1) on nine |
| 2 | (16, 10) × 9, (18, 12) × 6, (20, 14) × 18, (26, 20) × 54, (30, 24) × 18 | includes A4 ⊕ 2u(1) = su(5)+2U(1) |
| 3 (a line) | (20, 14), (28, 22), (36, 30), (46, 40) | includes D5 ⊕ u(1) = so(10)+U(1) |

This is the textbook E6 ⊃ SO(10)×U(1) ⊃ SU(5)×U(1)² ⊃ SM×U(1)² breaking pattern, as it must be: the 2T-torus is a
4-dimensional subtorus of a Cartan subalgebra, and centralizers of torus elements are Levi subalgebras.

## The rank theorem (why "14, not 12" was never a near miss)

The centralizer of any semisimple element of e6 contains a Cartan subalgebra and therefore has rank 6. The Standard
Model algebra su(3) ⊕ su(2) ⊕ u(1) has rank 4. **So no point of C, and no point of any torus in e6, has the 12-dimensional
Standard Model algebra as its centralizer.** The smallest SM-containing centralizer is su(3) ⊕ su(2) ⊕ u(1)³, dimension 14,
rank 6, and it occurs on six hyperplanes of C. The record's B950 correction ("14, not 12; two U(1)s left over") is
therefore not a computational refinement but a structural necessity of the method: measuring by centralizers in a
rank-6 algebra can never remove the two extra U(1)s. Removing them needs a mechanism outside this construction.
*(The record proved this same theorem on 2026-08-08, B952, and names the mechanism: §7.)*

## Choice inventory for the chain E6 → 14

1. E6 itself (a field lookup: T05; constant on the object's commensurability class).
2. The sl(2) class carrying 2T (principal; e6 has 20 nonzero nilpotent orbits, hence 20 sl(2) classes; a different
   class gives a different centralizer of 2T).
3. The identification of the abstract 2T with a subgroup of SL(2) (up to conjugacy: unique).
4. Which point of C is "measured": the generic stratum gives su(3)+4U(1); one of six hyperplanes gives SM+2U(1); one of
   nine others gives A3+3U(1); codimension-2 strata give SU(5)+2U(1); lines give SO(10)+U(1). **Nothing in the object
   selects among these.** The record's "tuned first charge" landed on the SO(10) line; its second charge is chosen by a
   perfect-square gate so that the plane through the first charge lies in one of the six hyperplanes.
5. Real form and compactness: the functionals on the real torus are complex (compact and split directions mix), and
   the compact real form of the result needs an external conjugation (the record's own B565/B1127/B1134).
6. Hypercharge, generations, chirality: outside this construction (B1160/B1170: anomaly arithmetic on an SM-shaped
   15-plet, arena-generic; B565/B1086: vector-like on any closed assembly; B1033/B1161: one generation per 27).

## Verdict on the chain E6 → 14 (§1–§6, unchanged)

Before values, the record's chain from E6 does reach a Standard-Model-shaped algebra, but only as one of several
Levi subalgebras of E6 sitting on six hyperplanes of a torus fixed by a field-level finite group, with two extra U(1)s
forced by rank, and with the object playing no part after the field. The construction is correct mathematics; its
endpoint is the standard E6 breaking pattern, chosen by hand at every branch. The seat's earlier statement that the
chain "leaves the object at the shape field" stands; this audit adds that even on the field's side, the SM algebra
itself is unreachable by the method, and SM+2U(1) is reachable only by a choice.

**Correction to the first version of this verdict.** The first version ended "…which the record does not have."
That was written before the record's downstream work was read whole, and it is wrong as stated: the record names the
mechanism (B964), computes two candidates for it (§7), and prices what each leaves open. The corrected statement is §9.

---

## §7. What the record already holds downstream of the 14 — read whole, 2026-09-02

Read myself, in full: `docs/THE_SM_VERDICT.md`, `SM_SPECIFICATION_LEDGER.md`, `GUT_REQUIREMENTS_LEDGER.md`,
`THE_CLAIM.md`, `THE_END_TO_END_CHAIN.md`, `THE_LADDER.md`, `THE_FRAMEWORK.md`, `WHAT_WOULD_COUNT.md`,
`GRAND_COMPUTATION_v0.md` + `_LEDGER.md`, `STRUCTURE_TO_NATURE_MASTERPLAN.md`, `PRICED_DOORS.md`, `MASTERPLAN.md`,
`STRATEGIC_SYNTHESIS.md`, `THE_CAMPAIGN.md`, `THE_PATTERN_MEDITATION.md`, `LITERATURE_GAUGE_SM_2026-07-13.md`,
`UNIFIED_STATE.md`, the P3 paper spec, claim pool and tex, the structure-paper skeleton and abstract, PC26/PC27, the
outside-bench chain documents, cc3's Part-0 correction, and the arcs B862–B864, B876, B884, B892, B897, B928,
B950–B955, B959–B964, B970, B978, B987, B992, B994, B1033, B1092–B1102, B1112, B1118, B1119, B1135, B1145,
B1160–B1162, B1170, B1185, B1195, B1205–B1208, B1216, B1220, B1225, B1226, B1229–B1240. (The reading ledger, one line
per document, is in the origin-axiom seat report `READING_2026-09-02_unification_record.md`.)

### 7.1 The record's own position on 12 vs 14

| item | where | status |
|---|---|---|
| 14 ≠ 12 caught | B950 (2026-08-08); B892's claim line corrected by B1237 (2026-09-02) | done |
| the rank obstruction as a theorem | B952 §D (= this audit's rank theorem) | done, same theorem |
| the two unsheddable units are U(1)_ψ, U(1)_χ; u(1)³ = span(Y, χ, ψ) | B953, B992 | computed |
| the cascade is an adjoint Higgs mechanism; the object supplies the rank-preserving (adjoint) VEVs and lacks the rank-reducing 27 VEV | B964 | the missing piece, named |
| toral, finite-image, abelian-holonomy and adjoint-form routes to rank 4 | B955, B959, B960, B1079/B1094 | all closed (agrees with A6) |
| the 27-VEV route | B962 (direction an input everywhere; F₄ = generic-VEV stabilizer); B1025 (the directions are canonical multiplicity-one lines: the SO(10)-singlet of the 27, the SU(5)-singlet ν^c of the 16); B1092 (the second VEV must be a pure spinor: an 11-dim cone, Spin(10)-transitive — a condition, not a point); B990/B1225 (no invariant readable off the object picks a point of its own orbit); Route A: h = h⁺ = 1 proved (B1093), coarse orbit unique at det 5 (B1099), the K-refined orbit count is frontier mathematics | OPEN at the point level; NEEDS-SPECIALIST |
| the non-abelian hatch | B1094 (named), B1098 (20 sl(2) classes; A2 → su(3)⊕su(3) rank 4; A1 → su(6) rank 5; 2A1 excluded), B1100 (27 complex at A2), B1102 (18 exact hypercharge directions, none colour-commuting), B1112 (A2 the unique projective SM-compatible landing; A1 needs the spin lift), B1118 (the 18 = two orbits fused by the mirror), B1236 (A1: SM-shaped 27 at exact multiplet grade, Y unique within support, the extra U(1) exhibited) | rank 4 reached; product structure not |
| the record's own summary sentence | GUT ledger §D (corrected R48-F2): "the MEASUREMENT lane is rank-6 by theorem; the object's own non-abelian lane reaches rank 4; what no lane yet supplies is the color-commuting product with the exact values"; THE_FRAMEWORK: "two steps from the SM's own twelve, not zero" | honest |

### 7.2 The record's own position on the Higgs sector

| item | where | status |
|---|---|---|
| the Higgs doublets | the (1,2)⊕(1,2) of the 10 ⊂ 27 (B884 grading; B987 dissolves ladder rung X8); one generation's 27 = 16 + 10 + 1 carries H_u, H_d, D, D^c, S (B970/B978; B1236 at exact multiplet grade) | derived structure — standard E6 content, reproduced |
| Yukawa support | the unique cubic in Sym³(27); 11 coupled cells = 16·16·10 ⊕ 1·10·10, every zero charge-forced (B884); no adjoint VEV can mass a 27 fermion, 78 ∉ 27⊗27 (B978) | derived structure |
| the GUT-scale (rank-reducing) Higgs | external; see the 27-VEV row above | open / external |
| doublet–triplet splitting | external, needs a colour choice (B298/B299); GUT ledger row 6 ABSENT | not addressed |
| the electroweak Higgs line | the down-readout's ℙ³ = ℙ(B₀): up to three continuous observer parameters; the cut ledger runs 3 → 2 → 1 and stops **one linear condition short** of points (B1205/B1206); the three named candidates negative and their space closed (B1208, B1220: CLOSED-PERMANENT at current knowledge); one live overturner, the lepton-leg ℤ/12 character fork, pending codex R030; the consistency turn names a finiteness test via Cardy states (B1229) | one condition short |
| EWSB, ⟨H⟩, m_H | dimensionful: the reader's by the scale wall (B1226 box C; B1012/B1088); every hatch arc: "EWSB remains outside" | external by theorem |
| Yukawa values and textures | three genuinely distinct suppression mechanisms (B1185): μ_u = 0 on the heterotic dressing but 6 nonzero components on the object channel; the down block an object-forced 3×3×4 tensor with an exact skew zero; the Family Rank Theorem (rank 2, kernel = the Higgs's own family, E8 fence); the hierarchy carried by the twist D₂ (B923/B928); the blind shape sheet vs measured mixing: MISS (outside bench, Stage 2 unsealed) | structure only; values withheld |
| "Higgs-synthesis OP theorem" | B859: ½[j(τ) − j(−4/τ)] = 818626500√3 at the cusp shape — a modular identity carrying the word "Higgs" as a label, not a Higgs mechanism | label only |

## §8. Two things the audit adds to the record's downstream work

**8.1 The hatch's rank drop is bought by breaking colour (A7, read physically).** At the A2 landing the 27 restricts
to su(3)₁ ⊕ su(3)₂ as (3, 3̄) ⊕ 3·(3̄, 1) ⊕ 3·(1, 3), the "3·" being the eaten factor's weights. A7 finds the 18
exact-hypercharge directions and shows every one mixes both surviving Cartans (a direction pure on one factor gives
≥ 9 zero eigenvalues; the SM multiset has 2). Read with the branching, the representative Y with
a = (−1/6, −1/6, 1/3) on su(3)₁ and b = (1/3, 1/3, −2/3) on su(3)₂ gives: on the (3, 3̄) nonet the charges
{−1/2 ×4, 1/2 ×2, 0 ×2, 1 ×1} — the lepton–Higgs nonet (L, H_d, H_u, ν^c, S, e^c); on the 3·(3̄, 1) the charges
{1/6 ×6, −1/3 ×3} — the quark doublet and D, **with the triple multiplicity as the colour index**; on the 3·(1, 3) the
charges {1/3 ×6, −2/3 ×3} — d^c, D^c, u^c. So the exact SM hypercharge at this landing forces the *eaten* factor to be
colour, and the eaten factor is broken to the real so(3) of its principal sl(2): **the A2 landing is trinification
SU(3)_L × SU(3)_R with SU(3)_C → SO(3)_C.** That is why no colour-commuting hypercharge exists there (B1102 side 2),
and it answers B1102's follow-up (ii): no identification of colour inside the centralizer is possible, because colour
is the factor the object's holonomy ate. The rank drop 6 → 4 at this stratum sheds two Cartan directions of colour,
not U(1)_ψ and U(1)_χ. Follow-up (i) is answered by B1098's own table: the only rank-4 centralizers among the twenty
sl(2) classes are a₂⊕a₂ (this one) and b₃⊕u(1) (excluded by B1098); the A1 stratum keeps the full SM product (B1236)
at rank 5, with one extra U(1) and one spin bit. **So within the sl(2)-factored hatch there is no stratum with an
unbroken su(3)_C ⊕ su(2)_L ⊕ u(1)_Y at rank 4.**

**8.2 The fourteen is a stratum, not a wall point (A6).** The record's language "the wall point y*" describes one
point on one of six hyperplanes of C, selected by a perfect-square gate; the 14-dimensional algebra is the generic
centralizer on each of those hyperplanes. This changes no banked number; it changes what the selection is a selection
of (a hyperplane, not a point) and makes the "second measurement" one of six equivalent choices.

## §9. Reconciled verdict: are "12 vs 14" and "the Higgs sector" solved in the record?

**12 vs 14 — located and priced, not derived.** The record found the gap first, proved why it is structural (B952 =
the rank theorem here), named the missing piece exactly (the rank-reducing 27 VEV, B964), and computed two candidate
mechanisms. Neither delivers su(3)_C ⊕ su(2)_L ⊕ u(1)_Y from the object: the 27-VEV route delivers the product
structure but needs a point on an 11-dimensional cone the object provably cannot select (B990/B1092/B1225; Route A's
counter is open mathematics), and the holonomy hatch delivers rank 4 only by breaking colour (§8.1), or keeps the
product at rank 5 with an extra U(1) and a spin bit (A1). In the record's own words (GUT ledger §D, THE_FRAMEWORK):
two steps from the SM, not zero. This audit agrees and adds the colour reading.

**The Higgs sector — the representation is derived structure; the mechanism is external, part by theorem.** The
doublets sit in the 10 ⊂ 27 with a charge-forced Yukawa support (B884/B987): standard E6 content, correctly labelled
"reproduced, not predicted" in the record. Everything that makes it a *mechanism* — the GUT VEV point, doublet–triplet
splitting, the electroweak Higgs line (one condition short), the VEV and mass (dimensionful, hence the reader's by the
scale wall) — is outside the object, and the record says so at every hatch. The one live lever the record names is the
lepton-leg character fork (B1208 (b)), which would close the Higgs line to a finite point set, not to a prediction.

**On "you should be upgrading the project."** The upgrade this audit can honestly offer is: (i) the stratification (A6)
replacing "wall point"; (ii) the colour reading of the hatch (§8.1), which turns B1102's obstruction into a physics
statement and closes two of its three follow-ups; (iii) independent recomputation of the two hatch centralizers and the
18 hypercharge directions (A7, A8). What it cannot offer is a derivation the record lacks: the SM product at rank 4
needs the 27-VEV point, and the record has proved that point is not the object's to give (B1225). The honest name for
the downstream chain is standard E6 grand-unified model-building — adjoint VEVs plus a 27 VEV — with the object
supplying the field-level E6 and the arena, exactly as `THE_END_TO_END_CHAIN.md` links 4, 8 and 11 already state.
