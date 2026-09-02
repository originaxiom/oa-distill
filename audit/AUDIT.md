# Audit: the record's chain from E6 to "Standard Model structure", recomputed from scratch

**Date:** 2026-09-02. **Question (owner):** did the seat go all the way along the chain to the SM structure, before values?
Answer before this audit: no. This audit does it. Scripts `A1`–`A6` in this directory; every number below is from them.
Nothing here is a theorem module of this tree: the tree does not claim a gauge group (NOT_CLAIMED.md). This is an audit of
what the record's construction produces when rebuilt independently, with every choice named.

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
rank-6 algebra can never remove the two extra U(1)s. Removing them needs a mechanism outside this construction (a
Higgs sector, or a non-semisimple "measurement"), which the record does not have.

## Choice inventory for the chain E6 → 14

1. E6 itself (a field lookup: T05; constant on the object's commensurability class).
2. The sl(2) class carrying 2T (principal; e6 has 21 nilpotent orbits, hence many sl(2) classes; a different class gives
   a different centralizer of 2T).
3. The identification of the abstract 2T with a subgroup of SL(2) (up to conjugacy: unique).
4. Which point of C is "measured": the generic stratum gives su(3)+4U(1); one of six hyperplanes gives SM+2U(1); one of
   nine others gives A3+3U(1); codimension-2 strata give SU(5)+2U(1); lines give SO(10)+U(1). **Nothing in the object
   selects among these.** The record's "tuned first charge" landed on the SO(10) line.
5. Real form and compactness: the functionals on the real torus are complex (compact and split directions mix), and
   the compact real form of the result needs an external conjugation (the record's own B565/B1127).
6. Hypercharge, generations, chirality: outside this construction entirely (R46: 36 anomaly-free hypercharges per
   frame; B565: vector-like; B1033: generation count not derived).

## Verdict

Before values, the record's chain from E6 does reach a Standard-Model-shaped algebra, but only as one of several
Levi subalgebras of E6 sitting on six hyperplanes of a torus fixed by a field-level finite group, with two extra U(1)s
forced by rank, and with the object playing no part after the field. The construction is correct mathematics; its
endpoint is the standard E6 breaking pattern, chosen by hand at every branch. The seat's earlier statement that the
chain "leaves the object at the shape field" stands; this audit adds that even on the field's side, the SM algebra
itself is unreachable by the method, and SM+2U(1) is reachable only by a choice.
