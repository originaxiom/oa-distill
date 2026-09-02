# Statements
Each theorem's statement and comparison class, taken verbatim from its script's docstring.

## T01_rule_torsor.py

T01 -- THE RULE.  Statement: the substitution sigma: a->ab, b->a is one point of a free K4-torsor of equivalent
rules (reversal-conjugation and swap-conjugation); its word-reversal is an INNER automorphism; the fixed-point language
is closed under reversal (no reading arrow); the only intrinsic arrow is that the rule is not surjective on the free
monoid (the word 'bb' has no preimage).  Comparison class: the four Fibonacci-type rules on {a,b}.

## T02_double_tick.py

T02 -- THE DOUBLE TICK.  Statement: the incidence matrix F of sigma has det -1, so sigma is orientation-reversing on
the once-punctured torus; the mapping torus of sigma is the Gieseking manifold (non-orientable); the orientation double
cover of the Gieseking manifold is m004, the mapping torus of sigma^2 = A = [[2,1],[1,1]]; F commutes with A, so sigma acts
fibrewise on m004 as an orientation-reversing symmetry.  Comparison class: the one-tick and two-tick mapping tori.

## T03_symmetry_and_mirror.py

T03 -- SYMMETRY AND MIRROR.  Statement: Sym(m004) has order 8 (D4), four of its isometries reverse orientation
(cusp-map determinant -1), and the group is amphichiral; A -> A^-1 is realised by an orientation-preserving fibre map
(the bundle is invertible); the discrete faithful character (tr a, tr b, tr ab) lies in Q(sqrt-3)^3 and its complex
conjugate is realised by explicit endomorphisms a->u, b->v of pi_1 (relator -> +-I): the mirror acts on characters as
complex conjugation.  Comparison class: the isometries of m004; words of length <= 7.

## T04_closing_lattice.py

T04 -- THE CLOSING LATTICE.  Statement: (i) for any involutive symmetry iota of the object, every function T of the
object decomposes uniquely as T = T_even + T_odd with T_odd(iota x) = -T_odd(x); the sign of T_odd on 'the object' versus
'the mirror' is fixed only by naming one of them (a choice); (ii) the symmetries of m004 acting on the cusp form
(Z/2)^2 = <mirror, flow reversal>, so a scalar tracker carries one of four linear characters (plus the SL/PSL lift sign,
which is a symmetry of the representation); (iii) hence no object-canonical datum can orient m004.
Comparison class: all functions of the object; all isometries.

## T05_field_lookups.py

T05 -- THE FIELD LOOKUPS, STATED AS LOOKUPS.  Statement: the route (field discriminant -> conductor N -> SL(2,Z/N)
-> finite subgroup of SU(2)? -> McKay label) emits a label only for N in {1, 3, 5} (checked N <= 24); imaginary quadratic
conductors reach only N = 3, so on hyperbolic manifolds the route says E6 or nothing.  Among the first 1200 orientable
cusped census manifolds, 14 have shape field Q(sqrt-3); their volumes are 12, 24, 30 times the covolume of PSL(2,O_-3).
Comparison class: all conductors; the census.

## T06_object_beyond_field.py

T06 -- THE OBJECT BEYOND ITS FIELD.  Statement: within the 14-member Q(sqrt-3) class, H1, the cusp shape and the
number of covers by degree single out m004 (covers separate all 14); the sister m003's (1,q) closings differ from m004's
at every slope; the two 2-tetrahedron triangulations have different gluing matrices and different deformation curves,
each verified against a filled structure; hence the DGG theories T[m004], T[m003] differ.  Comparison class: the 14.

## T07_locking_dictionary.py

T07 -- THE LOCKING DICTIONARY.  Statement: a world's discrete choices (h: weak handedness, m: matter label, t: arrow)
mod CPT form a rank-2 space; a nonzero measured sign is CPT-even and of one of four types by (eP, eC, eT); the object's
closing lattice (T4, T6 mirror-odd; T7 = T3 flow-odd; no (odd,odd) axis) admits six dictionaries to {P, C, T}, of which
three pass the axes' own semantics and one (c = P, gamma5 = T) is populated by measured signs.  Predictions follow.
Comparison class: all dictionaries.
