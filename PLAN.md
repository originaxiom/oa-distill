# Plan (executed in order; each module lands only when its acceptance test passes)

Principle: every module = statement + comparison class + script + frozen test. No module is written for a result it
does not reach. Status is updated in this file when a module lands or is blocked.

| # | module | question (one line) | comparison class | acceptance | status |
|---|---|---|---|---|---|
| T08 | rule spectrum | Do the rule's chain gaps carry the ℤ + ℤ/φ labels, and is the trace map the rule's Fricke action? | K4 rules vs the 4-letter rule | labels within 1/N; Fricke orbit = traces; invariant conserved | landed (on-site model; spectrum and trace map are one model) |
| T09 | the object's point | Where does the mapping torus sit on the rule's character variety? Is ℚ(√−3) derivable from the rule + the puncture alone? | metallic rules m = 1..6 | fixed curve of σ² computed symbolically; Markoff (κ = 0) gives x² − 3x + 3; a holonomy with a conjugating T (parabolic, commuting with [A,B]) exists numerically; field per m tabulated | landed |
| T10 | chain vs object level sets | Does any potential strength put the chain on the object's level set? Does any Dehn filling put the object on the chain's? | all V; fillings (p,q), \|p\|,\|q\| ≤ 8 | κ(chain) = 4 + 4V² ≥ 4 vs κ(object) = 0 proved; the set of fillings with real fibre-boundary trace computed | landed |
| T11 | 3d index of T[m004] vs T[m003] | Does the object-level quantum invariant differ for the sister? | m004 vs m003 | BLOCKED: no trusted reference value on this bench to validate an implementation against; do not land an unvalidated formula | blocked |
| T12 | locking predictions, pre-registered | Which signs does the dictionary predict before the measurements? | — | honest finding: the dictionary fixes relative signs only; a pre-registration must say so; no module until an absolute prediction exists | recorded in PREREGISTRATION.md; no module |

Rule for adding rows: a row needs a question with a class and an acceptance test before any script is written.
