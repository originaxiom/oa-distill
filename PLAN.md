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
| T13 | metallic fibre fields, exactly | From the rule side, which metallic rules σ_m give a conductor-3 fibre field? (closes T09's gap and T05's "label comes from the field") | m = 1..5, all four lift-sign twists | fixed Markoff points at 60 digits; minimal polynomials by integer relation; m = 1 reproduces x² − 3x + 3 | landed; m = 2 matches ℚ(i), m = 3 the octic; **m = 4 open**: rule-side candidates are a cubic-in-x² and an elliptic (non-geometric) point, the geometric one is not identified — needs the bundle's fibre generators |
| T14 | a measurable law of the rule's chain | How does the gap at label n open with the potential V? | labels n = ±1..±5; V ∈ [0.02, 0.4] | log-log slopes computed with residuals; the law stated only as computed; experiment class named | landed |
| T15 | the metallic chains | Do the metallic rules' chains carry labels in ℤ + ℤ·ω_m, ω_m = (√(m²+4) − m)/2? | m = 1, 2, 3 | labels within 3/N for each m; the wrong ω fails | landed |
| T16 | the 3d index at the fork | Same as T11, with an internal acceptance test that needs no literature value: invariance under re-triangulation of m004 and difference from m003 | randomized triangulations | series agree across triangulations to the computed order; differ for m003 | not started (research; only lands if the invariance test passes) |

Rule for adding rows: a row needs a question with a class and an acceptance test before any script is written.
