# Comparator reference statement versus upstream Formal Conjectures (checked 2026-09-09)

Upstream: google-deepmind/formal-conjectures, FormalConjectures/Millenium/NavierStokes.lean at commit 8bf45ed70d48b2b2a501de9c00b26bfa38c573ee (296 lines).
Repository copy: openai/NavierStokesAndEuler, ComparatorChallenges/NavierStokes.lean (286 lines).

`diff -w` shows only: import changed from FormalConjecturesUtil to Mathlib; namespace NavierStokes -> NavierStokes.Comparator; two local notations inlined; `@[category ...]`/`AMS` metadata attributes removed; the (A)/(B) theorems and the doc-comment sentences about errata removed; header comment added. Every definition, including `integrable : forall t >= 0, MemLp (norm (v . t)) 2`, `globally_bounded_energy`, `isOnePeriodic_pressure`, the decay structures and the derivWithin encoding, is byte-identical modulo whitespace. Upstream's doc comment (lines 49-53) states: "The Clay PDF also includes errata; in particular, we include spatial 1-periodicity of the pressure in the periodic case."
