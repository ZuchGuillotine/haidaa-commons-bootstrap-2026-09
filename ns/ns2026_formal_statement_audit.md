# Audit: does the Lean formal statement faithfully encode Fefferman's (C) and (D)?

Repository audited: `/private/tmp/claude-501/-Users-benjamincox-dsm-haidaa/da69bd54-2977-4882-88ec-9babfab4e89d/scratchpad/ns-paper/lean-repo`
(single squashed git commit `8937a8f`, so there is no history to diff against.)

Toolchain: `leanprover/lean4:v4.34.0-rc2`, Mathlib pinned at `v4.34.0-rc2`, Comparator pinned at
`v4.34.0-rc2` (`lakefile.toml`). No build was run for this audit (a build is running separately);
every claim below is from source reading and grep.

Scope note: the reference statement lives in
`/.../lean-repo/ComparatorChallenges/NavierStokes.lean` (286 lines) and is duplicated verbatim,
minus the two theorems, in `/.../lean-repo/NavierStokes/ComparatorDefinitions.lean` (245 lines).
A byte diff of the two definition blocks (`ComparatorDefinitions.lean:34-245` vs
`ComparatorChallenges/NavierStokes.lean:61-272`) shows **exactly one** difference: the challenge
file lacks the `open ContDiff Set InnerProductSpace MeasureTheory` line at that offset because it
appears earlier in the file (challenge line 61 is the `open` line itself; both files have it).
Otherwise the definition text is identical character for character. This matters: Comparator's
check is structural equality of the elaborated constants, so the two copies must agree.

---

## 1. The exact Lean definitions

All line numbers below are from `ComparatorChallenges/NavierStokes.lean` unless stated; the
identical text sits in `NavierStokes/ComparatorDefinitions.lean` at the offsets noted in
parentheses.

### Ambient notation (`ComparatorChallenges/NavierStokes.lean:71-76`, defs `:40-43`)

```lean
local notation "ℝ^" n:65 => EuclideanSpace ℝ (Fin n)
local notation "ℝ³" => EuclideanSpace ℝ (Fin 3)
variable {n : ℕ}
```

with `open ContDiff Set InnerProductSpace MeasureTheory` and `open scoped Laplacian` at
lines 61-62 (defs `:34-35`).

`ℝ³ = EuclideanSpace ℝ (Fin 3)` carries the Euclidean inner product and, via its `MeasureSpace`
instance, Lebesgue measure. Both millennium theorems are stated at `ℝ³`; the definitions are
dimension-generic in `n`.

### Divergence (`:80-81`, defs `:54-55`)

```lean
noncomputable
def divergence (v : ℝ^n → ℝ^n) (x : ℝ^n) : ℝ := (fderiv ℝ v x).trace ℝ (ℝ^n)
```

with local notation `∇⬝`. The docstring at `:74-78` explicitly acknowledges the junk value:
"If `v` is not differentiable at `x`, then `fderiv` is the zero map, so this definition has the
corresponding junk value 0." A `@[simp]` lemma `divergence_of_not_differentiableAt` (`:85-88`)
records this.

The project separately proves this coincides with the coordinate definition
`∑ᵢ ∂ᵢ vᵢ` in `NavierStokes/ComparatorBridge.lean:28-32` (`divergence_eq`), via
`LinearMap.trace_eq_sum_inner` on `EuclideanSpace.basisFun`. So the trace encoding is the
mathematically intended divergence.

### Periodicity (`:118-122`, defs `:95-96`)

```lean
def IsOnePeriodic {α : Sort*} (f : ℝ^n → α) : Prop :=
  ∀ x i, f (x + EuclideanSpace.single i 1) = f x
```

Quantifying over **all** `x` makes the negative shifts follow, so this is invariance under the
full lattice ℤⁿ with period 1 — Fefferman's `u°(x + eⱼ) = u°(x)`.

### `InitialVelocityCondition` and `InitialVelocityConditionDecay` (`:134-155`, defs `:108-129`)

```lean
structure InitialVelocityCondition (u₀ : ℝ^n → ℝ^n) : Prop where
  div_free : ∀ x, ∇⬝ u₀ x = 0
  smooth : ContDiff ℝ ∞ u₀

structure InitialVelocityConditionDecay (u₀ : ℝ^n → ℝ^n) : Prop extends
    InitialVelocityCondition u₀ where
  decay : ∀ m : ℕ, ∀ K : ℝ, ∃ C : ℝ, ∀ x, ‖iteratedFDeriv ℝ m u₀ x‖ ≤ C / (1 + ‖x‖) ^ K
```

Divergence-freeness of `u₀` **is** required (`div_free`, inherited). `ContDiff ℝ ∞` with
`open ContDiff` in scope means `∞ = ((⊤ : ℕ∞) : WithTop ℕ∞)`
(`Mathlib/Analysis/Calculus/ContDiff/FTaylorSeries.lean:120`), i.e. genuine C^∞ (all finite
orders), **not** the analytic order `ω = (⊤ : WithTop ℕ∞)`. This is the correct reading and is
called out in the project's own note at `NavierStokes/ProblemStatement.lean:18-19`.

`(1 + ‖x‖) ^ K` with `K : ℝ` is `Real.rpow` on a base `≥ 1`, so `C / (1+‖x‖)^K = C (1+‖x‖)^{-K}`.

### `ForceCondition` / `ForceConditionDecay` / `ForceConditionPeriodic` (`:175-206`, defs `:149-179`)

```lean
structure ForceCondition (f : ℝ^n → ℝ → ℝ^n) : Prop where
  smooth : ContDiffOn ℝ ∞ (↿f) (Set.univ ×ˢ Set.Ici 0)

structure ForceConditionDecay (f : ℝ^n → ℝ → ℝ^n) : Prop extends ForceCondition f where
  decay : ∀ m : ℕ, ∀ K : ℝ, ∃ C : ℝ, ∀ x, ∀ t ≥ 0,
    ‖iteratedFDerivWithin ℝ m (↿f) (Set.univ ×ˢ Set.Ici 0) (x, t)‖ ≤ C / (1 + ‖x‖ + t) ^ K

structure ForceConditionPeriodic (f : ℝ^n → ℝ → ℝ^n) : Prop extends ForceCondition f where
  isOnePeriodic : ∀ t ≥ 0, IsOnePeriodic (f · t)
  decay : ∀ m : ℕ, ∀ K : ℝ, ∃ C : ℝ, ∀ x, ∀ t ≥ 0,
    ‖iteratedFDerivWithin ℝ m (↿f) (Set.univ ×ˢ Set.Ici 0) (x, t)‖ ≤ C / (1 + t) ^ K
```

`↿f` is `Function.uncurry`, so the domain is `ℝ^n × ℝ` with the *space* coordinate first;
`Set.univ ×ˢ Set.Ici 0` is therefore ℝⁿ × [0,∞). Correct orientation.

### `NavierStokesExistenceAndSmoothness` (`:217-241`, defs `:191-209`)

```lean
structure NavierStokesExistenceAndSmoothness
    (nu : ℝ) (u₀ : ℝ^n → ℝ^n) (f : ℝ^n → ℝ → ℝ^n)
    (v : ℝ^n → ℝ → ℝ^n) (p : ℝ^n → ℝ → ℝ) : Prop where
  navier_stokes : ∀ x, ∀ t ≥ 0,
    derivWithin (v x ·) (Set.Ici 0) t + fderiv ℝ (v · t) x (v x t) =
      nu • Δ (v · t) x - gradient (p · t) x + f x t
  div_free : ∀ x, ∀ t ≥ 0, ∇⬝ (v · t) x = 0
  initial_condition : ∀ x, v x 0 = u₀ x
  velocity_smooth : ContDiffOn ℝ ∞ (↿v) (Set.univ ×ˢ Set.Ici 0)
  pressure_smooth : ContDiffOn ℝ ∞ (↿p) (Set.univ ×ˢ Set.Ici 0)
```

- **Time derivative**: `derivWithin (v x ·) (Set.Ici 0) t`, i.e. the derivative of the time slice
  `s ↦ v x s` relative to `[0,∞)`. At `t = 0` this is the one-sided derivative; no extension to
  negative time is differentiated. Documented at `:56-59` of the challenge header.
- **Advection**: `fderiv ℝ (v · t) x (v x t)` = `Dₓv(x,t)[v(x,t)]` = `(v·∇)v`. Correct.
- **Laplacian**: `Δ` is Mathlib's `InnerProductSpace.instLaplacian`
  (`Mathlib/Analysis/InnerProductSpace/Laplacian.lean:137-139`), defined as the second
  *unrestricted* iterated Fréchet derivative applied to the canonical covariant tensor; the
  orthonormal-basis formula `Δ f = fun x ↦ ∑ i, iteratedFDeriv ℝ 2 f x ![v i, v i]`
  (`Laplacian.lean:171-175`) makes it `∑ᵢ ∂ᵢ∂ᵢ`, applied componentwise to the ℝ³-valued slice.
  The project's `laplacian_eq` (`NavierStokes/ComparatorBridge.lean:41-53`) proves it equals its
  own coordinate Laplacian.
- **Pressure gradient**: `gradient (p · t) x` = `InnerProductSpace.gradient`, matched to the
  project's coordinate gradient in `ComparatorBridge.lean:34-39`.
- **Smoothness**: `ContDiffOn ℝ ∞` on ℝ³ × [0,∞) for *both* `v` and `p`. Fefferman (6)/(11)
  require exactly `p, u ∈ C^∞(ℝ³ × [0,∞))`. Pressure smoothness **is** required.
- The equation is imposed at `t = 0` as well (`∀ t ≥ 0`), matching Fefferman's (1) on `[0,∞)`.

### `NavierStokesExistenceAndSmoothnessRn` (`:245-254`, defs `:219-227`)

```lean
structure NavierStokesExistenceAndSmoothnessRn
    (nu : ℝ) (u₀ : ℝ^n → ℝ^n) (f : ℝ^n → ℝ → ℝ^n)
    (v : ℝ^n → ℝ → ℝ^n) (p : ℝ^n → ℝ → ℝ) : Prop
  extends NavierStokesExistenceAndSmoothness nu u₀ f v p where
  integrable : ∀ t ≥ 0, MemLp (‖v · t‖) 2
  globally_bounded_energy : ∃ E, ∀ t ≥ 0, (∫ x : ℝ^n, ‖v x t‖ ^ 2) < E
```

`MemLp (fun x => ‖v x t‖) 2` w.r.t. the default `volume` on `EuclideanSpace ℝ (Fin n)`; the
integral is the Bochner/Lebesgue integral `MeasureTheory.integral volume`.

### `NavierStokesExistenceAndSmoothnessPeriodic` (`:262-269`, defs `:236-243`)

```lean
structure NavierStokesExistenceAndSmoothnessPeriodic ... : Prop
  extends NavierStokesExistenceAndSmoothness nu u₀ f v p where
  isOnePeriodic_velocity : ∀ t ≥ 0, IsOnePeriodic (v · t)
  isOnePeriodic_pressure : ∀ t ≥ 0, IsOnePeriodic (p · t)
```

with the docstring "The pressure is also required to be 1-periodic, following the errata appended
to the Clay problem statement" (`:258-260`). No energy condition — correct for (D).

### The two theorems (`:272-285`)

```lean
theorem navier_stokes_breakdown_R3 (nu : ℝ) (hnu : nu > 0) :
    ∃ (u₀ : ℝ³ → ℝ³) (f : ℝ³ → ℝ → ℝ³),
    InitialVelocityConditionDecay u₀ ∧ ForceConditionDecay f ∧
    ¬ (∃ v p, NavierStokesExistenceAndSmoothnessRn nu u₀ f v p) := by
  sorry

theorem navier_stokes_breakdown_periodic (nu : ℝ) (hnu : nu > 0) :
    ∃ (u₀ : ℝ³ → ℝ³) (f : ℝ³ → ℝ → ℝ³),
    InitialVelocityConditionPeriodic u₀ ∧ ForceConditionPeriodic f ∧
    ¬ (∃ v p, NavierStokesExistenceAndSmoothnessPeriodic nu u₀ f v p) := by
  sorry
```

The submission `NavierStokes/ComparatorSolution.lean:16-27` restates these **character for
character** (same binder names `nu`, `hnu`, `u₀`, `f`, `v`, `p`) and discharges each with a single
`exact` into `ComparatorBridge.navier_stokes_breakdown_R3` / `..._periodic`
(`NavierStokes/ComparatorR3Theorem.lean:36-42` and `NavierStokes/ComparatorTheorem.lean:48-54`).
Lines 31-32 of `ComparatorSolution.lean` are `#print axioms` for both names.

`ComparatorSolution.lean` imports only `NavierStokes.ComparatorR3Theorem` and
`NavierStokes.ComparatorTheorem`; `ComparatorBridge.lean:1-2` imports only
`NavierStokes.ComparatorDefinitions` and `NavierStokes.ProblemStatement`. **No file in the
repository imports `ComparatorChallenges`** (verified: `grep -rn "import ComparatorChallenges"`
returns nothing), so the sorried reference is not in the submission's import closure.

---

## 2. Where the encoding is weaker, stronger, or different from Fefferman

### Junk values — none of them bite, and this is provable rather than accidental

Every operator in the statement (`fderiv`, `derivWithin`, `Δ`, `gradient`, the Bochner `∫`) has a
junk value in Mathlib. The question is whether any junk value lets a condition be satisfied
vacuously. It does not, because each structure carries a smoothness/integrability field that
forces the honest value:

| Operator | Junk value | Neutralized by |
|---|---|---|
| `∇⬝ u₀` (`fderiv`) | 0 | `InitialVelocityCondition.smooth : ContDiff ℝ ∞ u₀` |
| `∇⬝ (v · t)`, `fderiv ℝ (v · t) x` | 0 | `velocity_smooth`; the space slice is globally `ContDiff ℝ ∞` for every `t ≥ 0` (this is proved in `ComparatorBridge.lean:176-181`, `smooth_space_slice`) |
| `Δ (v · t) x` (unrestricted second `fderiv`) | 0 | same — the slice is smooth on all of ℝ³ even at `t = 0`, so the *unrestricted* Laplacian is the genuine one |
| `gradient (p · t) x` | 0 | `pressure_smooth`, same slice argument |
| `derivWithin (v x ·) (Ici 0) t` | 0 | `velocity_smooth` gives `DifferentiableWithinAt` on `Ici 0`; at `t = 0` this is the honest one-sided derivative |
| `∫ x, ‖v x t‖ ^ 2` | **0 when non-integrable** | `integrable : ∀ t ≥ 0, MemLp (‖v · t‖) 2` |

**The `∫` convention specifically.** This is the one place where the convention could have
mattered, and the statement handles it correctly. Mathlib's Bochner integral of a non-integrable
function is `0`. Had `globally_bounded_energy` stood alone, a smooth solution with infinite energy
would satisfy `0 < E` and hence *count as a solution*; the solution class would be strictly larger
and `¬∃` would be a **strictly stronger** claim than Fefferman's (C) — the theorem would then
assert non-existence of *any* global smooth solution regardless of energy, which is not what (C)
says and is a claim the proof does not support (the proof genuinely uses finite energy). The
`integrable` field removes that reading: the class is exactly "smooth solutions whose velocity is
in L² at each time, with a uniform bound on the L² norm", which is precisely Fefferman's (7)
(a finite Lebesgue integral of a non-negative measurable function is the same thing as
integrability). So the convention cuts in neither direction here, and **the bounded-energy
hypothesis is not vacuous**.

That `integrable` is load-bearing rather than decorative is visible in the proof:
`NavierStokes/R3FiniteEnergyComparison.lean:28-29` converts `h.integrable` into
`Integrable (fun x => ‖v (t,x)‖ ^ 2)` via `memLp_two_iff_integrable_sq_norm` before the uniform
bound at `:30-33` is used, and both feed
`NavierStokesR3.WholeSpaceUniqueness.classical_uniqueness_on_Icc` at `:59-62`.

### Item-by-item comparison with Fefferman

| Fefferman | Lean | Verdict |
|---|---|---|
| (4) `\|∂^α u°(x)\| ≤ C_{αK}(1+\|x\|)^{-K}` ∀α,K | `∀ m K, ∃ C, ∀ x, ‖iteratedFDeriv ℝ m u₀ x‖ ≤ C/(1+‖x‖)^K` | **at least as strong**. `‖iteratedFDeriv ℝ m u₀ x‖` is the operator norm of the full order-`m` jet, and `\|∂^α u°\| ≤ ‖iteratedFDeriv‖` for every `\|α\| = m` (unit basis vectors). Lean also uses one constant per order rather than per multi-index, and quantifies `K` over all reals rather than integers. Since (4) is on the *existential* side, "at least as strong" means the theorem is at least as strong as (C). |
| (5) `\|∂_x^α ∂_t^m f\| ≤ C(1+\|x\|+t)^{-K}` | `‖iteratedFDerivWithin ℝ m (↿f) (univ ×ˢ Ici 0) (x,t)‖ ≤ C/(1+‖x‖+t)^K` for `t ≥ 0` | **at least as strong**, same argument; the order-`m` total jet in `(x,t)` dominates every mixed partial with `\|α\| + j = m`. |
| (6) `p, u ∈ C^∞(ℝ³ × [0,∞))` | `velocity_smooth` + `pressure_smooth`, `ContDiffOn ℝ ∞` on `univ ×ˢ Ici 0` | **faithful**. Pressure smoothness is required. `∞` is C^∞, not analytic. |
| (7) `∫\|u(x,t)\|² dx < C` ∀t | `integrable` + `globally_bounded_energy`, `∀ t ≥ 0` | **faithful** (see above). Required for all `t ≥ 0`, uniform constant. |
| (1) on `ℝ³ × [0,∞)` | `navier_stokes : ∀ x, ∀ t ≥ 0` | **faithful**, including `t = 0` with a one-sided time derivative. |
| (2) `div u = 0` | `div_free : ∀ x, ∀ t ≥ 0` | faithful. |
| (3) `u(x,0) = u°(x)` | `initial_condition : ∀ x, v x 0 = u₀ x` | faithful. |
| `u°` smooth, divergence-free | `InitialVelocityCondition` (both fields inherited by both variants) | faithful — divergence-freeness of `u₀` **is** required in (C) and in (D). |
| (8) periodicity of `u°`, `f` | `InitialVelocityConditionPeriodic.isOnePeriodic`, `ForceConditionPeriodic.isOnePeriodic` | faithful. |
| (9) `\|∂_x^α ∂_t^m f\| ≤ C(1+t)^{-K}` | periodic `decay` field | faithful (no spatial decay demanded, correctly). |
| (10) periodicity of `u` | `isOnePeriodic_velocity` | faithful. |
| (11) `p, u ∈ C^∞` | inherited `velocity_smooth`/`pressure_smooth` | faithful. |
| erratum: `p` periodic | `isOnePeriodic_pressure` | faithful to the **corrected** (D) — see caveat below. |
| ν a given positive coefficient | `(nu : ℝ) (hnu : nu > 0)`, universally quantified | **stronger**: one statement for every positive ν. |

### The three differences worth naming explicitly

1. **Pressure periodicity in (D) is a real restriction on the solution class.** Requiring
   `isOnePeriodic_pressure` *shrinks* the set of competitors, so `¬∃` is formally weaker than the
   pre-erratum reading of (D) (which would also have to exclude solutions with, e.g., a
   linear-in-x pressure). This matches the Clay erratum and matches upstream formal-conjectures'
   stated intent, so it is correct — but a reader should understand the theorem as
   "erratum-corrected (D)", not as the literal text of the original PDF.

2. **The witnesses are `u₀ = 0` with a non-trivial force.** `ComparatorTheorem.lean:38` and
   `ComparatorR3Theorem.lean:30-31` both instantiate `u₀ := fun _ => 0`. That is entirely legal
   for (C)/(D), which quantify existentially over `u°` satisfying (4)/(8). It is *not* (A)/(B):
   nothing here bears on the unforced problem. The repository README is explicit about this.

3. **The proof discharges a smaller hypothesis set than the statement offers.** The bridges use
   the PDE only at `t > 0` (`ComparatorBridge.lean:199`, `ComparatorR3Bridge.lean:30`:
   `navier_stokes : ∀ t : ℝ, 0 < t → ...`) and never use the force's spatial decay against the
   solution. Using fewer hypotheses makes the proof obligation *harder*, not easier, so this is a
   robustness margin, not a gap: the result would survive if the `t = 0` reading of (1) were
   contested.

No place was found where the Lean encoding demands *less* of the witnesses `(u₀, f)` than
Fefferman's (4)/(5)/(8)/(9), and no place where it admits *fewer* competitor solutions than
(1),(2),(3),(6),(7) / (1),(2),(3),(10),(11)+erratum. Subject to the caveats in §5, the encoding
looks faithful.

---

## 3. What the Comparator mechanism checks — and what it cannot

Config: `ComparatorChallenges/NavierStokes.json`

```json
{
  "challenge_module": "ComparatorChallenges.NavierStokes",
  "solution_module": "NavierStokes.ComparatorSolution",
  "enable_nanoda": true,
  "theorem_names": ["NavierStokes.Comparator.navier_stokes_breakdown_R3",
                    "NavierStokes.Comparator.navier_stokes_breakdown_periodic"],
  "permitted_axioms": ["propext", "Quot.sound", "Classical.choice"]
}
```

From `.lake/packages/Comparator/Comparator/Compare.lean` and `Comparator/Axioms.lean`, what is
mechanically guaranteed (given the README's assumptions — trusted challenge imports and lakefile,
working `landrun` sandbox, correct Lean kernel):

1. **Statement identity, structurally.** `compareAt` (`Compare.lean:67-108`) looks up each name in
   the *challenge* export and the *solution* export, insists both are `thmInfo` (or both
   `axiomInfo`), and requires `challengeConst != solutionConst` to fail — i.e. exact equality of
   `ConstantVal` (name, universe params, and the **type expression**). It then walks
   `challengeConst.type.getUsedConstants` and, in `Compare.loop` (`:37-59`), requires *every*
   transitively reachable constant to be `==`-equal between the two environments. So
   `NavierStokesExistenceAndSmoothnessRn`, `divergence`, `IsOnePeriodic`, `MemLp`, `Δ`, the
   structure projections, and everything they mention must be identical in both builds. A solution
   that quietly redefined `divergence` or weakened `NavierStokesExistenceAndSmoothnessRn` would be
   rejected here. This is the guarantee that makes the duplicated
   `ComparatorDefinitions.lean` safe.
2. **Axiom whitelist.** `checkAxioms` (`Axioms.lean:52-70`) walks the transitive constant closure of
   the *solution's proof terms* and throws on any `axiomInfo` outside `permitted_axioms`. `sorryAx`
   is not on the list, so any `sorry` anywhere in the closure is caught.
3. **Kernel acceptance.** The solution export is re-typechecked by Lean's kernel, and — because
   `enable_nanoda: true` — additionally by the independent `nanoda_bin` kernel.
4. **Sandboxing.** The solution is elaborated under `landrun` so it cannot tamper with the
   challenge or the checking environment.

What Comparator **cannot** check, and what therefore remains human work:

- **Whether the Lean statement means what Fefferman wrote.** Comparator compares the solution to
  the challenge; it has no opinion on whether the challenge is a good formalization. Everything in
  §2 above is outside its reach.
- **Whether the challenge file is the upstream formal-conjectures statement.** The header at
  `ComparatorChallenges/NavierStokes.lean:20-33` admits the file "has been modified": imports,
  metadata attributes, namespace, and notation adapted, utilities inlined. Nothing in the repo pins
  or verifies the upstream text. (See §5 and §6.)
- **Mathlib's definitions.** `Δ`, `gradient`, `MemLp`, `iteratedFDerivWithin`, and the `volume`
  instance on `EuclideanSpace` are taken on trust as encoding the intended mathematics.
- **Junk-value semantics.** A statement can be structurally identical in both files and still be
  vacuous. §2 is the check Comparator does not perform.
- **The `.lake` cache.** The Comparator README notes that `lake exe cache get` is acceptable only
  if you trust the cache; a tampered Mathlib `.olean` is outside the guarantee.

---

## 4. Axioms and other trust-affecting declarations

`ComparatorSolution.lean:31-32` emits `#print axioms` for both theorems, but the build was still
running during this audit, so **I did not observe the actual output** and will not assert it. What
I can assert from the source:

- **No `axiom`, `opaque`, `unsafe`, or `partial` declaration exists anywhere** under
  `NavierStokes/`, `Euler/`, `ComparatorChallenges/`, or the root `.lean` files.
  (`grep -rnE "(^|[^a-zA-Z_.])(axiom|opaque|unsafe|partial)[[:space:]]+[a-zA-Z_]"` returns only
  prose hits inside docstrings: "partial derivatives", "not an axiom", etc.)
- **No `native_decide`** anywhere (`grep -rnw native_decide` — zero hits, including Mathlib-free
  project files).
- **No `implemented_by`, `@[extern]`, `Lean.ofReduceBool`, or `Lean.trustCompiler`** anywhere.
- **No `sorry` anywhere under `NavierStokes/` or `Euler/`** (`grep -rnw sorry` — zero hits across
  643 files / ~404k lines under `NavierStokes/`). The only `sorry`s in the repository are the four
  intended challenge placeholders: `ComparatorChallenges/NavierStokes.lean:277` and `:284`, and
  `ComparatorChallenges/Euler.lean:88` and `:184`.
- **No `set_option` anywhere under `NavierStokes/`** — so no elaboration-affecting options
  (`maxHeartbeats`, `checkBinderAnnotations`, `debug.skipKernelTC`, …) are in play.
- `Classical.choice` is used explicitly and legitimately in several places
  (e.g. `NavierStokes/FinalSlowBase.lean:634`, `NavierStokes/MatchingDebtBounds.lean:463`); it is
  on the permitted list.

Consequently, **if** the build succeeds without errors, the axiom set of the two theorems can only
be a subset of `{propext, Classical.choice, Quot.sound}` — the three declared in
`ComparatorChallenges/NavierStokes.json` and in `formalization.yaml:56-59, 65-68`. Comparator
enforces this independently. This should still be confirmed against the build log, not assumed.

One operational note: `lakefile.toml` lists `ComparatorChallenges` in `defaultTargets`, so a plain
`lake build` will report `declaration uses 'sorry'` warnings for the four reference theorems. That
is expected and does not indicate a sorry in the submission; the `Euler` library is the only one
with `warningAsError = true`, and it contains no sorries.

---

## 5. Plain-language verdict

**If the Lean build succeeds and Comparator passes, this establishes Fefferman's alternatives (C)
and (D) as literally stated — with three qualifications, none of which is a hole in the encoding
itself.**

The formal statement is a careful, honest encoding. The three things that could have made it
hollow — a junk-value `fderiv` making divergence-freeness vacuous, a junk-value Bochner integral
making the energy bound vacuous, and a missing pressure-smoothness requirement — are all closed:
smoothness fields force the honest derivatives everywhere, `MemLp` forces the honest integral, and
`pressure_smooth` is present. The decay conditions are, if anything, marginally *stronger* than
Fefferman's, which strengthens rather than weakens the theorem. The result is proved for every
`ν > 0` rather than a fixed one.

Residual gaps, in descending order of importance:

1. **The reference statement has not been verified against upstream.** The whole design rests on
   `ComparatorChallenges/NavierStokes.lean` being the google-deepmind/formal-conjectures
   formalization at commit `8bf45ed70d48b2b2a501de9c00b26bfa38c573ee`, modified only in
   inessential ways. The file itself claims this (`:24-33`) but the repo has a single squashed
   commit and no pinned copy of upstream. A line-by-line diff against upstream is the single
   highest-value remaining check — in particular whether `isOnePeriodic_pressure` and the
   `integrable` field are upstream's or this project's additions. If either were added here, the
   corresponding theorem would be weaker than the upstream challenge (though `isOnePeriodic_pressure`
   would still match the Clay erratum, and `integrable` still matches Fefferman's (7)).
2. **(D) is the erratum-corrected version.** Solutions with non-periodic pressure are not excluded.
   Correct per the erratum, but worth stating plainly.
3. **These are (C) and (D), not (A) and (B).** Blowup is *forced*; the initial data is zero. The
   Millennium Prize problem asks for one of (A)–(D), and (C)/(D) do qualify, but this is not a
   construction of unforced Navier–Stokes blowup.
4. **The formalization is trusted at the Mathlib boundary.** `Δ`, `gradient`, `MemLp`, and the
   Lebesgue measure on `EuclideanSpace ℝ (Fin 3)` are assumed to mean what they say. The project
   mitigates this by re-deriving each against its own coordinate definitions
   (`ComparatorBridge.lean:28-53`), which is good practice but only pushes the trust one level.
5. **The ~404k lines of proof under `NavierStokes/` were not read.** This audit covers the
   *statement*, not the argument. That is the correct division of labour — the kernel and Comparator
   check the argument — but it means this report says nothing about whether the build succeeds.

---

## 6. Three concrete verification tasks for newcomers

### A. 10–30 minutes — read the statement against the Clay PDF, then diff the two definition copies

No Lean toolchain required.

1. Open the Clay problem PDF and put conditions (1)–(11) side by side with
   `ComparatorChallenges/NavierStokes.lean:134-269`. Check off each condition against the table in
   §2 above and try to break it: find a competitor `(v, p)` that Fefferman would count as a
   solution but the Lean structure rejects, or vice versa.
2. Run the mechanical check that the submission's definitions are the challenge's:
   ```sh
   cd lean-repo
   diff <(sed -n '34,245p' NavierStokes/ComparatorDefinitions.lean) \
        <(sed -n '61,272p' ComparatorChallenges/NavierStokes.lean)
   ```
   Expect only the `open ...` offset difference and the theorem tail. (Comparator enforces this
   semantically anyway; this is the cheap human-readable version.)
3. Confirm no sorries or exotic declarations reach the submission:
   ```sh
   grep -rnw sorry NavierStokes/ Euler/
   grep -rnE "(^|[^a-zA-Z_.])(axiom|opaque|unsafe|partial)[[:space:]]+[a-zA-Z_]" NavierStokes/
   grep -rn "import ComparatorChallenges" .
   ```
   All three should be empty.

### B. 1–3 hours — build, read `#print axioms`, and run Comparator

1. `elan` + `lake exe cache get`, then `lake build NavierStokes` (expect a long build; the repo is
   ~404k lines under `NavierStokes/` alone). Capture the log.
2. Find the `#print axioms` output for
   `NavierStokes.Comparator.navier_stokes_breakdown_R3` and `..._periodic`
   (emitted by `ComparatorSolution.lean:31-32`). Confirm it is exactly
   `[propext, Classical.choice, Quot.sound]` and in particular contains no `sorryAx`.
3. `lake build ComparatorChallenges` and confirm the *only* `declaration uses 'sorry'` warnings are
   the four reference theorems named in §4.
4. Install `landrun`, `lean4export`, and `nanoda_bin` and run
   ```sh
   lake exe comparator ComparatorChallenges/NavierStokes.json
   ```
   preferably under the `systemd-run` wrapper from the Comparator README (Linux only). Note for the
   record which of the README's five trust assumptions hold in your environment — on macOS,
   `landrun` is unavailable, so you are relying on the "pre-built `.lake` obtained without
   compromising the environment" escape hatch instead of the sandbox.
5. Independently fetch google-deepmind/formal-conjectures at
   `8bf45ed70d48b2b2a501de9c00b26bfa38c573ee` and diff
   `FormalConjectures/Millenium/NavierStokes.lean` against
   `ComparatorChallenges/NavierStokes.lean`. This closes residual gap #1.

### C. Multi-session — write an independent statement and prove it equivalent

The strongest available check on the encoding, because it tests the statement rather than the proof.

1. In a fresh file (importing Mathlib and `NavierStokes.ComparatorDefinitions`, **not** the
   challenge), write your own formalization of Fefferman's conditions from the PDF, deliberately
   choosing different primitives: divergence as `∑ i, fderiv ℝ (fun y => v y i) x (e i)` instead of
   a trace; the Laplacian as an explicit `∑ i, iteratedFDeriv ℝ 2` sum; energy as
   `∫⁻ x, ‖v x t‖₊ ^ 2 ∂volume < ⊤` in `ℝ≥0∞` (which needs **no** integrability side-condition, so
   it independently settles the junk-value question); the decay conditions per multi-index
   `α : Fin 3 → ℕ` with per-α constants rather than per-order jets.
2. Prove `MyNavierStokesSolutionRn nu u₀ f v p ↔ NavierStokesExistenceAndSmoothnessRn nu u₀ f v p`
   and the analogous statement for the periodic case, and the corresponding equivalences for
   `InitialVelocityConditionDecay` / `ForceConditionDecay`. The multi-index ↔ jet-norm direction is
   the substantive one (finite-dimensional norm equivalence); the `∫⁻` ↔ `MemLp` + `∫` direction is
   the one that formally certifies the bounded-energy hypothesis is not vacuous.
3. Derive `∃ u₀ f, MyInitialConditions u₀ ∧ MyForceConditions f ∧ ¬∃ v p, MyNavierStokesSolutionRn …`
   from the project's theorem, and check `#print axioms` on your derived statement.
4. Stretch goal: state (C) with the energy condition dropped entirely and see whether the project's
   argument still closes. Per §2 it should not — the finite-energy hypothesis is genuinely used at
   `R3FiniteEnergyComparison.lean:28-33, 59-62` — and confirming that it does not is a useful
   negative control on whether the statement is doing real work.
