# OpenAI, "Finite time blowup for Navier–Stokes" — Clay compliance and concurrent-research map

Analyst report. Date of analysis: 2026-09-08.
Target document: `/Users/benjamincox/dsm:haidaa/navier-stokes.pdf` (166 pp., author line "OPENAI").
Pages read directly with the PDF reader: 1–3, 116–126, 165–166 (references).

**Epistemic key used throughout:**
- **[V]** = verified by me against a primary document I read in full or in the relevant part (the OpenAI PDF; Fefferman's official Clay statement, whose text I extracted with `pdftotext` and cross-checked against the rendered pages).
- **[V-abs]** = verified against an abstract / publisher landing page / arXiv listing, not the full paper.
- **[S]** = reported by a web search summary or a secondary news/blog source only; **not** verified against a primary document.
- **[NOT FOUND]** = I searched and could not find it. Stated explicitly rather than guessed.

Nothing in this report is cited from memory alone. Where I could not obtain an identifier (notably the Alpöge–Buckmaster arXiv numbers) I say so instead of supplying one.

---

## Part 0. What Theorem 1.1 and the relevant Section 10 results actually say

### Theorem 1.1 (p. 1) — transcribed [V]

> **Theorem 1.1.** *For every ν > 0 there exist a force f ∈ C_c^∞(ℝ³ × (0,∞); ℝ³), a compact set K ⊂ ℝ³, and smooth velocity and pressure fields u, p on ℝ³ × [0,1) satisfying*
>
> (1.1)  ∂_t u + (u·∇)u − ν Δu + ∇p = f,  ∇·u = 0,  u(·,0) = 0,
>
> *such that supp u(·,t) ∪ supp p(·,t) ⊂ K for every 0 ≤ t < 1,*
>
> sup_{0≤t<1} ‖u(t)‖_{L²(ℝ³)} < ∞,  limsup_{t↑1} ‖u(t)‖_{L^∞(ℝ³)} = ∞.
>
> *Consequently, there is no smooth solution (u,P) on ℝ³ × [0,∞) with the same force and initial datum whose kinetic energy is uniformly bounded sup_{t≥0} ½∫_{ℝ³}|u(x,t)|² dx < ∞.*
>
> "This establishes alternative (C) in the Millennium problem statement for Navier–Stokes as stated by Fefferman in [13]. Compact support also yields the corresponding construction on 𝕋³ = ℝ³/ℤ³, establishing alternative (D) in [13]; see Corollary 10.6."

Reference **[13]** in the paper's bibliography (p. 165) is [V]:

> Charles L. Fefferman, *Existence and smoothness of the Navier–Stokes equation*, Clay Mathematics Institute, Millennium Prize Problem statement, https://www.claymath.org/wp-content/uploads/2022/06/navierstokes.pdf

i.e. exactly the document analysed in Part 1. **There is no separate erratum entry in the bibliography** (22 entries, nos. 1–22) [V].

### The Section 10 machinery [V]

- **Prop. 10.1** (p. 117): produces compact 𝒦 ⊂ ℝ³ and smooth u, p on ℝ³×[0,1) with supp u(·,t) ∪ supp p(·,t) ⊂ 𝒦 for 0 ≤ t < 1, div u = 0, and u(·,t) = p(·,t) = 0 for all sufficiently small t ≥ 0 (10.2). Built by multiplying the vector potential 𝒜, azimuthal coefficient ℬ and pressure by cutoffs c = χ_x(x)χ_t(t): u = curl(c𝒜) + cℬe_θ, p = c p^loc (10.4). Taking the curl preserves incompressibility; the temporal cutoff gives zero initial datum.
- **f is defined as the residual**: f := ℛ(u,p) = ∂_t u + (u·∇)u − Δu + ∇p on 0 ≤ t < 1 (10.5), at viscosity one. So the equation holds *by construction*; the entire difficulty is making f smooth through t = 1.
- **Lemma 10.2** (p. 118): every ∂_x^α ∂_t^j f converges uniformly on ℝ³ as t ↑ 1, with limits F_j ∈ C_c^∞(ℝ³;ℝ³) supported in 𝒦 and ∂_x^α F_j(0) = 0 (10.6).
- **Lemma 10.3** (p. 120): f extends to f ∈ C_c^∞(ℝ³ × (0,∞); ℝ³) with support in 𝒦 × [0,2], via f(x,1+σ) = Σ_{j≥0} χ_0(b_jσ)(σ^j/j!)F_j(x) with rapidly increasing b_j. The proof ends with the explicit decay statement [V]:

  > "To state the decay bounds, we choose R with 𝒦 ⊂ B(0,R) and write M_{α,m} = ‖∂_x^α ∂_t^m f‖_∞. Compact support gives, for every integer k ≥ 0, |∂_x^α ∂_t^m f(x,t)| ≤ M_{α,m}(3+R)^k (1+|x|+t)^{−k}."

  and the paragraph immediately preceding it states, in the paper's own words [V]:

  > "Such a smooth compactly supported force satisfies the decay conditions in [13, (5)]; the zero datum satisfies its initial-data conditions as well."

- **Lemma 10.4 — the energy lemma** (p. 121) [V]:

  > **Lemma 10.4.** *Let ℱ(t) = ∫_0^t ‖f(s)‖_2 ds for 0 ≤ t ≤ 1. Then ℱ(1) < ∞, and*
  >
  > (10.13)  ‖u(t)‖_2² + 2∫_0^t ‖∇u(s)‖_2² ds ≤ ℱ(t)²  (0 ≤ t < 1).
  >
  > *In particular the kinetic energy is uniformly bounded and the total dissipation on [0,1) is finite.*

  Proof: on each [0,T], T < 1, the localized fields are smooth with fixed compact support, so ½ d/dt ‖u‖_2² + ‖∇u‖_2² = ⟨f,u⟩ (10.14) exactly; divide by (‖u‖_2²+δ²)^{1/2}, discard dissipation, Cauchy–Schwarz, integrate from the zero datum, δ ↓ 0 ⇒ ‖u(t)‖_2 ≤ ℱ(t); then integrate (10.14) ⇒ ‖u(t)‖_2² + 2∫‖∇u‖_2² ≤ 2∫_0^t ℱ′ℱ = ℱ(t)². Monotone convergence gives finite integrated dissipation on [0,1).

  Note the *role*: 10.13 bounds the **constructed** local solution (supplying `sup_{t<1}‖u(t)‖_{L²} < ∞` in Theorem 1.1). It is not a bound on hypothetical competitors — those are handled by (7)/bounded-energy hypothesis in Lemma 10.5.

- **Lemma 10.5 — the comparison / weak-strong uniqueness lemma** (p. 121) [V]:

  > **Lemma 10.5.** *Fix T < 1. If v, P is a smooth solution of (1.1) at viscosity one on ℝ³ × [0,T], with the force f of Lemma 10.3, zero initial velocity, and v ∈ L^∞([0,T]; L²(ℝ³)), then v = u on that interval, where u is the localized velocity of Proposition 10.1.*

  The proof is the technically load-bearing part, and it is set up precisely to match Fefferman's hypotheses. Quoting the framing sentence [V]:

  > "The comparison argument uses smoothness and a uniform spatial L² bound on [0,T], where T < 1. These hypotheses leave the growth of spatial derivatives at infinity unrestricted. We recover the pressure gradient from the equation before passing to the limit in the localized energy identity."

  and, inside the proof, "The pressure term requires separate control because no spatial growth condition has been imposed on P." The argument: set w = v − u, π = P − p, g_ij = v_i v_j − u_i u_j; π_* := Σ R_i R_j g_ij (10.16) with Riesz multipliers; show ∇π = ∇π_* in the space–time interior by proving H_a := ∫_0^T a(∇π − ∇π_*) dt ∈ H^{−3} has ΔH_a = 0, so its Fourier transform is a weighted L² function supported at {0}, hence zero. Then a pressure-flux estimate on expanding balls (10.17)–(10.19), a difference-energy inequality ½E_R′ + ½A_R² ≤ ‖∇u‖_∞ E_R + C_T/R, and Grönwall give E_R(t) ≤ C_T′/R → 0 as R → ∞, so w ≡ 0 on [0,T].

- **Proof of Theorem 1.1** (p. 124) [V]: along x_τ = (√(2X_in τ),0,0), t = 1−τ (10.20), the inner growth asymptotic gives u_θ(x_τ,1−τ) = τ^{−A}(e_0 + O(τ^{2h})) → +∞ (10.21). "Suppose that a global smooth solution v, P with the same data had uniformly bounded kinetic energy. For every T < 1, Lemma 10.5 applies on the closed interval [0,T]. It follows that v = u throughout [0,1). Equation (10.21) then contradicts the boundedness of the smooth v on a compact neighborhood of (0,1). The force is nonzero, since (10.13) would otherwise give u = 0." General ν by the exact rescaling u_ν(x,t) = √ν u(x/√ν, t), p_ν = ν p(x/√ν,t), f_ν = √ν f(x/√ν,t) (10.22) — time unchanged, support 𝒦_ν = √ν 𝒦, ∂_x^α∂_t^m f_ν = ν^{(1−|α|)/2}(∂_y^α∂_t^m f)(x/√ν,t), energy/dissipation scaling (10.23), and a bounded-energy competitor at viscosity ν is converted back to a viscosity-one competitor v(y,t) = ν^{−1/2}v_ν(√ν y,t), which "has already been excluded."

- **Corollary 10.6 — the periodic statement** (p. 125) [V]:

  > **Corollary 10.6.** *For every ν > 0, there is a smooth force on 𝕋³ × [0,∞), compactly supported in time, for which the solution (U, P_per) of the periodic Navier–Stokes equations with viscosity ν and zero initial velocity is smooth on [0,1) and satisfies*
  >
  > limsup_{t↑1} ‖U(t)‖_{L^∞(𝕋³)} = ∞.
  >
  > *Its velocity and pressure have support in a fixed compact subset of the interior of Q_0 at every time before one. There is no global smooth periodic solution for the same datum and force.*

  Construction: rescale by λ > 1 so that λ^{−1}𝒦_ν ⋐ Q_0 = (−½,½)³, set t_0 = 1−λ^{−2}, ũ(x,t) = λu(λx, λ²(t−t_0)), p̃ = λ²p(...), f̃ = λ³f(...), extend by zero for t < t_0, then periodize by summing integer translates. Distinct translates have disjoint supports with a positive gap, so the nonlinearity is preserved translate-by-translate and ∂_tU + (U·∇)U − νΔU + ∇P_per = F_per, div U = 0, U(·,0) = 0 holds exactly.
  Uniqueness (p. 126): "on each [0,T] with T < 1, two smooth periodic solutions with the same force and datum have difference w obeying ½ d/dt‖w‖²_{L²(𝕋³)} + ν‖∇w‖²_{L²(𝕋³)} ≤ ‖∇U‖_∞‖w‖²_{L²(𝕋³)}. **Periodic integration removes both the transport and pressure terms.** Grönwall gives equality of the solutions through T."
  Decay of the periodic force: "The force is smooth, with time support contained in [t_0, 1+λ^{−2}]. On the torus this implies every decay estimate ‖∂_x^α ∂_t^m F_per(t)‖_∞ ≤ C_{α,m,N}(1+t)^{−N}."
  And the paper's own remark on the erratum: "The pressure is periodic as well, **as required in the erratum to the problem statement [13]**."

---

## Part 1. Clay compliance

### 1.1 Fefferman's official statement — exact text [V]

Source: https://www.claymath.org/wp-content/uploads/2022/06/navierstokes.pdf ("EXISTENCE AND SMOOTHNESS OF THE NAVIER–STOKES EQUATION", Charles L. Fefferman, 5 pp., the file linked as "Official Problem Description" from the Clay Navier–Stokes page). Text extracted with `pdftotext -layout` and cross-read against the rendered pages. Verbatim (mathematical notation rendered in plain text; the extraction is faithful):

**Growth/decay restrictions on the data (whole-space case):**

> Hence, we will restrict attention to forces f and initial conditions u° that satisfy
>
> (4)  |∂_x^α u°(x)| ≤ C_{αK}(1 + |x|)^{−K}  on ℝⁿ, for any α and K
>
> and
>
> (5)  |∂_x^α ∂_t^m f(x,t)| ≤ C_{αmK}(1 + |x| + t)^{−K}  on ℝⁿ × [0,∞), for any α, m, K.

**Admissibility ("physically reasonable") conditions, whole space:**

> We accept a solution of (1), (2), (3) as physically reasonable only if it satisfies
>
> (6)  p, u ∈ C^∞(ℝⁿ × [0,∞))
>
> and
>
> (7)  ∫_{ℝⁿ} |u(x,t)|² dx < C  for all t ≥ 0  (bounded energy).

**Periodic case:**

> Alternatively, to rule out problems at infinity, we may look for spatially periodic solutions of (1), (2), (3). Thus, we assume that u°(x), f(x,t) satisfy
>
> (8)  u°(x + e_j) = u°(x),  f(x + e_j, t) = f(x, t)  for 1 ≤ j ≤ n
>
> (e_j = j-th unit vector in ℝⁿ).
> In place of (4) and (5), we assume that u° is smooth and that
>
> (9)  |∂_x^α ∂_t^m f(x,t)| ≤ C_{αmK}(1 + |t|)^{−K}  on ℝ³ × [0,∞), for any α, m, K.
>
> We then accept a solution of (1), (2), (3) as physically reasonable if it satisfies
>
> (10)  u(x,t) = u(x + e_j, t)  on ℝ³ × [0,∞) for 1 ≤ j ≤ n
>
> and
>
> (11)  p, u ∈ C^∞(ℝⁿ × [0,∞)).

**The four alternatives, verbatim:**

> A fundamental problem in analysis is to decide whether such smooth, physically reasonable solutions exist for the Navier–Stokes equations. To give reasonable leeway to solvers while retaining the heart of the problem, we ask for a proof of one of the following four statements.
>
> **(A) Existence and smoothness of Navier–Stokes solutions on ℝ³.** Take ν > 0 and n = 3. Let u°(x) be any smooth, divergence-free vector field satisfying (4). Take f(x,t) to be identically zero. Then there exist smooth functions p(x,t), u_i(x,t) on ℝ³ × [0,∞) that satisfy (1), (2), (3), (6), (7).
>
> **(B) Existence and smoothness of Navier–Stokes solutions in ℝ³/ℤ³.** Take ν > 0 and n = 3. Let u°(x) be any smooth, divergence-free vector field satisfying (8); we take f(x,t) to be identically zero. Then there exist smooth functions p(x,t), u_i(x,t) on ℝ³ × [0,∞) that satisfy (1), (2), (3), (10), (11).
>
> **(C) Breakdown of Navier–Stokes solutions on ℝ³.** Take ν > 0 and n = 3. Then there exist a smooth, divergence-free vector field u°(x) on ℝ³ and a smooth f(x,t) on ℝ³ × [0,∞), satisfying (4), (5), for which there exist no solutions (p,u) of (1), (2), (3), (6), (7) on ℝ³ × [0,∞).
>
> **(D) Breakdown of Navier–Stokes Solutions on ℝ³/ℤ³.** Take ν > 0 and n = 3. Then there exist a smooth, divergence-free vector field u°(x) on ℝ³ and a smooth f(x,t) on ℝ³ × [0,∞), satisfying (8), (9), for which there exist no solutions (p,u) of (1), (2), (3), (10), (11) on ℝ³ × [0,∞).

**The nearest thing to a definition of "breakdown" (pp. 2–3), verbatim** [V]:

> For a given initial u°(x), the maximum allowable T is called the "blowup time." Either (A) and (B) hold, or else there is a smooth, divergence-free u°(x) for which (1), (2), (3) have a solution with a finite blowup time. For the Navier–Stokes equations (ν > 0), if there is a solution with a finite blowup time T, then the velocity (u_i(x,t))_{1≤i≤3} becomes unbounded near the blowup time.

### 1.2 Does Theorem 1.1 satisfy (C)? — item by item

**(a) Does a compactly supported smooth force satisfy the decay condition (5)? — YES, unconditionally.** [V]

(5) demands |∂_x^α∂_t^m f(x,t)| ≤ C_{αmK}(1+|x|+t)^{−K} on ℝ³×[0,∞) for *every* α, m, K. For f ∈ C_c^∞ with supp f ⊂ 𝒦 × [0,2] and 𝒦 ⊂ B(0,R): outside the support the left side is 0; inside, 1+|x|+t ≤ 1+R+2 = 3+R, so (1+|x|+t)^{−K} ≥ (3+R)^{−K} and the bound holds with C_{αmK} = ‖∂_x^α∂_t^m f‖_∞(3+R)^K. This is exactly the inequality the paper writes down at the end of the proof of Lemma 10.3. **Verified as both correct and explicitly claimed in the paper.** No subtlety here.

**(b) Is zero initial datum admissible under (C)? — YES.** [V]

This is a structural asymmetry between (A)/(B) and (C)/(D) that is easy to misread. (A) and (B) quantify **universally** over the datum ("Let u°(x) be **any** smooth, divergence-free vector field satisfying (4)") and **fix f ≡ 0**. (C) and (D) quantify **existentially** over both ("**there exist** a smooth, divergence-free vector field u°(x) … and a smooth f(x,t) … satisfying (4), (5)"). u° ≡ 0 is smooth, divergence-free, and satisfies (4) with C_{αK} = 0. So starting from rest is fully admissible for (C), and (C) places no lower bound on the datum, no nondegeneracy requirement, and no requirement that the force be small, decaying in a stronger sense, or physically motivated.

Two further admissibility points, both clean: (C) requires f smooth on ℝ³ × [0,∞), and f ∈ C_c^∞(ℝ³×(0,∞)) is smooth there and vanishes identically near t = 0 (matching (10.2)); and the force must be nonzero for the claim to be non-vacuous, which the paper notes explicitly ("The force is nonzero, since (10.13) would otherwise give u = 0").

**(c) Does "no smooth solution with uniformly bounded energy" match the required negation? — YES, it is verbatim the negation, and the proof's hypotheses are correctly aligned with (6)+(7).** [V]

The negation (C) demands is: *no* pair (p,u) satisfies (1),(2),(3) **together with** (6) [p,u ∈ C^∞(ℝ³×[0,∞))] **and** (7) [∫|u(x,t)|²dx < C for all t ≥ 0]. Theorem 1.1's consequence excludes exactly "smooth solution (u,P) on ℝ³×[0,∞) with the same force and initial datum whose kinetic energy is uniformly bounded". These coincide. Three points that could have gone wrong and do not:

1. **(6) requires *p* smooth too, and imposes no decay or growth condition on p at spatial infinity; (7) controls only u in L².** A weak-strong-uniqueness argument that quietly assumed decaying pressure, or π ∈ L², would *not* establish (C): a competitor could in principle carry an arbitrarily growing harmonic pressure. The paper identifies this and closes it: Lemma 10.5's hypotheses are exactly "smooth … and v ∈ L^∞([0,T];L²)", it states outright that "no spatial growth condition has been imposed on P" and "the growth of spatial derivatives at infinity [is] unrestricted", and it proves ∇π = ∇π_* by the Liouville argument on H_a (ΔH_a = 0, Fourier transform a weighted L² function supported at the origin, hence zero). This is the correct alignment with (6) rather than with a convenient stronger hypothesis. **I verified the logical alignment of hypotheses, not the correctness of the estimates (10.17)–(10.19) in detail.**
2. **Excluding only competitors that are smooth *and* have bounded energy is precisely what (C) asks** — it is a weaker exclusion than "no smooth solution whatsoever", but (C) does not ask for the stronger one. So the weakness is not a gap relative to (C).
3. **Lemma 10.5 is only local in time (T < 1) and only at viscosity one.** Both are handled: the T < 1 restriction suffices because a competitor equal to u on every [0,T], T < 1, inherits (10.21)'s divergence; and the viscosity is handled by the exact rescaling (10.22), including the *inverse* rescaling of a putative competitor at viscosity ν back to viscosity one. (C) says "Take ν > 0", i.e. for each fixed positive ν, and the paper delivers for every ν — slightly more than required.

One further alignment worth recording: (C) requires nonexistence of solutions **on ℝ³ × [0,∞)**, and the paper's local solution lives only on [0,1). This is not a mismatch: (C) is a nonexistence statement about global solutions, and the local object is only the tool. The paper additionally notes that H³(ℝ³) ↪ L^∞(ℝ³) excludes classical H³ continuation through t = 1 and that the maximal classical existence interval is [0,1) — extra information, not needed for (C).

**Verdict on (C):** as *stated*, Theorem 1.1 is a formally exact instance of alternative (C), with the admissibility conditions (4) and (5) satisfied and with the exclusion hypotheses matched to (6)+(7) rather than to stronger convenience hypotheses. I verified the *statement-level* compliance. I did **not** verify the 166-page proof, nor the Lean formalization, and I found no evidence of independent expert verification (see §1.5).

### 1.3 Is Corollary 10.6 consistent with (D)? — Mostly yes, with one precise gap and one unverifiable citation

**(a) The datum and (8):** u° ≡ 0 satisfies (8) trivially. ✓ [V]

**(b) The force and (9) — no problem, and note the printed (9) is *time*-decay only.** [V] This is the point where a careless reading of the Clay statement produces a spurious objection. (5) (whole space) decays in (1+|x|+t); **(9) (periodic) decays only in (1+|t|)**. Had (9) inherited the |x| factor, then (8)+(9) together would force f ≡ 0 and (D) would be vacuous. It does not: the printed (9) is `|∂_x^α ∂_t^m f(x,t)| ≤ C_{αmK}(1+|t|)^{−K}`. Corollary 10.6's force is compactly supported in time (support in [t_0, 1+λ^{−2}]), and the paper's own sentence — "On the torus this implies every decay estimate ‖∂_x^α∂_t^m F_per(t)‖_∞ ≤ C_{α,m,N}(1+t)^{−N}" — is exactly (9). ✓ **Verified against the exact text of (9).**

**(c) The gap: the pressure-periodicity mismatch.** [V — this is the sharpest finding in Part 1]

(D) asks for nonexistence of solutions (p,u) of (1),(2),(3),**(10)**,**(11)**. As printed:
- **(10) requires only u to be periodic** — `u(x,t) = u(x+e_j,t)`. It says nothing about p.
- **(11) requires only p,u ∈ C^∞(ℝⁿ×[0,∞))** — smoothness, no periodicity of p.

Corollary 10.6 instead concludes "There is no global smooth **periodic** solution for the same datum and force", and its uniqueness proof is a torus energy identity whose decisive step is "**Periodic integration removes both the transport and pressure terms**". That step requires the *difference of pressures* π = P − P_per to be periodic (so that ∫_{𝕋³} ∇π · w = 0). So what Corollary 10.6 actually excludes is: competitors with **periodic u and periodic p**. A hypothetical competitor with periodic u but non-periodic smooth pressure satisfies (10) and (11) as literally printed and is **not** excluded by the argument as written.

The paper is aware of this and addresses it by appeal to an erratum: "The pressure is periodic as well, **as required in the erratum to the problem statement [13]**."

**Status of that appeal — I could not substantiate it.** [V for the negative findings; NOT FOUND for the erratum]
- The document at the URL the paper cites as [13] — the current official statement, https://www.claymath.org/wp-content/uploads/2022/06/navierstokes.pdf — is 5 pages long, I read all 5 and extracted the full text, and it contains **no erratum, no addendum, no correction note, and no requirement that the pressure be periodic**. (10) and (11) read exactly as quoted above.
- The Clay Navier–Stokes problem page (https://www.claymath.org/millennium/navier-stokes-equation/) links only the "Official Problem Description" (that same 2022 PDF) and the prize rules. **It mentions no erratum.**
- The paper's bibliography contains no separate erratum entry.
- Web searches for an official Fefferman/CMI erratum on the periodic pressure turned up **no primary source**. The only material I found is a third-party remark [S] on a non-CMI site (navier-stokes.org) asserting that the uniqueness proof requires periodicity of p, that periodicity is not required in the problem statement, and that the author of that site informed Fefferman and CMI and "never received any answer" — which, if anything, is evidence *against* the existence of a published erratum.

**Therefore, stated exactly:** Corollary 10.6 establishes (D) **with (11) strengthened by the additional requirement that p be periodic** (equivalently, with a mean-zero-pressure-gradient normalisation on the torus). Under that strengthening the corollary is a correct instance of (D). Against (D) as literally printed at the cited URL, the nonexistence claim is not established by the argument as written, and the paper's justification points to a document I cannot find. Two mitigations, so as not to overstate this:
- The requirement is universally regarded as the intended reading — the periodic Navier–Stokes problem is ill-posed for uniqueness without it (a non-periodic pressure admits spatially-linear additions and the associated uniform drift family), so this is very likely an editorial defect in the printed statement rather than a mathematical gap in the paper. But *this judgement is mine, not something I verified from a CMI source.*
- The *construction* side of (D) is unaffected: the periodized U, P_per, F_per are genuinely periodic (locally finite sums of translates of compactly supported pieces), the equation holds exactly (disjoint supports keep (U·∇)U translate-diagonal), U(·,0)=0, U is smooth on [0,1), and U_θ(x̃_τ, t̃_τ) = λ√ν τ^{−A}(e_0 + O(τ^{2h})) → +∞ gives limsup_{t↑1}‖U(t)‖_{L^∞(𝕋³)} = ∞.

**(d) A minor formal point on (D):** (D) is phrased with u° and f on **ℝ³** subject to periodicity (8), whereas Corollary 10.6 is phrased on **𝕋³**. These are interchangeable under the standard identification, given a fixed pressure convention. Also, (D)/(10)/(11) impose **no** energy condition (finite energy is automatic on the torus), so the periodic negation required is simply "no global smooth (periodic-pressure) solution" — a cleaner match than in the whole-space case, where the bounded-energy hypothesis (7) is what makes the comparison lemma applicable at all.

### 1.4 Is there a gap between "no global smooth bounded-energy solution" and "breakdown" as Fefferman defines it?

**Formally: no gap.** [V] Fefferman never gives a formal definition of "breakdown". "Breakdown of Navier–Stokes solutions on ℝ³" is only the *title* of (C); the mathematical content of (C) *is* the nonexistence of a solution satisfying (1),(2),(3),(6),(7) on ℝ³×[0,∞). Theorem 1.1's "consequently" clause is that nonexistence, word for word. There is no additional obligation hidden in the word "breakdown" — no requirement to exhibit a blowing-up local solution, to identify a blowup time, to prove unboundedness of the velocity, or to characterise the singular set. (Theorem 1.1 in fact supplies more than (C) needs: an explicit local solution on [0,1) with `sup_{t<1}‖u‖_{L²} < ∞` and `limsup_{t↑1}‖u‖_{L^∞} = ∞`, i.e. a genuine realisation of Fefferman's informal picture that "the velocity … becomes unbounded near the blowup time".)

**Substantively: there are three real gaps, none of which is a defect in the claim but all of which bound what has been shown.** [V for 1–2; V for 3]

1. **The force is essential and does work on the fluid.** Fefferman's informal discussion of a "blowup time" ("Either (A) and (B) hold, or else there is a smooth, divergence-free u°(x) for which (1),(2),(3) have a solution with a finite blowup time") is stated in the context of the *unforced* problem — (A) and (B) fix f ≡ 0. (C) and (D) admit a force, and Fefferman is explicit that the four-way menu exists "to give reasonable leeway to solvers while retaining the heart of the problem." So a proof of (C) is, on the face of the official statement, a *formally sufficient* answer to the Millennium Problem; but it leaves (A) — whether smooth decaying data with **no** force can lose smoothness — entirely open, and (10.13) shows the force injects energy (ℱ(1) > 0; indeed ℱ(1) = 0 would force u ≡ 0). In the OpenAI construction the force is not a perturbation of an autonomously blowing-up flow: f is *defined* as the momentum residual of a prescribed flow (10.5), and the whole technical content is arranging enough cancellation that this residual extends smoothly through t = 1. The blowup is engineered by the force rather than generated by the equation's own dynamics from given data. (Contrast Córdoba–Martínez-Zoroa–Zheng, §2.3 below, who emphasise that placing the force in a *local well-posedness class* is what makes the blowup attributable to the dynamics.)
2. **Nothing is claimed about ν → 0 uniformity or about a universal mechanism.** The result is "for every fixed ν > 0" via an exact scaling of a single viscosity-one construction (10.22), not a ν-uniform or ν-robust statement.
3. **Nothing is claimed about the *class* of data that break down.** (C) is existential; this is one engineered example.

**On prize eligibility, separately from mathematical content** [S / not verified]: the Clay page links prize rules at https://www.claymath.org/millennium-problems/rules/. **I did not fetch or read the rules**, so I make no claim about their contents beyond noting that eligibility is a separate question from whether (C) is satisfied, and that at least one news outlet framed the result as "solving the wrong problem" (see §2.9).

### 1.5 Summary table — Clay compliance

| Requirement | Status | Basis |
|---|---|---|
| (C): u° smooth, div-free, satisfies (4) | ✓ satisfied by u° ≡ 0 | [V] (C) is existential in u° |
| (C): f smooth on ℝ³×[0,∞), satisfies (5) | ✓ satisfied by f ∈ C_c^∞(ℝ³×(0,∞)) | [V] explicit bound at end of Lemma 10.3 proof |
| (C): no (p,u) satisfying (1),(2),(3),(6),(7) | ✓ exactly the "consequently" clause of Thm 1.1 | [V] statement-level; proof not verified |
| — hypotheses aligned with (6) (pressure growth unrestricted) | ✓ handled explicitly in Lemma 10.5 | [V] |
| — every ν > 0, incl. inverse rescaling of competitors | ✓ | [V] (10.22)–(10.23) |
| "Breakdown" needs more than nonexistence? | ✗ no — (C)'s content *is* the nonexistence | [V] Fefferman gives no formal definition |
| (D): u° satisfies (8) | ✓ u° ≡ 0 | [V] |
| (D): F_per satisfies (9) (note: (9) is *time*-decay only) | ✓ | [V] exact text of (9) |
| (D): no (p,u) satisfying (1),(2),(3),(10),(11) | **⚠ partial** — established only for competitors with **periodic pressure**; (10)/(11) as printed do not require it | [V] Cor. 10.6's Grönwall step needs periodic π |
| Paper's justification: "the erratum to the problem statement [13]" | **⚠ unsubstantiated** — no erratum at the cited URL, none on the Clay page, none in the bibliography, none found by search | [V] negatives; [NOT FOUND] erratum |
| Independent verification of the 166-page proof / Lean artifact | **not established** as of 2026-09-08 | [S] see §2.9 |

---

## Part 2. Concurrent and prior research map

Ordered roughly by closeness to the OpenAI result. For each: exact statement type (forced/unforced, which equation, force regularity, self-similar or not, computer-assisted or analytic), then how the OpenAI construction differs. The OpenAI baseline to compare against is: **smooth compactly supported force (C_c^∞ in space *and* time); full Laplacian; every viscosity ν > 0; from rest (u° ≡ 0); uniformly bounded kinetic energy; L^∞ blowup in finite time; analytic proof (with a Lean artifact claimed).**

### 2.1 Córdoba & Martínez-Zoroa 2023 — forced 3D Euler, uniform C^{1,1/2−ε} force

- **Citation:** Diego Córdoba and Luis Martínez-Zoroa, *Blow-up for the incompressible 3D-Euler equations with uniform C^{1,1/2−ε} ∩ L² force*, arXiv:2309.08495 (2023). [V — this is exactly reference **[6]** of the OpenAI paper, p. 165, which lists it with that title, that arXiv id, year 2023, and no journal.] Reported by search [S] as "to appear in Duke Mathematical Journal" — **I did not verify the journal acceptance against a publisher page.**
- **Statement type** [V-abs, from the arXiv abstract]: **forced**, **3D incompressible Euler on ℝ³** (no viscosity). Solutions constructed in **C^{3,1/2} ∩ L²** on [0,T), T < ∞ finite, with a force **uniform in C^{1,1/2−ε} ∩ L²**. The blowup is ∫_0^t ‖∇u‖ ds → ∞ as t → T, with the solution remaining smooth away from the origin. **Not self-similar** — the abstract states the construction does not use self-similar coordinates. **Analytic**, not computer-assisted. Treats solutions above the C^{1,1/3+} threshold regularity for axisymmetric solutions without swirl.
- **How OpenAI differs:** (i) Navier–Stokes with the **full Laplacian** at every ν > 0, versus ν = 0; (ii) force is **C^∞ and compactly supported**, versus merely C^{1,1/2−ε} (only one Hölder derivative — a large regularity gap, and the single most important axis of improvement); (iii) solution is **C^∞** on [0,1), versus C^{3,1/2}; (iv) blowup in **L^∞ of the velocity**, versus the integral of ∇u; (v) starts **from rest**. The OpenAI paper's own historical section (p. 2) credits this line of work with establishing "a strategy for singularity formation based on amplification across scales while controlling the regularity of the external force," and describes the mechanism: "larger-scale strain amplifies smaller-scale vorticity, with leading self-interactions and feedback on the larger scales suppressed."

### 2.2 Córdoba & Martínez-Zoroa 2024/2025 — IPM with a *smooth* source (the key regularity precedent)

- **Citation:** Diego Córdoba and Luis Martínez-Zoroa, *Finite time singularities of smooth solutions for the 2D incompressible porous media (IPM) equation with a smooth source*, arXiv:2410.22920 (submitted 30 Oct 2024; revised 12 Nov 2024 and 13 Feb 2025). [V — reference **[7]** of the OpenAI paper, listed as "arXiv:2410.22920, 2024, Revised 2025"; and [V-abs] arXiv listing.]
- **Statement type** [V-abs / S]: **forced (source term)**, **2D IPM** (an active-scalar / Darcy-law model, not Navier–Stokes). Establishes smooth **finite-energy** solutions with a **compactly supported uniformly smooth source** that develop singularities in finite time. **Analytic**, not computer-assisted. Mechanism per the OpenAI paper's own summary (p. 2) [V]: "constructed singularities from smooth initial data by successive amplification of oscillatory layers, using approximations of increasing order to keep every spatial derivative of the source uniformly bounded."
- **Why this is the closest methodological ancestor:** it is the first result in this program to get the force/source all the way up to **uniformly smooth**, which is precisely the barrier that alternative (C) requires to be crossed. The OpenAI paper says its construction "also exploits dynamical amplification, with a different role for the amplified disturbances: oscillatory pulses generate a mean momentum flux that supplies the missing force on a collapsing background vortex" (p. 3) [V].
- **How OpenAI differs:** IPM in 2D is a far weaker equation than 3D Navier–Stokes (no nonlinear vortex stretching in the same sense, no viscous dissipation to defeat); OpenAI's is the full 3D NS system with the Laplacian, and additionally has compact support in space *and* time for the force and starts from rest.
- **Related, same program, authorship not independently verified** [S]: *Finite-time singularity via multi-layer degenerate pendula for the 2D Boussinesq equation with uniform C^{1,√(4/3)−1−ε} ∩ L² force*, arXiv:2505.20988, appearing in Advances in Mathematics (ScienceDirect S0001870825003780); and *Vorticity blow-up for the 2D incompressible non-homogeneous Euler equations with uniform C^{1,√(4/3)−1−ε} force*, arXiv:2605.29866. I found these in search listings but **did not verify their author lists**.

### 2.3 Córdoba, Martínez-Zoroa & Zheng — hypodissipative Navier–Stokes (ARMA 2026)

- **Citation:** Diego Córdoba, Luis Martínez-Zoroa and Fan Zheng, *Finite time blow-up for the hypodissipative Navier–Stokes equations with a force in L¹_t C^{1,ε}_x ∩ L^∞_t L²_x*, Archive for Rational Mechanics and Analysis **250** (2026), art. 38, doi:10.1007/s00205-026-02198-0; preprint arXiv:2407.06776 (2024). [V — reference **[8]** of the OpenAI paper gives exactly this volume, article number, year and DOI; [V-abs] the Springer landing page and arXiv listing corroborate title, authors and venue.]
- **Statement type** [V-abs / S]: **forced**, **fractionally dissipative ("hypodissipative") Navier–Stokes** with dissipation |∇|^α for **α ∈ [0, α₀), where α₀ = (22 − 8√7)/9 > 0** (numerically ≈ 0.093). Solutions on ℝ³ × [0,T], T < ∞, with u ∈ C^∞ ∩ L² for 0 ≤ t < T and ∫_0^t ‖∇u‖ ds → ∞ as t → T. Force in **L¹_t([0,T]) C^{1,ε}_x ∩ L^∞_t L²_x**. **Analytic.** Billed as the first blowup result for hypodissipative Navier–Stokes **in a well-posedness class** — the significance being that, because the force lies in a class where the problem is locally well posed, "the blow-up is generated by the dynamics of the equation and not by the force itself."
- **How OpenAI differs:** this is the most direct predecessor *on the viscous side*, and the contrast is stark on two axes. (i) **Dissipation strength:** α < α₀ ≈ 0.093 is *very* weak fractional dissipation; the full Laplacian is α = 2. Getting from α ≈ 0.09 to α = 2 is not a quantitative improvement but a qualitatively different regime, because at α = 2 dissipation is strong enough that Ladyženskaja/Prodi–Serrin-type criticality considerations bite and viscous energy loss can overwhelm cross-scale amplification (this is the objection Stan Palasek is reported [S] to have raised about unforced approaches — §2.9). (ii) **Force regularity:** L¹_t C^{1,ε}_x is one Hölder derivative in space and only L¹ in time, versus C_c^∞ in space *and* time. (iii) OpenAI is at **every** ν > 0, from rest, with uniformly bounded kinetic energy (10.13).

### 2.4 Albritton, Brué & Colombo 2022 — non-uniqueness of forced Leray solutions

- **Citation:** Dallas Albritton, Elia Brué and Maria Colombo, *Non-uniqueness of Leray solutions of the forced Navier–Stokes equations*, Annals of Mathematics **196** (2022), no. 1, 415–455, doi:10.4007/annals.2022.196.1.3. [V — reference **[1]** of the OpenAI paper, matching volume/issue/pages/DOI; [V-abs] Project Euclid landing page corroborates.]
- **Statement type** [V-abs / S; and [V] for the OpenAI paper's own characterisation]: **forced**, **full 3D Navier–Stokes**, but the result is **non-uniqueness, not blowup**. Exhibits **two distinct Leray(–Hopf) solutions with zero initial velocity and the same body force**. **Self-similar**: the background solution is unstable for the Navier–Stokes dynamics in similarity variables, with similarity profile a smooth compactly supported vortex ring whose cross-section modifies Vishik's unstable 2D vortex. **Analytic.** The OpenAI paper's own description (p. 2) [V]: "Albritton, Brué, and Colombo [1] subsequently constructed distinct suitable Leray–Hopf solutions with zero initial velocity and the same force, using an unstable vortex in similarity variables. **Their force lies in L¹_t L²_x and is singular at the initial time.**"
- **How OpenAI differs:** (i) the *conclusion* is different in kind — non-uniqueness of weak solutions versus finite-time L^∞ blowup of a smooth solution and nonexistence of any global smooth bounded-energy solution; (ii) their force is only L¹_t L²_x and, critically, **singular at t = 0**, so it does not satisfy (5) and cannot be used for (C); (iii) OpenAI's force is C_c^∞ and vanishes near t = 0; (iv) OpenAI's is not a similarity-variable instability argument (though the leading profile is self-similar, the terminal-time smoothness of the force is achieved by oscillatory pulses cancelling the singular residual, plus further corrections).
- **Related background cited by OpenAI** [V, ref. **[4]**]: Tristan Buckmaster and Vlad Vicol, *Nonuniqueness of weak solutions to the Navier–Stokes equation*, Annals of Mathematics **189** (2019), no. 1, 101–144, doi:10.4007/annals.2019.189.1.3 — convex integration, non-uniqueness among finite-energy weak solutions, unforced.

### 2.5 Tao 2016 — averaged Navier–Stokes

- **Citation:** Terence Tao, *Finite time blowup for an averaged three-dimensional Navier–Stokes equation*, Journal of the American Mathematical Society **29** (2016), no. 3, 601–674, doi:10.1090/jams/838; preprint arXiv:1402.0290 (v1 3 Feb 2014, v3 1 Apr 2015). [V — reference **[22]** of the OpenAI paper, matching volume/issue/pages/DOI; [V-abs] arXiv listing corroborates dates.]
- **Statement type** [V-abs]: **unforced**, a **modified ("averaged") 3D Navier–Stokes equation** in which the bilinear operator is replaced by an averaged version that still obeys the same cancellation condition and the same energy identity. Constructs a smooth solution of that averaged equation which blows up in finite time. **Not** the true Navier–Stokes equation. **Analytic.** Its stated purpose is to formalize the **supercriticality barrier**: global regularity cannot be established by any "abstract" method using only upper-bound function-space estimates on the nonlinearity together with the energy identity.
- **How OpenAI differs:** OpenAI treats the **true** Navier–Stokes nonlinearity, not an averaged surrogate, but pays for it with a nonzero force. In a sense the two occupy complementary corners of the same trade-off: Tao keeps the datum unforced and weakens the equation; OpenAI keeps the equation exact and adds a force. Neither reaches (A).

### 2.6 Chen & Hou — 3D Euler blowup with boundary (2022–2025)

- **Citations:** Jiajie Chen and Thomas Y. Hou, *Stable nearly self-similar blowup of the 2D Boussinesq and 3D Euler equations with smooth data I: Analysis*, arXiv:2210.07191 (2022); *… II: Rigorous Numerics*, arXiv:2305.05660 (2023); published as Jiajie Chen and Thomas Y. Hou, *Singularity formation in 3D Euler equations with smooth initial data and boundary*, **PNAS 122** (2025), no. 27, e2500940122, doi:10.1073/pnas.2500940122 (open access, CC BY). [V-abs for all — arXiv listings and the PNAS landing page; PubMed 40577113.] **Not** cited in the OpenAI bibliography [V].
- **Statement type** [V-abs / S]: **unforced**, **3D axisymmetric-with-swirl incompressible Euler in a periodic cylinder — i.e. with a smooth boundary** — and **2D Boussinesq**, from **smooth initial data of finite energy**. Proves finite-time **boundary** singularity. **Nearly self-similar** (not exactly self-similar): the framework proves nonlinear stability of an approximate self-similar profile constructed numerically via the dynamic rescaling formulation. **Computer-assisted**: Part II supplies rigorous numerics for the profile; Part I supplies the analytic stability framework, using a combination of weighted L^∞ and weighted C^{1/2} energy estimates plus sharp functional inequalities exploiting kernel symmetry and optimal-transport techniques.
- **How OpenAI differs:** this is the strongest *unforced* smooth-data blowup theorem in the fluid literature, but (i) it is **Euler (ν = 0)**, not Navier–Stokes — viscosity is exactly what the Millennium Problem is about; (ii) it **requires a boundary**, so it does not apply on ℝ³ or 𝕋³ as (A)–(D) demand; (iii) it is **computer-assisted**, whereas OpenAI's is presented as analytic (with a separate Lean artifact); (iv) OpenAI is **forced**, which is the price of removing the boundary and adding the full Laplacian.

### 2.7 Elgindi 2021 — 3D Euler blowup for C^{1,α} velocity, no boundary

- **Citation:** Tarek M. Elgindi, *Finite-time singularity formation for C^{1,α} solutions to the incompressible Euler equations on ℝ³*, Annals of Mathematics **194** (2021), no. 3, 647–727, doi:10.4007/annals.2021.194.3.2; preprint arXiv:1904.04795 (v1 9 Apr 2019, v2 4 May 2020). [V-abs — Project Euclid landing page and arXiv listing.] **Not** cited in the OpenAI bibliography [V].
- **Statement type** [V-abs]: **unforced**, **3D incompressible Euler on ℝ³**, **no boundary**, but only for velocities of regularity **C^{1,α}** (Hölder-continuous gradient) with suitable decay at infinity — i.e. **not C^∞ data**. Self-similar. **Analytic.**
- **Related** [S, not verified against a publisher page]: T. M. Elgindi, T.-E. Ghoul and N. Masmoudi, *On the stability of self-similar blow-up for C^{1,α} solutions to the incompressible Euler equations*, Cambridge Journal of Mathematics (2021). Also [S] T. M. Elgindi, *Singularity formation in the incompressible Euler equation in finite and infinite time*, arXiv:2203.17221 (survey).
- **How OpenAI differs:** Elgindi trades away *data regularity* (C^{1,α} instead of C^∞) to keep the problem unforced and boundary-free; OpenAI trades away *unforcedness* to keep both the data and the solution C^∞ and to add the full Laplacian at every ν.

### 2.8 Hou's numerical evidence for potential 3D Euler/NS singularities

- **Citations** [V-abs]: Thomas Y. Hou, *Potentially singular behavior of the 3D Navier–Stokes equations*, Foundations of Computational Mathematics **23** (2023), 2251–2299, doi:10.1007/s10208-022-09578-4 (published online Sept 2022); preprint arXiv:2107.06509. Companion: Thomas Y. Hou, *Potential singularity of the 3D Euler equations in the interior domain*, Foundations of Computational Mathematics (2023), doi:10.1007/s10208-022-09585-5. Earlier boundary scenario: G. Luo and T. Y. Hou (2014) [S, not verified]. Also [S] *Nearly self-similar blowup of generalized axisymmetric Navier–Stokes equations*, arXiv:2405.10916.
- **Statement type** [V-abs / S]: **unforced**, **3D axisymmetric Navier–Stokes** with smooth finite-energy data; **numerical evidence only — not a proof**. Reports nearly self-similar singular scaling with maximum vorticity amplified by a factor ~10⁷, and attributes the behaviour to a potential finite-time Euler singularity established numerically in the companion paper. The interior-domain Euler scenario is distinct from the 2014 Luo–Hou boundary scenario. **Not** cited in the OpenAI bibliography [V].
- **How OpenAI differs:** this is the only item in this map that concerns *unforced* Navier–Stokes with the full Laplacian and no boundary — i.e. genuinely alternative (A)'s regime — and it is **not a theorem**. The OpenAI result is a theorem but in the forced regime. The two are complementary: Hou's numerics bear on (A); OpenAI's theorem bears on (C).

### 2.9 The 2025 Google DeepMind / PINN work on unstable self-similar singularities

- **Precursor** [V-abs]: Yongji Wang, Ching-Yao Lai, Javier Gómez-Serrano and Tristan Buckmaster, *Asymptotic self-similar blow-up profile for three-dimensional axisymmetric Euler equations using neural networks*, **Physical Review Letters 130** (2023), 244002, doi:10.1103/PhysRevLett.130.244002; preprint arXiv:2201.06780. (Covered in Quanta, "Deep Learning Poised to 'Blow Up' Famed Fluid Equations", 12 Apr 2022 [S].)
- **The 2025 result** [V-abs]: *Discovery of Unstable Singularities*, **arXiv:2509.14185** (Sept 2025). Author list per the arXiv HTML and the Princeton/Stanford announcements [V-abs/S]: Javier Gómez-Serrano, Tristan Buckmaster, Gonzalo Cao-Labora, Yao Lai, Yongji Wang, together with Google DeepMind authors (reported as Mehdi Bennani, James Martens, Sébastien Racanière, Sam Blackwell, Alex Matthews, Stanislav Nikolov). Accompanying DeepMind blog: "Discovering new solutions to century-old problems in fluid dynamics" [V-abs].
- **Statement type** [V-abs / S]: **unforced**, **three equations — the incompressible porous media (IPM) equation, 2D Boussinesq, and 3D Euler with boundary**. Systematically discovers **new families of unstable self-similar blowup solutions** using physics-informed neural networks with curated architectures/training plus a high-precision Gauss–Newton optimizer, reaching accuracies near double-float machine precision; reports a simple empirical asymptotic formula relating blowup rate to order of instability. **Explicitly a numerical/discovery result, not a proof** — the stated aim is to reach the precision required to *enable* rigorous computer-assisted proofs. Follow-up [S, authorship not verified]: *Resolving Sharp Gradients of Unstable Singularities to Machine Precision via Neural Networks*, arXiv:2511.22819.
- **How OpenAI differs:** entirely different in kind. This is machine-assisted **discovery of candidate profiles** for unforced equations (and, for Euler, still with a boundary), producing high-precision numerics rather than theorems; OpenAI produces a theorem for forced Navier–Stokes. The overlap is sociological rather than mathematical (Buckmaster and Cao-Labora appear in both this line and the September 2026 events). Note that AI is used at opposite ends: PINNs to *find* a profile, versus LLM agents to *write and formalize* a proof.

### 2.10 Public reaction to the OpenAI paper (September 2026) — **reaction exists**

Direct answer to the instruction to note explicitly if no public reaction can be found: **public reaction does exist and is extensive, dated 7–8 September 2026.** I list it with reliability flags, because the bulk of it is secondary journalism and AI-written aggregator blogs, and my attempts to retrieve two primary sources failed.

**Primary source retrieved [V]:**
- **Terence Tao**, blog post *"Finite time blowup with smooth forcing term for the incompressible porous medium, Boussinesq, and incompressible Euler equations"*, What's New, **7 September 2026** (https://terrytao.wordpress.com/2026/09/07/...). I fetched this. Content: it discusses the **Alpöge–Buckmaster** results (finite-time blowup with smooth forcing for IPM — already handled by Córdoba–Martínez-Zoroa — plus **2D Boussinesq** and **3D incompressible Euler without boundary**, all with smooth forcing terms and smooth initial data). It explains the iterative strategy in the words "one has already managed to construct a low frequency solution … and would like to perturb it to create a new solution that adds a high frequency correction," notes that "their work has also been formalized in Lean," and states that the method "has a high likelihood of also extending to Navier–Stokes as well" while that goal remains incomplete. It **does not address Clay alternatives (C)/(D) explicitly**, and the **main text does not discuss OpenAI's announcement** (which came the following day, 8 Sept) — per the fetch, the OpenAI controversy appears only in the comment section. A later edit to the post mentions separate PINN work on **unforced** Euler by **Ganeshram, Duruisseaux and Anandkumar**, with the caveat that "actually establishing its stability … remains a major challenging task." Tao's framing sentence: "the actual solving of these problems is only a proxy goal for the primary goal of developing mathematical understanding."

**Primary sources I attempted and could NOT retrieve:**
- **Tristan Buckmaster's Mastodon post**, https://mastodon.social/@tristanbuckmaster/117233413705701198 — the fetch returned only the opening fragment "Today, Levent Alpöge and I have made public three…". **I could not read the rest**, so I cannot confirm the preprint titles, arXiv identifiers, or any allegation from the primary source.
- **arXiv identifiers for the three Alpöge–Buckmaster preprints: [NOT FOUND].** I searched arXiv listings, Tao's post, and two news write-ups; none gave the identifiers, and `arxiv.org/a/buckmaster_t_1` returns 404. I am deliberately not supplying numbers. What is corroborated across sources [S]: three preprints, posted 7–8 September 2026, by **Levent Alpöge** (Anthropic) and **Tristan Buckmaster** (NYU Courant), on **2D IPM on the torus with a uniformly smooth space-time force**, **inviscid 2D Boussinesq with smooth forcing in both equations**, and **3D incompressible Euler with smooth forcing (no boundary)**, with Lean formalizations at **https://github.com/tristanbuckmaster/fluid_lean**, and with disclosed heavy use of LLMs (reported as a mix of Claude, OpenAI Codex, and "Astra" models) over roughly a year. **Note that these results are for Euler/Boussinesq/IPM, not Navier–Stokes** — so they are not competing claims to (C)/(D), and Tao's assessment is that extending them to NS is plausible but not done.
- **Yahoo/Tech**, "OpenAI's 10,000-Agent Navier-Stokes Claim Solves the Wrong Problem — and the Right One Has a Provenance Controversy" — **HTTP 403, could not read**. The headline itself is the clearest public articulation of the (C)-versus-(A) point made in §1.4, but I cannot report its argument.

**Secondary sources read [S — journalism, treat quotes as unverified against primaries]:**
- **Scientific American**, "OpenAI claims blockbuster math breakthrough amid swirl of controversy" (~8 Sept 2026). Quotes **Diego Córdoba**: "We're a little bit in shock" and "If it's done, that will be a big surprise for us." Quotes **Luis Silvestre** (U. Chicago): "Yesterday and today are crazy days… We're all, in the community, discussing the implications of this." Reports **Sébastien Bubeck** (OpenAI) denying copying allegations. Reports the proof was "certified using the programming language Lean, which all but guarantees its correctness." Notably, per my fetch, the article **does not explicitly state whether the result resolves the Millennium Prize problem** or address the (C)/(D) versus (A)/(B) distinction.
- **kingy.ai**, "OpenAI's Navier–Stokes Proof Claim: Evidence and Dispute" — an analysis blog (likely AI-assisted; low reliability). Useful for hard artifacts it names: OpenAI announcement https://x.com/OpenAI/status/2097374646148481532; paper https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf; Lean repo https://github.com/openai/NavierStokesAndEuler; and it reports the bibliography was **expanded from 16 to 22 entries** in a revision to give Córdoba–Martínez-Zoroa fuller credit — **which is consistent with what I verified directly: the copy at hand has exactly 22 references, nos. 1–22, including [6], [7], [8] for the Córdoba line** [V]. It states the verification position as "OpenAI's claimed proof is public; independent acceptance is not established in this report," reports named technical objections from **Stan Palasek** (viscous energy loss overwhelming growth mechanisms in certain *unforced* approaches) and **Gonzalo Cao-Labora** (requesting derivative bounds for the separate Euler candidate, noting small residuals can coexist with larger derivative errors), and reports the Bubeck/Altman/Buckmaster exchange over authorship and timing (Buckmaster raising training-data timing from 28 August onward).
- **Fortune** (8 Sept 2026), **Unite.AI**, **XenoSpectrum**, **Time News**, **explainx.ai**, **glitchwire.com** — aggregator coverage of the same events, all [S], adding: an "88-hour autonomous run with 10,000 agents" figure attributed to OpenAI; Buckmaster's statement "I have not seen the proof" and "I am not accusing anyone of anything. I am stating what I was told"; and Bubeck characterising the allegations as "false and inflammatory."
- **MathOverflow: [NOT FOUND].** My searches surfaced no MathOverflow thread on the OpenAI paper. **Hacker News: [NOT FOUND]** as a direct hit; a **Lobsters** thread on Tao's post exists at https://lobste.rs/s/ki3ylq/finite_time_blowup_with_smooth_forcing [S].

**Bottom line on reaction:** as of 2026-09-08 there is heavy press coverage and an active priority dispute, one substantive positive expert assessment on record (Tao, but about **Alpöge–Buckmaster's Euler/Boussinesq/IPM results, not about the OpenAI Navier–Stokes paper**), and **no independent expert verification of the OpenAI 166-page proof or its Lean artifact that I could find**. Claims that Lean certification "all but guarantees correctness" are journalistic; verifying that a Lean development actually formalizes the *stated* Theorem 1.1 (rather than a weakened variant), and that its axioms/`sorry`s are clean, is itself a nontrivial check that I have not performed and found no record of anyone else performing.

---

## Summary of findings

1. **(C) is satisfied as stated, and the statement-level fit is careful rather than accidental.** [V] Compact smooth support gives (5) for free (the paper proves it); u° ≡ 0 is admissible because (C) is *existential* in u° (unlike (A), which is universal in u° with f ≡ 0); and the "consequently" clause is verbatim the negation "(1),(2),(3),(6),(7) has no solution". The comparison Lemma 10.5 is deliberately hypothesised on smoothness + a uniform L² bound only — matching (6)+(7) — with the unrestricted-pressure-growth difficulty handled by an explicit Liouville argument. I verified statement-level compliance only; I did not verify the proof.
2. **The one precise gap is in (D), not (C).** [V] (10) and (11) as printed at the URL the paper cites require only *u* periodic and *p,u* smooth; Corollary 10.6's Grönwall argument works only against competitors with **periodic pressure** ("Periodic integration removes both the transport and pressure terms"). So Corollary 10.6 establishes (D) with (11) strengthened by periodicity of p.
3. **The paper's justification for that strengthening — "the erratum to the problem statement [13]" — is unsubstantiated by the cited source.** [V negatives / NOT FOUND] The 5-page official statement at that exact URL contains no erratum and does not require periodic pressure; the Clay Navier–Stokes page links no erratum; there is no separate erratum entry in the 22-item bibliography; and I found no primary CMI/Fefferman erratum by search. The requirement is almost certainly the intended reading (the periodic problem is not uniqueness-well-posed without it), but that is my judgement, not a verified fact.
4. **There is no formal gap between "no global smooth bounded-energy solution" and "breakdown".** [V] Fefferman never defines "breakdown"; it is only the title of (C)/(D), whose content *is* the nonexistence. The substantive limitation is different and real: the force is essential and does work (ℱ(1) > 0 by 10.13), f is *defined* as the residual of a prescribed flow, and (A) — unforced smooth decaying data — remains untouched. Fefferman's own framing ("to give reasonable leeway to solvers while retaining the heart of the problem") makes (C) a formally sufficient answer, which is exactly why the public framing "solves the wrong problem" is about mathematical significance rather than about compliance.
5. **The literature position is coherent and the deltas are large.** The closest viscous predecessor (Córdoba–Martínez-Zoroa–Zheng, ARMA 2026) needs dissipation |∇|^α with α < (22−8√7)/9 ≈ 0.093 versus α = 2, and a force only in L¹_t C^{1,ε}_x ∩ L^∞_t L²_x versus C_c^∞. The closest smooth-force predecessor (Córdoba–Martínez-Zoroa, arXiv:2410.22920) is 2D IPM, not 3D NS. The closest unforced smooth-data theorems (Chen–Hou PNAS 2025; Elgindi Annals 2021) are Euler and pay with a boundary or with C^{1,α} data respectively. Albritton–Brué–Colombo (Annals 2022) is full NS but proves non-uniqueness with a force singular at t = 0. Tao (JAMS 2016) is unforced but for an averaged equation. Hou's FoCM 2022/23 and the DeepMind/PINN line (arXiv:2509.14185) are numerical, not proofs.
6. **Public reaction exists (7–8 Sept 2026) but independent verification does not.** Tao's 7 Sept post praises **Alpöge–Buckmaster** (Euler/Boussinesq/IPM with smooth forcing, Lean-formalized) and says extension to Navier–Stokes looks feasible; it is not an endorsement of the OpenAI paper. I could not obtain arXiv identifiers for the Alpöge–Buckmaster preprints and have not invented any. No MathOverflow thread found. No record of independent verification of the OpenAI proof or its Lean artifact.

---

## Sources

Primary documents I read directly:
- OpenAI, *Finite time blowup for Navier–Stokes*, local file `/Users/benjamincox/dsm:haidaa/navier-stokes.pdf`, pp. 1–3, 116–126, 165–166. Public copy reported at https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf
- Charles L. Fefferman, *Existence and Smoothness of the Navier–Stokes Equation*, Clay Mathematics Institute official problem statement (5 pp.), full text read and extracted: https://www.claymath.org/wp-content/uploads/2022/06/navierstokes.pdf
- Clay Mathematics Institute, Navier–Stokes problem page (checked for errata; none found): https://www.claymath.org/millennium/navier-stokes-equation/
- Clay Millennium Prize rules (link noted, **not read**): https://www.claymath.org/millennium-problems/rules/
- Terence Tao, *Finite time blowup with smooth forcing term for the incompressible porous medium, Boussinesq, and incompressible Euler equations*, What's New, 7 Sept 2026: https://terrytao.wordpress.com/2026/09/07/finite-time-blowup-with-smooth-forcing-term-for-the-incompressible-porous-medium-boussinesq-and-incompressible-euler-equations/

Literature (arXiv / publisher pages, abstract-level verification):
- Córdoba & Martínez-Zoroa, forced 3D Euler, C^{1,1/2−ε} ∩ L² force: https://arxiv.org/abs/2309.08495
- Córdoba & Martínez-Zoroa, 2D IPM with a smooth source: https://arxiv.org/abs/2410.22920
- Córdoba, Martínez-Zoroa & Zheng, hypodissipative NS, ARMA 250 (2026) 38: https://arxiv.org/abs/2407.06776 and https://link.springer.com/article/10.1007/s00205-026-02198-0
- Albritton, Brué & Colombo, Annals 196 (2022) 415–455: https://projecteuclid.org/journals/annals-of-mathematics/volume-196/issue-1/Non-uniqueness-of-Leray-solutions-of-the-forced-Navier-Stokes/10.4007/annals.2022.196.1.3.short
- Tao, averaged 3D NS, JAMS 29 (2016) 601–674: https://arxiv.org/abs/1402.0290
- Chen & Hou, Part I (Analysis): https://arxiv.org/abs/2210.07191 ; Part II (Rigorous Numerics): https://arxiv.org/abs/2305.05660 ; PNAS 122 (2025) e2500940122: https://www.pnas.org/doi/10.1073/pnas.2500940122
- Elgindi, Annals 194 (2021) 647–727: https://arxiv.org/abs/1904.04795 and https://projecteuclid.org/journals/annals-of-mathematics/volume-194/issue-3/Finite-time-singularity-formation-for-C1alpha-solutions-to-the-incompressible/10.4007/annals.2021.194.3.2.short
- Hou, potentially singular 3D NS, FoCM 23 (2023) 2251–2299: https://arxiv.org/abs/2107.06509 and https://link.springer.com/article/10.1007/s10208-022-09578-4
- Hou, potential singularity of 3D Euler in the interior domain, FoCM (2023): https://link.springer.com/article/10.1007/s10208-022-09585-5
- Wang, Lai, Gómez-Serrano & Buckmaster, PRL 130 (2023) 244002: https://arxiv.org/abs/2201.06780 and https://link.aps.org/doi/10.1103/PhysRevLett.130.244002
- Gómez-Serrano, Buckmaster, Cao-Labora, Lai, Wang & Google DeepMind, *Discovery of Unstable Singularities*, arXiv:2509.14185: https://arxiv.org/html/2509.14185v1 ; DeepMind blog: https://deepmind.google/discover/blog/discovering-new-solutions-to-century-old-problems-in-fluid-dynamics
- Follow-up PINN precision paper (authorship not verified): https://arxiv.org/abs/2511.22819
- Related Córdoba–Martínez-Zoroa-program preprints (authorship not verified): https://arxiv.org/abs/2505.20988 ; https://arxiv.org/abs/2605.29866

September 2026 commentary (secondary; quotes not verified against primaries):
- Alpöge–Buckmaster Lean formalizations (reported): https://github.com/tristanbuckmaster/fluid_lean
- Buckmaster Mastodon post (**only the first line retrievable**): https://mastodon.social/@tristanbuckmaster/117233413705701198
- Scientific American: https://www.scientificamerican.com/article/openai-claims-blockbuster-math-breakthrough-amid-swirl-of-controversy/
- Fortune, 8 Sept 2026: https://fortune.com/2026/09/08/openai-says-it-cracked-navier-stokes-math-grand-challenge-buckmaster-accusation-cheating-intimidation-tao-lament/
- Unite.AI: https://www.unite.ai/buckmaster-and-alpoge-post-ai-fluid-blowup-proofs-dispute-openai-contact/
- kingy.ai analysis blog: https://kingy.ai/blog/navier-stokes-ai-proof-claims-dispute/
- XenoSpectrum: https://xenospectrum.com/en/openai-navier-stokes-singularity-clay-dispute/
- Yahoo/Tech (**HTTP 403, not read**): https://tech.yahoo.com/ai/articles/openai-10-000-agent-navier-200053357.html
- OpenAI announcement and Lean repo, as reported by kingy.ai (**not independently opened**): https://x.com/OpenAI/status/2097374646148481532 ; https://github.com/openai/NavierStokesAndEuler
- Lobsters discussion of Tao's post: https://lobste.rs/s/ki3ylq/finite_time_blowup_with_smooth_forcing
