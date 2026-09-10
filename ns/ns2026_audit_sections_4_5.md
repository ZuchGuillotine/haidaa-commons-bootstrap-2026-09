# Sections 4-5 of "Finite time blowup for Navier-Stokes" (OpenAI, compiled 2026-09-08)

Reader's report on pages 24-62 of `/Users/benjamincox/dsm:haidaa/navier-stokes.pdf`.
All page numbers are printed page numbers, which coincide with PDF page indices in this
file. Math is transcribed in plain text; I preserve the paper's equation numbers.

Scope note: Sections 4 and 5 construct the *background* (leading profile plus all
positive-order axisymmetric corrections) and the *annular stress datum* that the wave
construction of Sections 6-8 is later required to realize. Neither section contains any
Navier-Stokes solution: everything here is an algebraic/ODE construction of profiles plus
bookkeeping of residuals. The actual solving of the equations happens elsewhere
(Sections 7-9 and Appendices A-C, which are outside my assignment).

---

## 1. Precise statements of the main results in range

### 1.0 Setup: similarity variables (Section 4.1, pp. 24-26)

Equation (4.1), p. 24, fixes the coordinates and exponents:

    tau = 1 - t,   A = 1/2 + h,   D = 1/2 - h,   0 < h < 1/2,
    z = q^D eta,   tau = q(1 - eta^2),   d = 1 - eta^2,
    L = 1 - 2 h eta^2,   s = r^2/2,   X = s/q.

Here `q = q(z,t) > 0` is defined implicitly by `tau = q(1 - (z q^{-D})^2)`; the text
argues (p. 25) that for `0 < h < 1/2` this determines a unique `q > |z|^{1/D}`, because
`q - z^2 q^{2h}` vanishes at that endpoint, tends to infinity, and has derivative
`1 - 2 h z^2 q^{2h-1} = L >= 1 - 2h > 0`. So `|eta| < 1` in the interior and `eta = ±1`
are one-sided limits. `D_X f = X d_X f` throughout (logarithmic radial derivative).

**Lemma 4.1 (p. 25).** For a smooth profile `f` and any real `b`,

    d_t(q^b f) = q^{b-1} T_b f,     d_z(q^b f) = q^{b-D} Z_b f,

with

    T_b f = L^{-1}( -b f + D eta f_eta + D_X f ),
    Z_b f = L^{-1}( 2 b eta f + d f_eta - 2 eta D_X f ).            (4.2)

Proof is the chain rule from `q_t = -L^{-1}`, `eta_t = D eta/(qL)`, `X_t = X/(qL)`,
`q_z = 2 eta q^{1-D}/L`, `eta_z = d/(q^D L)`, `X_z = -2 eta X/(q^D L)`, using
`L - 2D eta^2 = d`.

Ansatz (4.3), p. 25:

    u_theta^{(0)} = q^{-A} E,   u_z^{(0)} = q^{-A} U,   r u_r^{(0)} = V_0,
    p^{(0)} = q^{-2A} Pi,       E = C^{-1} sqrt(2X) phi,

with a fixed amplitude normalization `C > 1` and `phi` smooth. Axis regularity
(Definition 3.2) reduces to

    E = sqrt(2X) F,  V_0 = X v_0,  F = phi/C,  and U, v_0, Pi in C^inf([0,X_c] x [-1,1]).   (4.4)

Cartesian form is displayed in (4.5).

Profile identities (4.7), p. 26 (incompressibility and the leading centrifugal-pressure
balance):

    V_0 = (X/L)( 2 eta U - 2 D eta A_X(U) - d d_eta A_X(U) ),      Pi_X = E^2/(2X),

where `A_X(f)(X,eta) = X^{-1} int_0^X f(x,eta) dx` is the radial average with
`A_X(f)(0,eta) = f(0,eta)` (4.6). Definitions (4.8):

    H = sqrt(2X) E,  F = E/sqrt(2X),  l = D_X log H,
    W = 1 - 2 D eta A_X(U) - d d_eta A_X(U),  H_c = D eta + d U.

The two scalar "integrated inviscid" profiles `Q_s, N_s` are defined as the
smoothly-extending solutions of the radial ODEs (4.9):

    D_X Q_s + (1 + l) Q_s = S_q,      D_X N_s + N_s = S_n,
    S_q = -W l - h(1 - 2 eta U) - H_c (log E)_eta,
    S_n = -W D_X U - A(1 - 2 eta U) U - H_c U_eta - d Pi_eta + 4 A eta Pi + 2 eta D_X Pi.

with integrating factors `XH` and `X`; the integration constants `C_Q = C_N = 0` are
forced by smooth extension to `X = 0` (p. 26), giving (4.10) and the rescaled forms
`Q_s = F(X,eta)^{-1} int_0^1 s F(sX,eta) S_q(sX,eta) ds`, `N_s = int_0^1 S_n(sX,eta) ds`.

The stress profile (4.11), p. 27, in the component order `(r theta, r z)`:

    T_0 = F (p_s - s),   p_s = ( X Q_s / L , X N_s / (L E) ),
    s = (a, -b_s),   a = 1 - 2 D_X log E = 2 - 2 l,   b_s = 2 D_X U / E,

with `s` normalized by
`( d_r u_theta^{(0)} - u_theta^{(0)}/r , d_r u_z^{(0)} ) = -q^{-A-1/2} F s`.
`F p_s` is the inviscid contribution and `-F s` the radial-viscosity contribution. The
physical stress to be produced by the waves is `T = q^{-A-1/2} T_0`.

Leading tangential residuals (4.12), p. 27:

    R_j^{(0)} = R(u^{(0)}, p^{(0)}) . e_j + d_z^2 u_j^{(0)},   j in {theta, z},

i.e. the full momentum residual with *axial* viscosity removed and radial viscosity kept.

### 1.1 Proposition 4.2 (p. 27) — the stress identity

*Statement.* The residuals `R_theta^{(0)}, R_z^{(0)}` in (4.12) are minus the cylindrical
divergences `d_r + 2/r` and `d_r + 1/r` of `q^{-A-1/2} T_0`. In particular the equations
for vanishing leading residual stress are

    -2L (X phi_XX + 2 phi_X)/phi = S_q,      -2L (X U_XX + U_X) = S_n.        (4.13)

*Proof sketch given.* The angular equation in angular momentum `r u_theta^{(0)} = q^{-h} H`
becomes (4.14)

    (d_t + u_r^{(0)} d_r + u_z^{(0)} d_z)(q^{-h} H)
      = (q^{-h-1}/L){ W D_X H + H_c H_eta + h(1 - 2 eta U) H } = -(q^{-h-1}/L) H S_q,

the axial material derivative plus pressure gradient is `-q^{-A-1} S_n/L`; integrating
factors for the two radial divergences are `r^2` and `r`, so radial integration of the
negative inviscid residual produces `F X Q_s/L` and `F X N_s/(LE)`. Radial viscosity
contributes `d_r u_theta^{(0)} - u_theta^{(0)}/r = -q^{-A-1/2} F a` and
`d_r u_z^{(0)} = -q^{-A-1/2} F b_s`. Setting the stress to zero gives
`X Q_s/L = a`, `N_s/L = -2 U_X`, and substituting into (4.9) yields (4.13).

### 1.2 The five cumulative radial integrals (4.15) and Lemma 4.3

Definition (4.15), p. 28:

    M = int_0^X U dx,   I = int_0^X H dx,   J = int_0^X U H dx,
    S = int_0^X (U^2 - E^2/2) dx,   C_p = int_0^X E^2/(2x) dx,   Pi = Pi(0,eta) + C_p.

**Lemma 4.3 (p. 28).** For the solutions of (4.9) with `C_Q = C_N = 0`,

    Q_s = -W + [ (1-h) I - D eta I_eta - d J_eta + 2(h-D) eta J ] / (X H),
    N_s = -W U + [ D(M - eta M_eta) + 4 h eta S - d S_eta ] / X + 4 A eta Pi - d Pi_eta.   (4.16)

Proof: use `d_X(XW) = 1 - 2 D eta U - d U_eta` to integrate `-W D_X H`, `-W D_X U`;
combine angular terms into `-d(UH)_eta`; integrate `Pi_X = E^2/(2X)` by parts; divide by
`XH` and `X`.

### 1.3 Lemma 4.4 (pp. 28-30) — joining/perturbation in the five integrals

*Statement.* Fix `h in (0,1/2)`. For `i = 1,2` let `(U_i, E_i)` be smooth on
`(0,inf) x [-1,1]` with `E_i > 0`, and suppose the five integrals (4.15) define smooth
functions there (including one-sided eta-derivatives at `eta = ±1`). Write
`H_i = sqrt(2X) E_i`, `m_i = (M_i, I_i, J_i, S_i, C_{p,i})`. Fix one smooth
`Pi_ax(eta)` on `[-1,1]` and set `Pi_i(X,eta) = Pi_ax(eta) + C_{p,i}(X,eta)`, so
`Pi_i(0,eta) = Pi_ax(eta)`. Define `V_{0,i}` by (4.7), `Q_{s,i}, N_{s,i}` by (4.16), and
`p_{s,i}, a_i, b_{s,i}, T_{0,i}` by (4.11). Let `Delta f = f_2 - f_1`. Then:

**(i)** If for some `X_h > 0` one has `(U_1,E_1) = (U_2,E_2)` for `X >= X_h` and
`Delta m(X_h, eta) = 0` for all `-1 <= eta <= 1`, then
`Pi_i, V_{0,i}, Q_{s,i}, N_{s,i}, p_{s,i}, a_i, b_{s,i}, T_{0,i}` all agree for
`X >= X_h`, `eta in [-1,1]`.

**(ii)** Fix `0 < X_0 < X_1 < inf`, integer `k >= 0`, `e_min > 0`. On
`R = [X_0,X_1] x [-1,1]` put `||f||_j = max_{0<=r<=j} sup_R |d_eta^r f(X,eta)|`. If
`E_i >= e_min` on `R`, then

    || Delta(Q_s, N_s, p_s) ||_k <= C_k ( || Delta(U,E) ||_{k+1} + || Delta m ||_{k+1} ).   (4.17)

`C_k` depends only on `k, h, X_0, X_1, e_min`, the `C^{k+1}`-norm of `Pi_ax`, and a common
bound for `||(U_i,E_i,m_i)||_{k+1}`. **No radial derivative of a difference occurs.**

*Proof.* Part (i): the five difference integrands in (4.18) vanish for `X >= X_h`, so
`Delta m(X,eta) = Delta m(X_h,eta) = 0` identically in `eta`, hence all `eta`-derivatives
agree too; the common pressure datum gives `Delta Pi = Delta C_p = 0`; then (4.7) gives
`Delta V_0 = 0` and (4.16) gives `Delta Q_s = Delta N_s = 0`; radial derivatives agree
because the profiles agree on an interval; (4.11) gives equality of
`p_s, a, b_s, T_0`. Part (ii): rewrite `H_i, W_i, Pi_i` in terms of `E_i, M_i, C_{p,i}`,
use `Delta(H^{-1}) = -Delta H/(H_1 H_2)`, `Delta(E^{-1}) = -Delta E/(E_1 E_2)`, apply the
product rule to (4.16) and then to `p_s` in (4.11); `X >= X_0`, `L >= 1 - 2h > 0`,
`H_i >= sqrt(2 X_0) e_min` supply denominators.

*Remark recorded on p. 30:* part (ii) is what Proposition C.2 (radial modulation) uses:
`U, E`, the five integrals and their eta-derivatives stay close even though the radial
shears `a, b_s` change substantially; the estimate controls `p_s`, while the full stress
`T_0 = F(p_s - s)` also depends on those shears.

A `rho`-rescaled ("hatted") version is set up in (4.19), p. 30, with
`X = rho x`, `H = sqrt(rho) H^`, `(M,I,J,S,C_p) = (rho M^, rho^{3/2} I^, rho^{3/2} J^, rho S^, C_p^)`,
`p_s^(x,eta) = rho^{-1} p_s(rho x, eta)`; substitution gives the same formulas with `X`
replaced by `x`, so the estimate holds on a fixed `x`-rectangle with constants
independent of `rho`.

### 1.4 The cone conditions (Section 4.3) and Lemma 4.5 (pp. 30-32)

For `a > 0`, definitions (4.20), p. 30:

    t_s = -b_s/a,   v_s = a(1 + t_s^2),   P_c = p_{s,1} + t_s p_{s,2},   J_c = p_{s,2} - t_s p_{s,1}.

Definition (4.21), p. 31:

    U(P_c, J_c) = P_c + J_c^2/4 - |J_c| sqrt( (P_c - 2)/2 + J_c^2/16 ),
    relaxed cone condition:      P_c > 2  and  v_s < U(P_c, J_c);
    admissible stress cone condition:  those two inequalities together with v_s > 2.

**Lemma 4.5 (p. 31).** Let `a > 0`, `b_s in R`, `p_s = (p_{s,1}, p_{s,2}) in R^2`, and
define `t_s, v_s, P_c, J_c` by (4.20). If `v_s > 2`, the admissible stress cone condition
is equivalent to

    P_c > v_s,      (v_s - 2) J_c^2 < 2 (P_c - v_s)^2.                     (4.22)

Second assertion: for every nonempty compact `K` contained in
`{(a,b_s,w) in R^3 : a > 0}` such that for every `(a,b_s,w) in K`

    a - b_s w > 0,     2 b_s w + b_s^2/a + (a-2) w^2 < 2,

there exists `P_K > 0` such that for every `(a,b_s,w) in K` and every `p_{s,1} >= P_K`,
setting `p_{s,2} = w p_{s,1}` gives the relaxed cone condition (4.21); if also `v_s > 2`
these data satisfy the admissible stress cone condition.

*Proof.* Both conditions imply `P_c > 2`. The polynomial
`2(P_c - v)^2 - (v-2) J_c^2` in `v` has roots
`v_± = P_c + J_c^2/4 ± |J_c| sqrt((P_c-2)/2 + J_c^2/16)`; for `P_c > 2`,
`2 < v_- <= P_c`, so on `2 < v < P_c` positivity holds exactly when `v < v_-`, and
`U(P_c,J_c) = v_-`. For the second assertion set `c = 1 - b_s w/a`, `j = w + b_s/a`,
`v = a + b_s^2/a = v_s`; then `P_c = p_{s,1} c`, `J_c = p_{s,1} j`, and

    (v_s - 2)(w + b_s/a)^2 - 2(1 - b_s w/a)^2
      = (1 + b_s^2/a^2)( (a-2) w^2 + 2 b_s w + b_s^2/a - 2 ).

Hypotheses give `c > 0` and `G := 2c^2 - (v-2) j^2 > 0` on `K`; compactness gives
`c >= c_K`, `G >= gamma_K`, `v <= B_K`, `c v <= B_K`; choose
`P_K = max{ (B_K + 2)/c_K , 8 B_K/gamma_K }`. Then `P_c >= B_K + 2 > max{2, v}` and

    (2(P_c - v)^2 - (v-2) J_c^2)/p_{s,1}^2 = G - 4 c v/p_{s,1} + 2 v^2/p_{s,1}^2 >= gamma_K/2 > 0.

If `v > 2` the first assertion gives the relaxed condition; if `v <= 2` it follows from
`P_c > 2` and `U(P_c,J_c) > 2`. One threshold works throughout `K`, including points with
`v_s = 2`.

In stress coordinates, (4.23), p. 32:

    T_{0,theta} + t_s T_{0,z} = F(P_c - v_s),     T_{0,z} - t_s T_{0,theta} = F J_c,

so for `v_s > 2` the admissible condition is
`T_{0,theta} + t_s T_{0,z} > 0` and
`(v_s - 2)(T_{0,z} - t_s T_{0,theta})^2 < 2 (T_{0,theta} + t_s T_{0,z})^2`.
These are homogeneous in `T_0`, so the unit stress direction can keep a strict margin
even where the stress magnitude tends to zero at an annular edge — this is what
Theorem 4.6(iii) asserts.

### 1.5 Collars and flatness (4.24), p. 32

For `0 < X_a < X_b` and `0 < eps < log(X_b/X_a)`:

    C_a(eps) = [X_a, X_a e^{eps}] x [-1,1],   C_b(eps) = [X_b e^{-eps}, X_b] x [-1,1].   (4.24)

`g` is *flat at the radial edge* `X_e in {X_a, X_b}` if `d_X^j g(X_e, eta) = 0` for every
`j >= 0` and every `eta in [-1,1]`.
# Sections 4-5, part 2: pages 32-62

Continuation of `sections-4-5.md`, which stopped at (4.24) on page 32. This file covers
Theorem 4.6 and its proof, Lemmas 4.7-4.11, Proposition 4.10, and the whole of Section 5,
of `/Users/benjamincox/dsm:haidaa/navier-stokes.pdf` (OpenAI, compiled 2026-09-08).
Printed page numbers coincide with PDF page indices. Math is transcribed in plain text.

Throughout I mark each step:

- **[PROVED]** - a complete argument is on the page in my range.
- **[REF]** - stated in Section 4/5, proof delegated by name to Appendix A, B or C.
- **[ASSERTED]** - stated with no argument and no reference.
- **[VERIFIED]** - I redid the computation by hand (and, where noted, numerically).
- **[FLAG]** - I judge the write-up under-justified, or I could not follow it.

---

## 0. Where the mathematical content actually lives

Section 4 has six subsections. 4.1-4.3 (pp. 24-32) are covered in part 1. The rest:

| Item | Pages | Status |
|---|---|---|
| 4.4 leading-profile properties + **Theorem 4.6** | 32-34 | statement only |
| 4.5 Lemma 4.7 (moment solve) | 34-35 | **[REF]** = Lemmas A.1, A.2; short proof reproduced |
| 4.5 Lemma 4.8 (outer family) | 35-36 | **[REF]** = Prop A.4, Lemmas A.5, A.6, Prop A.7 |
| 4.5 Lemma 4.9 (backward stress, outer edge) | 36 | **[REF]** = Lemma A.8, Prop A.10 |
| 4.5 Proposition 4.10 (inner extension) | 36-38 | **[REF]** = Props B.2, B.3, B.5, B.8, Cor B.6, B.10; but pp. 37-38 contain a *substantive* proof with real content (the parameter chain, the reference integrals (4.34), the restoration window) |
| 4.5 Lemma 4.11 (periodic shear) | 38-39 | **[PROVED]** on p. 39 from Lemma C.1's outputs (C.8)-(C.10) |
| 4.6 **Proof of Theorem 4.6**, Steps 1-4 | 39-45 | **[PROVED]** from the above inputs |
| 5.1 Lemma 5.1 (inner coefficient ODE system) | 46-49 | **[PROVED]** |
| 5.2 Lemma 5.2 (radial extension + five moments) | 49-54 | **[PROVED]** |
| 5.3 Proposition 5.3 (induction, finite residual) | 55-56 | **[PROVED]** |
| 5.4 **Lemma 5.4** (divergence-preserving summation) | 56-59 | **[PROVED]** |
| 5.5 Proposition 5.5 (the realized base field) | 60-62 | **[PROVED]** |

So the honest summary of "proved on the page" versus "asserted with reference":

- **Section 4 proves** Propositions 4.2 and Lemmas 4.3, 4.4, 4.5 (part 1), Lemma 4.11,
  and Theorem 4.6 *conditional on* Lemmas 4.7, 4.8, 4.9 and Proposition 4.10.
- **Section 4 asserts with reference** the four hard construction results 4.7-4.10, all
  of which are appendix theorems. That is the bulk of the actual work.
- **Section 5 is essentially self-contained**: its five results are proved in full from
  Theorem 4.6 plus two appendix items (Lemma A.9 for the flat-endpoint factorization used
  in (5.22), and (A.12)/(A.45) for the exterior heat multiplier and pure power).
- I found **nothing in pp. 32-62 that is asserted with no argument and no reference** and
  that is also load-bearing, with the two partial exceptions listed in section 6 below
  (the derivative-order closure convention on p. 34, and the final uniformity sentence on
  p. 34).

Neither section contains a Navier-Stokes solution. Section 4 builds one *leading profile*
plus its annular stress datum; Section 5 builds the full *background* `(u_B, p_B)` as an
asymptotic series in `q^{2h}` whose momentum residual is exactly minus the divergence of a
prescribed annular stress plus a flat error. The stress is the input to Sections 6-8.

---

## 1. Theorem 4.6 (pp. 32-34): the exact statement

> **Theorem 4.6.** There exist fixed `h in (0, 1/100)`, `lambda > 0`, `C > 1`, and
> `0 < X_a < X_b`, together with profiles `E, U, Pi` for `X >= 0`, `-1 <= eta <= 1`, having
> the following properties. The quantities `F, H, V_0, Pi, Q_s, N_s, p_s, a, b_s` and `T_0`
> are defined by Equations (4.7) to (4.11); where `a > 0`, `t_s, v_s` are defined in (4.20).
> We use `A, d` from (4.1). The closed annulus below means the profile rectangle
> `[X_a, X_b] x [-1,1]` in `(X, eta)`. All assertions hold for every `eta in [-1,1]`.

**(i) Regularity, positivity, analyticity, pressure normalization (p. 33).**
For every finite `R > 0` the scalars `F = phi/C`, `U`, `Pi`, `V_0/X` are smooth on
`[0,R] x [-1,1]` with uniform bounds for every fixed mixed `X, eta`-derivative; derivatives
at the boundary `eta = ±1` are one-sided. `phi > 0` and `E = sqrt(2X) F > 0` for `X > 0`.
There is a fixed `X_an in (X_a, X_b)` such that on `[0, X_an] x [-1,1]` these scalars and
their fixed radial derivatives are analytic in `eta` on one common complex neighborhood of
`[-1,1]`. The inner collar `[X_a, X_an] x [-1,1]` has the positive directional margin in
part (iii). The pressure is normalized to vanish at radial infinity:

    Pi(X, eta) = - int_X^infinity E(x,eta)^2 / (2x) dx.                        (4.25)

**(ii) Exact balances and stress support (p. 33).**
The radial pressure balance `Pi_X = E^2/(2X)` of (4.7), with its smooth extension
`Pi_X = F^2` at `X = 0`, and the tangential residual identities of Proposition 4.2 with
physical stress `q^{-A-1/2} T_0`, hold **exactly**. The stress profile `T_0` is zero for
`0 <= X <= X_a` and for `X >= X_b`, and is nonzero at every `X_a < X < X_b`. For
`0 <= X <= X_a` the profiles solve (4.13).

**(iii) Cone margin on the closed annulus (p. 33).**
`F`, `a`, `v_s - 2` have positive lower bounds on the closed annulus `[X_a,X_b] x [-1,1]`.
The unit direction `n = T_0/|T_0|`, initially defined on the open radial interior, extends
smoothly to the closed annulus including the radial boundaries. There is a fixed
`kappa in (0,2)` such that throughout this rectangle

    n_theta + t_s n_z >= kappa,
    (v_s - 2)(n_z - t_s n_theta)^2 <= (2 - kappa)(n_theta + t_s n_z)^2.        (4.26)

For each `eta`, `n(X_a, eta)` is parallel to `(a(X_a,eta), -b_s(X_a,eta))`; at `X = X_b`,
`n(X_b,eta) = (1,0)` and `b_s(X_b,eta) = 0`.

**(iv) Flat weight and derivative bounds (p. 33).**
There is a smooth weight `zeta(X)`, positive on `X_a < X < X_b` and zero outside, every
radial derivative of which vanishes at `X_a` and `X_b`. With
`delta = min{1, log(X/X_a), log(X_b/X)}` there is `c > 0` independent of the multi-index,
and for each fixed `alpha = (alpha_1, alpha_2)` finite constants `C_alpha, m_alpha`, with
`d^alpha = d_X^{alpha_1} d_eta^{alpha_2}`, such that on `(X_a,X_b) x [-1,1]`

    |T_0| >= c zeta,        |d^alpha T_0| <= C_alpha zeta delta^{-m_alpha}.    (4.27)

**Note both bounds carry the same weight `zeta`** - the lower bound and every derivative
bound. This is exactly the shape that the audit of (7.25) (p. 82) argued is missing there.

**(v) Total moments and the heat exterior (p. 33).**
The limits `M(inf,eta), J(inf,eta), S(inf,eta)` of the cumulative integrals (4.15) exist;
the angular moment is made convergent by subtracting `H_pow` before integration. These four
convergent quantities satisfy

    M(inf,eta) = J(inf,eta) = S(inf,eta) = 0,
    int_0^inf (H - H_pow) dX = 0,      H_pow = sqrt(2X) c_inf X^{-A},          (4.28)

with `c_inf > 0` independent of `eta`. For `X >= X_b`,

    U = V_0 = 0,     E = c_inf X^{-A} Hcal(2d/X),
    Hcal(Z) = (1/Gamma(1+h)) int_0^inf e^{-v} v^h (1 + Zv)^{-h} dv.            (4.29)

The physical swirl there is the exact radial heat solution `c_inf s^{-A} Hcal(2 tau/s)`.
`Hcal` is positive and smooth for `Z >= 0` including its one-sided endpoint at zero. For
some fixed `X_v in (X_a, X_b)`, `U = V_0 = 0` throughout `[X_v, inf)`, including the outer
collar `[X_v, X_b]`.

**(vi) Two reserved intervals (pp. 33-34).**
Two fixed disjoint intervals `I_pos` and `I_mean` remain available for later positive-order
background corrections (Lemma 5.2) and averaged-flow corrections (Lemma 8.7), with
`sup I_pos < inf I_mean`. Both lie in `(X_a, X_v)`. For every `X` in either interval,

    U = 0,      E = c_patch (1 + eta^2)^{-1} X^{-1/2 - lambda},                (4.30)

with `c_patch > 0` independent of `X, eta`. The leading-profile construction leaves both
intervals unchanged.

**Closing sentence of the theorem (p. 34), which I regard as one of the two most
load-bearing sentences in the whole of Sections 4-5:**

> *All constants in the theorem depend only on the fixed profile choices and, for the
> derivative bounds, on the derivative order. They are independent of physical `q` and of
> all later dyadic bands and correction stages.*

---

## 2. What the background profile is, and what the stress datum is

**Background profile** (the object Theorem 4.6 delivers). A triple of scalar functions
`(E, U, Pi)` of `(X, eta)` on `[0,inf) x [-1,1]`, from which everything else is algebra:

- `E` is the azimuthal profile; the *physical* swirl is `u_theta^{(0)} = q^{-A} E`, and
  `E = sqrt(2X) F` with `F = phi/C` smooth and positive at the axis.
- `U` is the axial profile; `u_z^{(0)} = q^{-A} U`.
- `V_0` is the radial-inflow profile, `r u_r^{(0)} = V_0`, *determined* by `U` through the
  incompressibility identity (4.7); it is not independent data.
- `Pi` is the pressure profile, `p^{(0)} = q^{-2A} Pi`, *determined* by `E` through the
  centrifugal balance `Pi_X = E^2/(2X)` plus the normalization (4.25).

So the free data is really `(E, U)` on the half-line, one radial function of `X` for each
`eta`. In radial layers, reading outward:

1. `0 <= X <= X_a`: the **analytic axis core**. `T_0 = 0`; the profiles solve (4.13), i.e.
   the leading residual stress vanishes identically. Analytic in `eta` on `[0,X_an]`.
2. `X_a < X < X_b`: the **active annulus**. `T_0` is nonzero, and its direction lies in the
   admissible cone with margin `kappa`. This is the only region that forces waves.
3. `X >= X_b`: the **heat exterior**. `U = V_0 = 0` and `E = c_inf X^{-A} Hcal(2d/X)`, i.e.
   the swirl is an exact solution of the radial heat equation for `u_theta`, with no
   nonlinearity and no residual at all.

Within the annulus, four reserved intervals `I_1 < I_2 < I_3 < I_4` are carved out by
Lemma 4.8(ii): `I_1` is consumed by the Theorem 4.6 moment restoration (Step 3), `I_2` by
the heat compensation of Proposition A.7, `I_3 = I_pos` by Lemma 5.2, `I_4 = I_mean` by
Lemma 8.7 (p. 45, explicitly). Four intervals, four consumers - this bookkeeping closes.

**Stress datum.** From (4.11), with `s = (a, -b_s)` the radial-viscosity shear vector and
`p_s = (X Q_s/L, X N_s/(L E))` the inviscid vector,

    T_0 = F (p_s - s),      physical stress   T = q^{-A-1/2} T_0,

a **two-component** object in the component order `(r theta, r z)`. Proposition 4.2 says the
leading tangential residuals are exactly minus the cylindrical divergences
`(d_r + 2/r)` and `(d_r + 1/r)` of that stress. So "the residual is the divergence of a
stress" is an identity, not an approximation. Sections 6-8 must then realize `T` as the
torus-average of a quadratic wave covariance. Theorem 4.6(iii) says the *direction* of `T`
lies strictly inside the cone of realizable covariances; 4.6(iv) says the *magnitude*
degenerates only through the flat weight `zeta`, with every derivative degenerating at the
same rate up to a finite inverse power of `delta`.

Section 5's Proposition 5.5 upgrades this from the leading order to the full background:
the physical stress is `T_phys`, the summed `sum_n q^{2nh} T_n`, supported in
`X_a <= X(r,z,t) <= X_b`, with (5.43)

    |D^I (q^{A+1/2} T_phys - T_0)| <= C_I q^{2h} zeta delta^{-N_I} <= C_I q^h zeta delta^{-N_I}.

So the higher-order corrections to the stress are `O(q^{2h})` relative *and carry the same
`zeta` weight*. That is the precise sense in which Section 7 may treat `T_0` as the leading
covariance target and everything else as signed corrections.

---

## 3. Where the admissible cone margin comes from

This is the question the charter asks, so I trace it exactly.

**Step A - the cone in stress coordinates.** (4.23), p. 32:
`T_{0,theta} + t_s T_{0,z} = F (P_c - v_s)` and `T_{0,z} - t_s T_{0,theta} = F J_c`.
Since `F > 0`, for `v_s > 2` the admissible condition is
`T_{0,theta} + t_s T_{0,z} > 0` and `(v_s-2)(T_{0,z} - t_s T_{0,theta})^2 < 2(T_{0,theta} + t_s T_{0,z})^2`.
Both are **homogeneous of degree 2 (resp. 1) in `T_0`**, so they are conditions on the unit
direction `n = T_0/|T_0|` alone. That is why a strict margin can survive at an annular edge
where `|T_0| -> 0`.

**Step B - the margin is a compactness constant (p. 44).** Set

    A_n = n_theta + t_s n_z,     B_n = n_z - t_s n_theta.

In the radial interior, (4.23) gives `A_n = (P_c - v_s)/|p_s - s| > 0`.
At the two edges, `B_n = 0`, while `A_n = sqrt(1 + t_s^2)` at `X_a` and `A_n = 1` at `X_b`.
So `A_n > 0` on the **closed** rectangle. Define `G_n = 2 - (v_s - 2) B_n^2 / A_n^2`. By
(4.22), `G_n = 2 - (v_s-2)J_c^2/(P_c-v_s)^2 > 0` in the interior, and `G_n = 2` at both
edges. Then

    kappa = min{ 1, min_{[X_a,X_b]x[-1,1]} A_n, min_{[X_a,X_b]x[-1,1]} G_n } > 0,

which gives (4.26), with `kappa < 2`. **[PROVED]**, and **[VERIFIED]**: I checked that at
`X_a`, where `n` is parallel to `(a,-b_s)`, one gets
`A_n = (a + b_s^2/a)/sqrt(a^2+b_s^2) = v_s/(a sqrt(1+t_s^2)) = sqrt(1+t_s^2)` and
`B_n = (-b_s - t_s a)/sqrt(a^2+b_s^2) = 0`, exactly as claimed.

**Consequences worth stating plainly.**

- `kappa` is produced by a **compactness argument on a fixed compact rectangle**. It is a
  minimum of two continuous positive functions. **No explicit formula, no explicit lower
  bound, and no dependence on `h` is ever displayed.** This is the single most important
  non-explicit constant in Sections 4-5.
- The *order of choices* is nevertheless unambiguous and favourable. `kappa` is a constant
  of the profile, fixed once the profile parameters are fixed, in profile coordinates,
  before physical `q`, before any dyadic band index `l`, before any correction stage. That
  is exactly what the theorem's closing sentence asserts and what Step 4 of the proof
  re-affirms ("All choices were made in profile coordinates before physical `q` or any
  later correction stage is introduced, as required", p. 45).
- **This resolves the external hypothesis (a) flagged by the Sections 6-8 audit.** That
  audit noted that Proposition 7.5's positivity requires the cone margin to exceed
  `C S_*^{-1/2}` with `S_* = l^2`, since freezing `s` at `sigma u_*` in (7.28) is accurate
  only to `O(S_*^{-1/2})`, and observed that the paper only says "the errors in (7.28)
  preserve this strict inequality". Because `kappa` is fixed *before* `l_0` is chosen and
  `S_*^{-1/2} = 1/l -> 0` as the band index grows, `kappa > C S_*^{-1/2}` holds for all
  `l >= l_0(kappa, C)`, and folding that into the existing "sufficiently large lower band
  index `l_0`" of p. 64 is legitimate. **There is no circularity here.** What remains
  genuinely unquantified is the *size* of `kappa`, hence of `l_0`, not its existence.

**Where each piece of the margin is manufactured, upstream:**

- *interior of the annulus*: Lemma 4.11's `mu_L > 0` (a minimum of the four components of
  `Psi` over a compact set, again compactness), transported through (4.41).
- *inner edge* `X_a`: Proposition 4.10(ii) - `a(X_a,eta) > 0` and `v_s(X_a,eta) > 2 + c_ex`,
  plus the first-collar factorization (4.33) whose leading factor
  `B_a(0,eta) = F(X_a,eta) s(X_a,eta) != 0` fixes `n(X_a,eta) parallel (a,-b_s)`.
- *outer edge* `X_b`: Lemma 4.9 - `b_s = 0`, `2+h < a <= 2+2h`, `T_theta > 0`, and
  `2 - (a-2)(T_z/T_theta)^2 >= kappa_o > 0`, plus the factorization (4.32)
  `T_{0,theta} = e^{-4/y_b^2} y_b^{-3} b_theta(y_b,eta)` with `inf b_theta > 0`.

---

## 4. The proof of Theorem 4.6 (pp. 39-45), step by step

### Step 1 (pp. 39-41): choose parameters and join the profiles

Outer parameters `M_d, T_d, P_*, lambda, h, Pi_0` from Lemma 4.8, with the endpoint
property of Lemma 4.9. Proposition 4.10 is applied with the *same* `Pi_0` - this is the
device that breaks the apparent circularity between the outer construction (which fixes the
axis pressure datum) and the inner construction (which needs it). `Lambda` and `T_sh` are
fixed before the amplitude normalization `C`.

The parameter chain (p. 40), with `<<` the smallness convention of Definition 3.3 (choices
proceed right to left):

    0 < C^{-1} << T_sh^{-1} << Lambda^{-1} << sigma_* << delta_* << j_0
        << eps_m << h << lambda << P_*^{-1} << M_d^{-1} << 1,
    0 < omega_fin << t_1 << kappa_0 << C^{-1},

together with the derived scales

    T_d = e^{M_d} + 10,     P_* > e^{T_d},     X_R = 110 (C P_*)^{10},

and, from Lemma 4.8, `0 < h < min{1/100, lambda, e^{-T_d}}`. Then, and only then,
`N^{-1} << omega_fin` fixes the modulation frequency. **See section 6.1 below for what this
chain actually implies about `h`.**

The join. `X_h = X_R e^{-5}`; `(U_0, E_0) = (U_in, E_in)` for `X <= X_h` and
`(U_out, E_out)` for `X >= X_h`. Both definitions equal
`(4 eta, P_*(1+eta^2)^{-1}(X/X_R)^{1/10})` on a neighborhood of `X_h`, so the join is
smooth. **[VERIFIED]** that the five cumulative integrals then agree for all `X >= X_h`:
Proposition 4.10(iii) gives `m_0(X_h,eta) = m_out(X_h,eta)`, the five integrands agree for
`X >= X_h` (the displayed column vector of differences is identically zero there), so the
difference is constant `= 0`. Lemma 4.4(i) then transfers equality to
`V_0, Q_s, N_s, p_s, a_0, b_0, T_0` on `X >= X_h`.

The single coverage claim, on p. 40-41:

> *Let `I_1` be the first reserved correction interval from Lemma 4.8, and choose a closed
> interval `I = [X_-, X_+] subset (X_a, inf I_1)` containing every point of `(X_a,X_b)`
> where the joined pair may fail the admissible condition.*

The Appendix audit listed this as its **highest-severity** flag ("never stated as a single
claim"). **In Section 4 it *is* stated as a single claim**, and it is supported, though the
supporting sentence is not written out. I reconstruct it: admissible failures require
`v_s <= 2`. On `[X_good, e^{1/2} X_tail]` the outer pair is admissible by Lemma 4.8(iii);
on `[e^{1/2} X_tail, X_b)` it is admissible by Lemma 4.9; the joined pair equals the outer
pair for `X >= X_h` and `X_h = X_R e^{-5} < X_R < X_good`. So all failures lie in
`(X_a, X_good) subset (X_a, inf I_1)`. Near `X_a` the inner pair is admissible on the first
collar `(X_a, X_an]` by Proposition 4.10(ii), so `X_-` may be taken in `(X_a, X_an)` as the
text requires. **I judge the claim correct but the assembly left to the reader** - a
downgrade of the Appendix audit's flag from "high" to "presentational".

Concrete evidence that failures really occur: p. 38 shows that during moment restoration,
around `x = e^{-6}` (i.e. `X approx X_R e^{-6} < X_h`), one has `.7 <= a <= .9` and
`|b_s| <= .1/(1+w_*)`, hence

    v_s = a + b_s^2/a <= .9 + .01/.7 < 1 < 2,

so only the **relaxed** cone holds there. **[VERIFIED]** by direct arithmetic. Also
**[VERIFIED]**: `G = Q_s - b_s N_s/(a E) >= (6/7) Q_s`, since
`|b_s N_s/(aE)| <= [.1/(1+w_*)] w_* Q_s / .7 <= (1/7) Q_s`; and
`P_c = X Q_s/L - (b_s/a) X N_s/(LE) = (X/L) G = (X_R x/L) G`, matching the displayed
`P_c = X_R x G/L > 2 for every sufficiently large X_R`.

### Step 2 (pp. 41-42): realize the periodic shear at a finite frequency

Lemma 4.11 supplies a period-one loop `(a_L, -b_L)(X,eta,phi)` with mean `(a_0, -b_0)`,
every value of which satisfies the admissible condition with margin `mu_L`. Zero-mean
`phi`-primitives `Ascr, Bscr` are defined by

    d_phi Ascr = -(1/2)(a_L - a_0),   d_phi Bscr = (1/2) E_0 (b_L - b_0),
    int_0^1 Ascr dphi = int_0^1 Bscr dphi = 0,

extended by zero off `I`, and then, for an integer `N >= 1`,

    E_N = E_0 exp( Ascr(X,eta,N log X)/N ),    U_N = U_0 + Bscr(X,eta,N log X)/N.  (4.38)

**[VERIFIED]** the two shear identities on p. 41, which are the crux of the whole device.
Writing `D_X = X d_X` for the *full* derivative and using `D_X(N log X) = N`:

    D_X log E_N = D_X log E_0 + (D_X Ascr)/N + d_phi Ascr
                = D_X log E_0 + (D_X Ascr)/N - (1/2)(a_L - a_0),
    a_N = 1 - 2 D_X log E_N = a_0 + (a_L - a_0) - 2 (D_X Ascr)/N = a_L - 2 D_X Ascr / N.

    D_X U_N = D_X U_0 + (D_X Bscr)/N + d_phi Bscr = E_0 b_0/2 + E_0 (b_L - b_0)/2 + (D_X Bscr)/N
            = E_0 b_L / 2 + (D_X Bscr)/N,
    b_N = 2 D_X U_N / E_N = e^{-Ascr/N} ( b_L + 2 D_X Bscr /(N E_0) ).

Both reproduce the printed formulas exactly. **The mechanism is now transparent:** the
`O(N)` term produced by differentiating the fast phase is *designed* to convert the mean
shear `(a_0,-b_0)` into the loop value `(a_L,-b_L)` at leading order, and the genuine `X`
derivative only contributes at `O(1/N)`.

The estimate (4.39) `||U_N-U_0, E_N-E_0||_k <= A_k/N` and `||a_N-a_L, b_N-b_L||_k <= B_k/N`
is then immediate because `|| . ||_k` is a maximum over **eta-derivatives only** and the
phase `N log X` has no `eta` dependence. **[VERIFIED]**. This is the point of Lemma
4.4(ii)'s design: "no radial derivative of a difference occurs". The paper says it
explicitly on p. 42 ("the estimate uses no smallness of radial derivatives of `E_N - E_0`
or `U_N - U_0`"). This is genuinely elegant and, as far as I can see, correct.

The stability of the cone under the `O(1/N)` perturbation is (4.41), using the smooth map
`Psi(a,b,p) = (a, v-2, c-v, 2(c-v)^2 - (v-2)j^2)` of (4.35) and its Lipschitz constant
`C_Psi` on a fixed compact neighborhood, with `mu_L, mu_R > 0` the compact minima.
**[VERIFIED]** that positivity of all four components of `Psi(a,b_s,p_s)` is exactly
`a > 0`, `v_s > 2` and (4.22), as the text says: `v = a + b^2/a = a(1+t_s^2) = v_s`,
`c = p_1 - (b/a)p_2 = p_{s,1} + t_s p_{s,2} = P_c`, `j = p_2 + (b/a)p_1 = J_c`.

### Step 3 (pp. 42-43): solve the five moment equations

On `(Y_0,Y_1) subset I_1`, where the modulation has not touched anything, `U_N = U_0 = 0`
and `E_N = E_0 = K(eta) X^{-1/2-lambda}`. Two bumps `beta_1,beta_2` and three bumps
`gamma_1,gamma_2,gamma_3` with ordered disjoint supports give
`u_c = sum alpha_i(eta) beta_i(X)`, `e_c = sum xi_j(eta) gamma_j(X)`, and the five exact
changes at `Y_1` in the order `(M, J, I, S, C_p)` are (4.42)

    int u_c dX
    int (H_0 u_c + sqrt(2X) U_0 e_c + sqrt(2X) u_c e_c) dX
    int sqrt(2X) e_c dX
    int (u_c^2 - E_0 e_c - e_c^2/2) dX
    int (E_0 e_c / X + e_c^2/(2X)) dX          = B(eta) c + Q_eta(c,c) = d_N(eta).

**[VERIFIED]** all five rows against the definitions (4.15):
`Delta M = int u_c`; `Delta I = int sqrt(2X) e_c`;
`Delta J = int (U_f H_f - U_N H_N) = int (H_0 u_c + sqrt(2X) U_0 e_c + sqrt(2X) u_c e_c)`;
`Delta S = int (2 U_0 u_c + u_c^2 - E_0 e_c - e_c^2/2)`, and the `2 U_0 u_c` term is absent
from the display precisely because `U_0 = 0` on `(Y_0,Y_1)`;
`Delta C_p = int (2 E_0 e_c + e_c^2)/(2X) dX`. All match.

**[VERIFIED]** the block-diagonal structure and the exponent lists. Because `U_0 = 0`
there, rows 1-2 involve only `alpha` and rows 3-5 only `xi`:

    B_U rows:  (int beta_i,  int H_0 beta_i)             weights X^0, X^{-lambda}
    B_E rows:  (int sqrt(2X) gamma_j, -int E_0 gamma_j, int E_0 gamma_j / X)
                                                        weights X^{1/2}, X^{-1/2-lambda}, X^{-3/2-lambda}

using `H_0 = sqrt(2X) E_0 = sqrt 2 K X^{-lambda}` and `E_0 = K X^{-1/2-lambda}`. Exponents
`{0,-lambda}` and `{1/2,-1/2-lambda,-3/2-lambda}` are distinct exactly when `lambda > 0`,
so Lemma 4.7 (= Lemma A.1) applies to both blocks. This matches the printed exponent lists
character for character. **The decoupling into two Vandermonde-type blocks is a consequence
of `U_0 = 0` on the reserved interval, and is never said in those words**, but it is what
makes the `2x2 + 3x3` split legitimate rather than a `5x5` solve.

The quadratic smallness is handled by Lemma 4.7's second assertion: with
`||d_N||_{C^k} <= D_k/N`, choosing `N >= 8 max_{k=0,2} nu_k^2 q_k D_k` gives `8 nu_j^2 q_j
||d||_{C^j} <= 1` for `j = 0, 2`. **[VERIFIED]** arithmetic. The inverse norms `nu_k, q_k`
are finite because `lambda`, the radii and all bump shapes were fixed before `N`.

### Step 4 (pp. 43-45): verify the six conclusions

The final profiles satisfy `(U_f,E_f) = (U_0,E_0)` for `X <= X_-` and `X >= Y_1`, and
`m_f(X) = m_0(X)` for `X >= Y_1`, so by Lemma 4.4(i) the pressures, `V_0`, `Q_s`, `N_s`,
`p_s`, `s_0`, `T_0` also agree there (4.43). The six parts are then checked in order:

- (i) regularity: changes in Steps 2-3 are compactly supported above `X_an`, so the
  analytic axis rectangle survives; positivity of `E_f` from the exponential in (4.38).
- (i)/(ii) pressure: `C_p(inf,eta) = -Pi_0(eta)` gives
  `Pi(X,eta) = C_p(X,eta) - C_p(inf,eta) = -int_X^inf E^2/(2x) dx`, which is (4.25), and
  `Pi_X = E^2/(2X) = F^2` extends smoothly to the axis. Proposition 4.2 then gives the
  tangential identities. **[VERIFIED]**.
- (ii) support: `T_0 = 0` for `X <= X_a` from (4.13) at the first edge, and for `X >= X_b`
  by (4.43) plus Lemma 4.9. Nonvanishing on the open annulus from (4.23):
  `T_{0,theta} + t_s T_{0,z} = F(P_c - v_s) > 0`. **Note** this needs only `P_c > v_s`,
  which follows from the *relaxed* condition alone, since Lemma 4.5's proof shows
  `U(P_c,J_c) = v_- <= P_c` whenever `P_c > 2`. So the coverage argument closes even on the
  sub-interval where only the relaxed condition is known. **[VERIFIED]**; the paper does
  not point this out.
- (iii) the `kappa` construction, as in section 3 above.
- (iv) the weight: `zeta(X) = exp(-c_a/y_a^2 - 4/y_b^2)` on `(X_a,X_b)`, zero outside, with
  `c_a = t_1^2` from (4.33) and the `4` matching (4.32). On a fixed inner collar
  `zeta/e^{-c_a/y_a^2}` is bounded above and below; likewise `zeta/e^{-4/y_b^2}` on an outer
  collar. The nonzero factors in (4.33), (4.32) then give (4.27).
- (v) the moment identities survive both exact connections; `X_v` precedes the terminal
  collar and lies after all four reserved intervals; `V_0 = (2 eta X U - 2 D eta M -
  d M_eta)/L = 0` for `X >= X_v`. **[VERIFIED]** consistency with (4.7).
- (vi) `I_pos = I_3`, `I_mean = I_4`; `I_2` was used by the heat compensation and `I_1` by
  Step 3, so exactly the two untouched intervals remain.

**Overall judgement of Section 4.6.** Conditional on the four appendix results, this proof
is complete and I could follow every step. The only place where I had to reconstruct an
argument the paper leaves implicit is the coverage claim of Step 1, and the block
decoupling in Step 3.

---

## 5. Section 5 in detail (pp. 45-62)

### 5.0 The expansion and its bookkeeping (pp. 45-47)

Physical coefficients (5.1), `lambda_n = 2nh`:

    u_theta,n = q^{-A+lambda_n} E_n = r q^{-A-1/2+lambda_n} phi_n / C,   E_n = sqrt(2X) phi_n / C,
    u_z,n = q^{-A+lambda_n} U_n,    r u_r,n = q^{lambda_n} V_n,    p_n = q^{-2A+lambda_n} Pi_n.

Two independent reasons the spacing is `2h` and not something else, both stated on p. 46
and both **[VERIFIED]**:

1. *Axial viscosity.* `d_z(q^b f) = q^{b-D} Z_b f` by (4.2), so two axial derivatives cost
   `q^{-2D} = q^{-1} q^{2h}` (using `1 - 2D = 2A - 1 = 2h`). Two radial derivatives cost
   `q^{-1}` alone via `X = r^2/(2q)`. So axial viscosity at order `n` is fed by order
   `n - 1`. Hence the `Z^{[2]}_{b,n-1} phi_{n-1}` and `Z^{[2]}_{c,n-1} U_{n-1}` terms in
   (5.3), (5.4).
2. *Pressure.* The leading radial balance `r d_r p^{(0)} = (u_theta^{(0)})^2` has both sides
   at power `q^{-2A}`, whereas the remaining radial-momentum terms multiplied by `r` start
   at `q^{-1} = q^{-2A} q^{2h}`. Hence the order shift in (5.5) and the collection of the
   remaining radial terms into `Omega_k` in (5.6), one order later.

The order-`n` system is (5.3)-(5.6), with (5.2) the incompressibility relation
`d_X V_n = -Z_{-A,n} U_n`. The system is **linear in the order-`n` unknowns** because every
transport sum `sum_{i+j=n}` splits into `(i,j) = (0,n), (n,0)` (linear) and `1 <= i,j < n`
(known). **[VERIFIED]** for the pressure row: the paper's display
`Pi_n' = 2 C^{-2} phi_0 phi_n + C^{-2} sum_{i=1}^{n-1} phi_i phi_{n-i} - Omega_{n-1}/(2X)`
is exactly (5.5) with the `i+j=n` sum split.

### 5.1 Lemma 5.1 (pp. 47-49): the inner solve

> **Lemma 5.1.** There exists an interval `0 <= xi <= a`, `xi = sqrt X`, extending from the
> axis into the inner collar of Theorem 4.6(i) where the stress cone inequalities hold with
> a uniform positive margin, on which every positive order has a unique solution of (5.2)
> to (5.6) with `phi_n(0,eta) = U_n(0,eta) = Pi_n(0,eta) = 0`. The profiles are smooth in
> `X` and, for each fixed radial derivative, holomorphic on a neighborhood of `[-1,1]` in
> `eta`. That neighborhood and its bounds may depend on the order and the radial
> derivative; **`a` does not**.

The proof is the most technically interesting argument in my range, and it is complete.

*Step 1* recasts the system as a first-order `6 x 6` singular ODE in `xi = sqrt X` (5.7):

    K_n = A_X(U_n) - U_n,   W_n = (phi_n, U_n, K_n, Pi_n, d_xi phi_n, d_xi U_n)^T,
    d_xi W_n + xi^{-1} diag(0,0,2,0,3,1) W_n = A_0 W_n + A_1 d_eta W_n + f_n,  W_n(0) = 0.

The auxiliary `K_n` replaces the radial average `A_X(U_n)` by a first-order equation, so the
system contains no unevaluated radial integral of the current unknowns.

*Step 2* is the key. The `d_eta` coefficient `A_1` has the sparse block form displayed on
p. 48: its only possibly nonzero entries are `(A_1)_{51}, (A_1)_{52}, (A_1)_{53},
(A_1)_{62}, (A_1)_{63}, (A_1)_{64}` - i.e. **rows 5,6 and columns 1-4 only**. Hence `A_1`
maps the first four coordinates into the last two and annihilates the last two.
Consequently, for any intervening diagonal kernel `D_0`,

    A_1(xi) D_0 A_1(s) = 0,     A_1(xi) D_0 d_eta A_1(s) = 0.

**[VERIFIED]**: `(A_1 D_0 A_1)_{ij} = sum_k A_1[i,k] D_0[k,k] A_1[k,j]` needs
`k in {1,..,4}` (from the left factor) and `k in {5,6}` (from the right factor)
simultaneously - impossible. So a nonzero composition of `k` factors of
`K = G(A_0 + A_1 d_eta)` contains at most `p_k = ceil(k/2)` parameter derivatives.

Splitting the Cauchy radius loss `Delta = rho - rho'` among at most `p_k` derivatives, with
diagonal kernels of norm at most one and simplex volume `a^{k+1}/(k+1)!`, gives (5.8)

    ||K^k G f_n||_{S_rho'} <= C_n^{k+1} a^{k+1} / (k+1)! * (max{1, p_k/Delta})^{p_k}.

**[VERIFIED]** that the `k`-th root tends to zero: the factorial part contributes
`~ C a e / k`, the derivative-loss part contributes `~ (k/(2 Delta))^{1/2}`, so the product
is `~ C a e / sqrt(2 Delta k) -> 0`. The Picard series converges **for every finite `a`,
regardless of `C_n`** - which is exactly the assertion that the radial interval does not
shrink with the order. Had the nilpotency failed and `p_k = k`, the derivative-loss factor
would be `(k/Delta)^k`, whose `k`-th root is `k/Delta -> infinity`, and the product with
`Cae/k` would converge to `Cae/Delta`: convergence would then require `a` small depending on
`C_n`, i.e. a shrinking interval. **So the `ceil(k/2)` really is load-bearing, and it is
proved.** This is the best-argued single point in my range.

*Step 3* recovers smoothness at the axis using `(Gg)_i = xi int_0^1 t^{c_i} g_i(t xi) dt`
and its differentiated form, and even/odd parity to see that the profiles are smooth in
`X = xi^2` on the half-interval.

The `eta`-analyticity neighborhood `S_rho` is allowed to shrink with the order `n` and with
the radial derivative; only `a` is order-independent. That is consistent with (5.17)'s
`C_{n,m} < infinity` with constants "which may grow with `n` and `m`", and with the
diagonal cutoff device of Lemma 5.4 which is designed to absorb exactly such growth.

### 5.2 Lemma 5.2 (pp. 49-54): radial extension and the five moments

The inner profiles do not yet have compact support; radial integration of a cut-off `E_n,
U_n` leaves a nonzero exterior pressure constant and a nonzero exterior stress tail. Five
total moments are imposed:

    m_{n,1} = int_0^inf R U_n dR,   m_{n,2} = int_0^inf R^2 E_n dR,   m_{n,3} = int_0^inf d_R Pi_n dR,
    m_{n,4} = int_0^inf R^2 sum_{i+j=n} U_i E_j dR,
    m_{n,5} = int_0^inf ( R sum_{i+j=n} U_i U_j - (1/2) R^2 d_R Pi_n ) dR,      (5.10),(5.11)

with `R = sqrt(2X)`. The solve is on the reserved patch `I_pos`, where `U_0 = 0` and
`E_0 = e_* f(eta) R^{-1-2lambda}`.

**[VERIFIED]** the change of variable: `E_0 = c_patch f X^{-1/2-lambda}` from (4.30) and
`X = R^2/2` give `E_0 = c_patch f 2^{1/2+lambda} R^{-1-2lambda}`, so
`e_* = c_patch 2^{1/2+lambda} > 0`, exactly as claimed.

**[VERIFIED]** both moment matrices and their exponent lists:

    (B_U)_{ij} = int R^{p_i} b^U_j dR,   (p_1,p_2) = (1, 1-2lambda),
    (B_E)_{ij} = int R^{s_i} b^E_j dR,   (s_1,s_2,s_3) = (2, -2-2lambda, -2lambda).

Row by row. `m_{n,1}`: linear part `int R sum alpha_j b^U_j dR`, exponent `1`.
`m_{n,4}`: with `U_0 = 0` the only linear contribution is `U_n E_0`, giving
`int R^2 b^U_j e_* f R^{-1-2lambda} dR`, exponent `2-1-2lambda = 1-2lambda`.
`m_{n,2}`: `int R^2 b^E_j dR`, exponent `2`. `m_{n,3}`: using
`d_R Pi_n = 2 E_0 E_n / R + R^{-1}(sum_{i=1}^{n-1} E_i E_{n-i} - Omega_{n-1})`, the linear
part is `2 e_* f int R^{-2-2lambda} b^E_j dR`, exponent `-2-2lambda`. `m_{n,5}`: the
pressure term contributes `-(1/2) int R^2 d_R Pi_n dR`, whose linear part is
`-e_* f int R^{-2lambda} b^E_j dR` after `dX = R dR`; exponent `-2lambda`. All five match
the printed lists, and the printed normalizations of

    d_{U,n} = (m^0_{n,1}, m^0_{n,4}/(e_* f))^T,
    d_{E,n} = (m^0_{n,2}, m^0_{n,3}/(2 e_* f), -m^0_{n,5}/(e_* f))^T,

with `B_U alpha_n = -d_{U,n}`, `B_E beta_n = -d_{E,n}`, are consistent with those linear
parts including the signs. **[VERIFIED]**. Exponents are distinct for every `lambda > 0`, so
Lemma A.1 applies. `1/f = 1 + eta^2` is a polynomial, so all `eta`-derivatives of the
coefficients are bounded on `[-1,1]`, as the text says.

**Step 3** closes velocity and pressure outside the correction region: `m_{n,1} = 0` gives
`F_n(X,eta) = 0` for `X >= X_+`, hence `A_X(U_n) = V_n = 0` and `U_0 = V_0 = 0` there;
`m_{n,3} = 0` gives `Pi_n = 0` for `X >= X_+`. (5.18): `E_n = U_n = F_n = V_n = 0` on
`I_mean` for every `n >= 1` - this is what Lemma 8.7 later relies on, and it is proved here.

**Step 4** cancels the total tangential residual integrals using conservative forms. The
angular one at `n = 1` needs the compensated identity (5.19)
`int_0^inf r^2 (u_{theta,0}(r,t) - P(r)) dr = 0`. **[VERIFIED]** that this is exactly
Theorem 4.6(v)'s `int_0^inf (H - H_pow) dX = 0`: the change of variable `X = r^2/(2q)` gives
`int r^2 u_{theta,n} dr = q^{3/2-A+lambda_n} int R^2 E_n dR = q^{3/2-A+lambda_n} int H_n dX`,
so the two statements coincide with `P <-> H_pow`.

The consequence is the backward representation on p. 54,

    T_{n,theta}(R) = R^{-2} int_R^inf rho^2 r_{theta,n}(rho) drho,
    T_{n,z}(R) = R^{-1} int_R^inf rho r_{z,n}(rho) drho,

which is what proves compact support of the stress. **[VERIFIED]** that this follows from
the vanishing of the total integrals plus (5.9).

**Step 5** bounds the one exceptional term `T_1` near the outer edge, using the exact
terminal multiplier `f_o` with `1 - f_o = e^{-4/delta_b^2}` times a smooth factor, and
Lemma A.9 to gain three powers of `delta_b` from the backward integration, giving (5.22)
`|T_1| <= C e^{-4/delta_b^2} delta_b^{-3}` and `|d^I T_1| <= C_I e^{-4/delta_b^2}
delta_b^{-N_I}`. This is the only place in Section 5 where an appendix lemma is used
quantitatively. **[REF]**.

### 5.3 Proposition 5.3 (pp. 55-56): the induction

> `|F_slow(U^{[N]})|_m <= C_{N,m} q^{2h(N+1) - K_m}`, `K_m` **independent of `N`**.  (5.25)

The augmented residual `F_slow(u,p,T) = R(u,p) + (d_r + 2/r) T_theta e_theta +
(d_r + 1/r) T_z e_z` (5.24) adds back the tangential divergence of the retained stress, so
that the residual to be estimated is the *uncancelled* part only.

The `N`-independence of `K_m` is the structural point and it comes from (5.26): a
physical derivative of a profile `g(x_perp/sqrt q, eta)` costs

    |d_{x_perp}^beta d_z^k d_t^l [q^b g]| <= C_m (1+|b|)^m ||g||_{C^m} q^{b - a/2 - Dk - l},

with `a = |beta|`, so the loss `a/2 + Dk + l` depends only on the derivative multi-index and
not on which coefficient it hits. The `n`-dependence enters only through `(1+|b|)^m` with
`b = -A + lambda_n`, i.e. polynomially in `n`, and is absorbed into `C_{N,m}`. **[VERIFIED]**
the chain rule statement against `d_t = q^{-1} L^{-1}(-q d_q + X d_X + D eta d_eta)` and
`d_z = q^{-D} L^{-1}(2 eta q d_q - 2 eta X d_X + d d_eta)` of (4.2). The claim
"transverse differentiation introduces a factor `q^{-1/2}`" is the `X = r^2/(2q)` scaling.

Also stated and used: on the common stress support `X >= X_a > 0`, so `r^{-1} =
q^{-1/2}(2X)^{-1/2}` also has `N`-independent derivative losses.

### 5.4 Lemma 5.4 (pp. 56-59): the summation lemma

This is the result the Sections 9-10 audit flagged as the single most important thing it
could not check. **It is proved in full on pp. 57-59, and the concerns raised there are
answered.** I record the answer carefully because it matters.

*Hypotheses.* A domain `Omega` with `0 < q < q_0 <= 1`, `q = q(z,t)` smooth with
`|d_z^a d_t^b q| <= C_{a,b} q^{1-aD-b}` (5.28); a divergence-free base `U_0` with
`|U_0|_m <= C_m q^{-K_m} Lambda_log(q)^{P_m}`; smooth physical tuples
`Z_j = (A_j, B_j, p_j, T_j)`, `U_j = (curl A_j + B_j, p_j, T_j)` with

    |Z_j|_m <= C_{j,m} q^{g_j - l_m} Lambda_log(q)^{P_{j,m}},   0 < g_1 <= g_2 <= ... -> inf,  (5.30)

**where the losses `l_m` are independent of `j`**; coefficient bounds (5.32) with constants
independent of the truncation index; and the residual hypothesis (5.33)

    |F(U^{[J]})|_m <= C_{J,m} q^{rho_J - K^F_m} Lambda_log(q)^{P_{J,m}} + E_{J,m},
    rho_J -> infinity,   K^F_m independent of J,   0 <= E_{J,m} <= C_{J,m,N} q^N for every N.

*Conclusion.* There are numbers `a_{j+1} >= 2 a_j` with `a_1^{-1} < q_0` and a fixed cutoff
`chi` (one on `[0,1/2]`, zero on `[1,inf)`) such that

    U = U_0 + sum_{j>=1} ( curl(chi(a_j q) A_j), chi(a_j q) B_j, chi(a_j q) p_j, chi(a_j q) T_j )  (5.34)

is locally finite and smooth with `div u = 0`, and

    |U - U^{[J]}|_m <= 2^{-J} q^{g_{J+1}/2 - l'_m}   for 0 < q < 1/(2 a_J), J >= max(1,m),  (5.35)

with **`l'_m` independent of `J`**; and, under (5.33),

    for all m, N there exist delta, C > 0 with |F(U)|_m <= C q^N on 0 < q < delta.   (5.36)

**Answers to the two questions the Sections 9-10 audit posed.**

1. *"The `a_j` are never displayed."* They are displayed, not in closed form but as an
   explicit finite selection criterion, in Step 2 (p. 58):

       choose a_{j+1} >= 2 a_j, each cutoff chi(a_j q) supported in q < q_0, and
       Chat_{j,m} Lambda_log(q)^{Pbar_{j,m}} q^{g_j/2} <= 2^{-j}    (0 < q <= a_j^{-1}, 0 <= m <= j).  (5.37)

   The condition is imposed only for `m <= j` - the **diagonal device** - so at step `j`
   only finitely many requirements are imposed, each on constants `Chat_{j,m},
   Pbar_{j,m}` already known, and each holds for `q` small because `g_j > 0`. So `a_j` is
   chosen recursively and the recursion is well-founded. **This is exactly the standard
   Borel/Whitney diagonalization the 9-10 audit predicted would be needed, and it is
   written down.** I regard the earlier "never displayed" as **refuted**.

2. *"The tail must be flat to all orders simultaneously because the `a_k` were fixed before
   `N` was named."* That is precisely what (5.35) delivers, and I **[VERIFIED]** the
   summation:

       sum_{j > J} 2^{-j} q^{g_j/2 - l'_m} <= q^{g_{J+1}/2 - l'_m} sum_{j>J} 2^{-j}
                                            = 2^{-J} q^{g_{J+1}/2 - l'_m},

   using `q < 1` and `g_j` nondecreasing. For `j` with `a_j^{-1} < q` the corresponding term
   is identically zero because `chi(a_j q) = 0` there, so the restriction `0 < q <= a_j^{-1}`
   in (5.37) is not violated. Since `g_j -> infinity` and `l'_m` does not depend on `J`, the
   right-hand side of (5.35) is `O(q^N)` for any `N` once `J` is large. **The tail is flat
   to all orders simultaneously.**

3. *Quantifier order.* Step 3 states it explicitly and correctly (p. 59): *"The order of
   choices is `(m,N)`, then `J`, then the neighborhood and constants; the comparison uses
   only the remainder `E_{J,m}` at that fixed truncation."* The explicit requirement is

       g_{J+1}/2 - l'_{m+s} >= N + H_m + (d-1)(K_{m+s} + 1),      rho_J - K^F_m >= N + 1,

   with `s, d` the fixed order and degree of the differential polynomial (5.31) and `H_m`
   independent of `J`. So `J` is chosen *after* `(m,N)`, which is legitimate. **[VERIFIED]**
   that this is consistent: `g_{J+1} -> infinity` and `rho_J -> infinity` with everything
   else fixed.

4. *Is flatness conditional?* **Yes, and correctly so.** (5.36) is stated "under the
   residual hypothesis (5.33)". Lemma 5.4 is a general summation lemma with two customers:
   Proposition 5.5 (which verifies (5.33) from Proposition 5.3's (5.25), with
   `rho_N = 2h(N+1)` and *no* remainder `E`) and Proposition 9.9 (which verifies it from
   Lemmas 9.7, 9.8). That is ordinary mathematics, not a gap. The Sections 9-10 audit's
   statement that "flatness in Lemma 5.4 is conditional" is **confirmed**, but the
   implication it drew - that the headline local result is not checkable - is now only about
   whether *Section 9's* hypothesis verification is right, not about Lemma 5.4 itself.

*One genuine subtlety I want on record.* The comparison in Step 3 is with **one finite
partial sum**, not with a summed series of residuals:

    |F(U) - F(U^{[J]})|_m <= C_m q^{-H_m} |e|_{m+s} (1 + |U^{[J]}|_{m+s} + |e|_{m+s})^{d-1},
    e = U - U^{[J]}.                                                              (5.39)

This is why `K_m` and `H_m` can be `J`-independent: the "bad" factor `|U^{[J]}|_{m+s}` is
bounded by `C_{J,m} q^{-K_m} Lambda_log^{P}` with `K_m` independent of `J` (5.38) - which in
turn holds because *every increment has `g_j > 0`*, so adding terms never worsens the
negative power. That is the architectural reason the whole device works, and the paper says
it in one clause ("with `K_m >= 0` independent of `J`: every increment has `g_j > 0`").
**[VERIFIED]** and I consider it correct.

*Extension (5.40).* "At step `j` one may add any finite list of requirements
`B_{j,k} Lambda_log(q)^{P_{j,k}} q^{gamma_{j,k}}/2 <= 2^{-j}`, `gamma_{j,k} > 0`." This is
what Proposition 5.5 uses to obtain (5.46). It is legitimate for exactly the same reason as
(5.37): finitely many conditions at each step, each satisfiable by taking `a_j` larger.

### 5.5 Proposition 5.5 (pp. 60-62): the realized base field

> **Proposition 5.5.** There are smooth axisymmetric fields `(u_B, p_B)` for `q > 0` and two
> physical tangential stress components `T_phys`, supported where `X_a <= X(r,z,t) <= X_b`,
> with `div u_B = 0` and
>
>     R(u_B, p_B) = -(d_r + 2/r) T_phys,theta e_theta - (d_r + 1/r) T_phys,z e_z + E_B.   (5.41)
>
> For every fixed `0 < X_max < infinity`, uniformly on `[0,X_max] x [-1,1]`, for every
> Cartesian space-time derivative `d^alpha` and every `M > 0`,
> `|d^alpha E_B| <= C_{alpha,M} q^M` as `q -> 0`.

Plus: the normalized comparisons (5.42) `D^I(q^A u_{theta,B} - E_0) = O_I(q^{2h})` etc.; the
weighted stress comparison (5.43); vanishing of the positive-order corrections beyond `X_+`;
and (5.44), the exact agreement of the base velocity with the leading profile on `I_mean`.

The proof is four steps and is complete. Two points worth recording:

- **The divergence-free property is preserved exactly**, because the cutoffs are applied to
  the *vector potentials* `Acal_n = S_n e_theta / r` before taking the curl (5.27), not to
  the velocities. (5.45) displays the resulting radial term explicitly. Since
  `Acal_n = (1/2) q^{-A+lambda_n}(F_n/X)(-x_2,x_1,0)` and `F_n/X` is smooth at `X = 0`, the
  potential extends smoothly across the axis for `q > 0`. **[VERIFIED]** that
  `curl(S e_theta / r) = -(d_z S/r) e_r + (d_r S / r) e_z` reproduces the radial and axial
  velocities of (5.27).
- The exponent loss for the coefficient tuple is `l_m = 2A + m`, "since `D < 1/2`", which is
  the `j`-independent loss required by (5.30). **[VERIFIED]** dimensionally against (5.26).

Step 4 recovers the full expansion after truncation: for fixed `N, m` choose
`M >= max(N,m)` with `g_{M+1}/2 >= g_{N+1}`, split at `J = max{2N+2, m+1}`, and use the
geometric bound `sum_{n >= J} 2^{-n} q^{nh} <= 2^{1-J} q^{Jh} <= 2^{1-J} q^{2h(N+1)}`.
**[VERIFIED]** the arithmetic with `J >= 2N+2`.

---

## 6. Findings: gaps, suspected errors, and things I could not follow

I found no outright error in pp. 32-62. The following are, in decreasing order of how much
they would matter if they were wrong.

### 6.1 The stated range `0 < h < 1/100` is not the operative constraint [FLAG, medium]

Three different ranges for `h` appear:

- p. 4 (physical description) and p. 16 (Theorem 3.1): `0 < h < 1/100`.
- p. 24, (4.1): `0 < h < 1/2`, with the note on p. 24 "The coordinate identities in this
  subsection hold for `0 < h < 1/2`; Theorem 4.6 will choose `h < 1/100`, as in Theorem 3.1."
- p. 34, **Lemma 4.8**: `T_d = e^{M_d} + 10`, `P_* > e^{T_d}`, and
  `0 < h < min{1/100, lambda, e^{-T_d}}`.

Combining with the chain on p. 40 (`h << lambda << P_*^{-1} << M_d^{-1} << 1`, with `M_d`
large), the operative bound is

    h < e^{-T_d} = e^{-(e^{M_d} + 10)},

i.e. `h` is bounded by a *doubly exponentially small* quantity in a parameter `M_d` that is
itself only constrained to be "large". Nothing is inconsistent - Theorem 4.6 and Theorem 3.1
are existence statements about *some* `h` - and I found no circularity in the chain (reading
`<<` right-to-left: `M_d`, then `P_*`, then `lambda`, then `h`, matches Lemma 4.8's stated
order of choices exactly, and `lambda << P_*^{-1} < e^{-T_d}` makes `h < e^{-T_d}`
automatic). But three consequences deserve to be stated:

1. **No admissible parameter tuple is exhibited anywhere in the paper.** Twelve nested
   `<<` relations plus `omega_fin << t_1 << kappa_0 << C^{-1}` and a later `N^{-1} << omega_fin`
   are all governed by Definition 3.3's "sufficiently small ... depending on all data
   already fixed". Every such relation is standard practice; the concern is their number
   and the double exponential.
2. **The blow-up growth exponent `A = 1/2 + h` is therefore unquantified**, and so is the
   `L^3` divergence rate `tau^{-4h/3}` recorded in the paper overview. The construction is
   "Type II" only by an amount that cannot be named.
3. Presentationally, `0 < h < 1/100` in the physical description reads like a range one may
   choose from; it is not.

### 6.2 The cone margin `kappa` has no explicit lower bound [FLAG, medium]

As traced in section 3. Every ingredient (`mu_L` from Lemma 4.11, `min A_n`, `min G_n`,
`c_ex` from Proposition 4.10, `kappa_o` from Lemma 4.9) is a compactness minimum. Nothing
in Sections 4-5 lets a reader estimate `kappa`, and hence nothing lets a reader estimate
`l_0` in Sections 6-8. The *existence* and the *order of quantifiers* are fine (see section
3); only the size is unavailable. Since the paper's whole architecture is asymptotic, this
does not threaten correctness, but it does mean the claim "the asymptotic regime begins at
`tau` of order `1e-100`" recorded elsewhere in this Common is itself unjustified from these
sections - the true threshold could be far smaller.

### 6.3 Closure of the list of `eta`-derivative orders is a convention, not a count [FLAG]

p. 34: *"A parameter derivative in these statements means an `eta`-derivative, with all
construction constants held fixed. ... Only finitely many such norms must be small when
choosing parameters. Once those parameters are fixed, every other fixed derivative order has
a finite bound."* Proposition 4.10 correspondingly begins "Fix the outer-profile data of
Lemma 4.8 **and any prescribed finite collection of `eta`-derivative orders**."

Each use of Lemma 4.4(ii) costs one order (`(4.17)` bounds `|| . ||_k` by `|| . ||_{k+1}`).
Within Section 4's own proof the accounting is visible and closes: (4.40) uses order `k+1`
inputs for an order-`k` output, and (4.41) needs only `|| . ||_0`; p. 38 notes explicitly
that "the extra derivative required by (B.35) is included among the prescribed orders". But
the *global* count across A.7 -> A.8 -> C.2 -> C.3 is never carried out, and Section 4's
paragraph is what stands in for it. This is the Appendix audit's flag #2, and I confirm it
from the Section 4 side: **Section 4 asserts closure by convention and does not count.**

### 6.4 The uniformity sentence on p. 34 is asserted, but structurally supported [FLAG, low]

*"All constants ... are independent of physical `q` and of all later dyadic bands and
correction stages."* The supporting argument is one sentence at the end of Step 4 (p. 45):
all choices were made in profile coordinates, before `q` is introduced. That is a real
structural argument and I believe it, but it is a one-liner carrying the weight of the
repeated "constants are uniform in band, label, rectangle copy and point" assertions on
pp. 69, 73, 78, 90, 97. It is the Section 4 half of the Sections 6-8 audit's item 7, and it
**partially discharges** that item: the profile-side constants really are fixed first.

### 6.5 Two places where the paper leaves the reader to assemble the argument [FLAG, low]

- The coverage claim of Step 1 (section 4 above): correct, assembled from Lemma 4.8(iii),
  Lemma 4.9 and Proposition 4.10(ii), but the assembly is not written.
- The block-diagonal structure of (4.42): it holds *because* `U_0 = 0` on `I_1`, which is
  never said. A reader who does not notice this will not see why a `5x5` solve splits into
  `2x2` and `3x3` blocks with distinct exponent lists.

### 6.6 `-Z Hcal'/Hcal < h/4` (A.56) quantified [resolved, was Appendix flag #8]

Lemma 4.9 asserts `2 + h < a <= 2 + 2h` on the outer collar. I derived exactly

    a = 2 + 2h + 2 Z Hcal'(Z)/Hcal(Z),      Z = 2d/X,

so `a <= 2 + 2h` is automatic (since `Hcal' < 0`), and `a > 2 + h` is *equivalent* to
`-Z Hcal'/Hcal < h/2`. Numerically (artifact `ns2026_heat_shear_threshold.py`) the
threshold for the stronger asserted form `-Z Hcal'/Hcal < h/4` sits at `Z* ~ 0.42`, almost
independent of `h` (`Z* = 0.4230` at `h = 1e-4`, `0.4190` at `h = 1e-2`, `0.4036` at
`h = 0.05`), because `-Z Hcal'/Hcal = h Z (1 + O(Z))`. So the asserted inequality is the
concrete requirement `X > about 4.73 d`, a *fixed* lower bound on radius and **not** a
further smallness condition on `h`. With `X >= e^{1/2} X_tail > X_R = 110 (C P_*)^{10}` and
`P_* > e^{e^{M_d}+10}`, it holds with enormous room. This converts an "asserted in one
sentence" step into a checked one.

---

## 7. What I verified by hand, listed

Numbered so records can cite individual items.

1. `Hcal` of (4.29) satisfies exactly `Z^2 Hcal'' + ((2+2h)Z + 1) Hcal' + h(1+h) Hcal = 0`,
   which is precisely the reduction of `d_t K = (d_rr + r^{-1}d_r - r^{-2})K` for
   `K = c_inf s^{-A} Hcal(2 tau / s)`, `s = r^2/2`, `A = 1/2 + h`. Proof: substitute, group
   the three `(1+h)`-terms into `(1+h) v^h w^{-h-2}(w - Zv)^2` with `w = 1 + Zv`, note
   `w - Zv = 1`, and close with the integration by parts of `d/dv[e^{-v} v^{h+1} w^{-h-1}]`.
   Numerically confirmed to `1e-11` relative (artifact `ns2026_heat_profile_check.py`).
2. `Hcal(0) = 1` (the `Gamma(1+h)` normalization), `Hcal > 0`, `Hcal' < 0`, and
   `A Hcal + Z Hcal' = Gamma(1+h)^{-1} int e^{-v} v^h w^{-h-1}[A + (A-h) Z v] dv > 0` since
   `A - h = D = 1/2 - h > 0`; hence `K_s = -c_inf s^{-A-1}(A Hcal + Z Hcal') < 0` and
   `K_r = r K_s < 0`, which is Lemma 4.8(iii)'s monotonicity claim.
3. The identity `2 d/X = 2 tau/s` from (4.1), so (4.29)'s profile and the physical heat
   solution are the same object.
4. `a = 2 - 2 l` with `l = D_X log H`, `H = sqrt(2X) E`, is consistent with
   `a = 1 - 2 D_X log E` of (4.11); on the reference inner branch
   `E = P_* f(eta) x^{1/10}` one gets `H propto X^{3/5}`, `l = 3/5`, `a = 4/5` **exactly**,
   matching the value `a(X_i,eta) = 4/5` on p. 37 and the window `.7 <= a <= .9` on p. 38.
5. On the heat exterior `a = 2 + 2h + 2 Z Hcal'/Hcal <= 2 + 2h`, matching Lemma 4.9's upper
   endpoint exactly; the lower endpoint `2 + h` is equivalent to `-Z Hcal'/Hcal < h/2`
   (see 6.6).
6. The restoration-window arithmetic on p. 38: `v_s <= .9 + .01/.7 < 1 < 2`;
   `|b_s N_s/(aE)| <= (1/7) Q_s` hence `G >= (6/7) Q_s >= 3 Q_min/7`; and
   `P_c = (X/L) G = X_R x G/L`.
7. Lemma 4.11's proof: `int_0^{2pi} (dphi/dtheta') dtheta' = (a/v_l)(1 + t_s^2 + rho_l/a) =
   (a(1+t_s^2) + rho_l)/v_l = (v_s + rho_l)/v_l = 1`, using `<t_l> = t_s`,
   `<(t_l - t_s)^2> = rho_l/a`, `v_l = v_s + rho_l`; and
   `int_0^1 (a_L,-b_L) dphi = (a/2pi) int_0^{2pi} (1, t_l) dtheta' = (a, a t_s) = (a, -b_s)`.
   Also `-b_L/a_L = t_l` and `a_L(1 + t_l^2) = v_l` from (4.37).
8. The map `Psi` of (4.35): `v = a + b^2/a = v_s`, `c = P_c`, `j = J_c`, so positivity of
   all four components is exactly `a > 0`, `v_s > 2` and (4.22).
9. The shear identities (Step 2 of the Theorem 4.6 proof): `a_N = a_L - 2 D_X Ascr/N` and
   `b_N = e^{-Ascr/N}(b_L + 2 D_X Bscr/(N E_0))`, both derived above; and the reason (4.39)
   costs no radial derivative.
10. All five rows of (4.42) against the definitions (4.15), including the two terms absent
    because `U_0 = 0` on `(Y_0,Y_1)`; the resulting block-diagonal structure; and the
    exponent lists `{0,-lambda}` for `B_U` and `{1/2,-1/2-lambda,-3/2-lambda}` for `B_E`.
11. The smallness threshold `N >= 8 max_{k=0,2} nu_k^2 q_k D_k` against Lemma 4.7's
    hypothesis `8 nu_j^2 q_j ||d||_{C^j} <= 1`.
12. Theorem 4.6(iii)'s edge values: at `X_a`, `A_n = sqrt(1+t_s^2)` and `B_n = 0`; at `X_b`,
    `A_n = 1`, `B_n = 0`, `t_s = 0`.
13. `T_{0,theta} + t_s T_{0,z} > 0` follows from the *relaxed* condition alone, because
    `U(P_c,J_c) = v_- <= P_c` when `P_c > 2`; so the "nonzero at every interior point"
    conclusion of part (ii) does not need admissibility everywhere.
14. Lemma 5.1's nilpotency `A_1 D_0 A_1 = 0` from the sparse block form (rows 5-6, columns
    1-4 only), hence `p_k = ceil(k/2)`; and that `(5.8)^{1/k} ~ C a e/sqrt(2 Delta k) -> 0`,
    so the Picard series converges for *every* finite `a`. With `p_k = k` instead, the
    `k`-th root would tend to `C a e/Delta` and `a` would have to shrink with the order.
15. Lemma 5.2's five moment exponents `(1, 1-2lambda)` and `(2, -2-2lambda, -2lambda)`,
    the change of variable `e_* = c_patch 2^{1/2+lambda}`, and the sign/normalization of
    `d_{U,n}, d_{E,n}` against `B_U alpha_n = -d_{U,n}`, `B_E beta_n = -d_{E,n}`.
16. The equivalence of (5.19) with Theorem 4.6(v)'s `int (H - H_pow) dX = 0` via
    `int r^2 u_{theta,n} dr = q^{3/2-A+lambda_n} int H_n dX`.
17. Lemma 5.4's tail summation `sum_{j>J} 2^{-j} q^{g_j/2 - l'_m} = 2^{-J} q^{g_{J+1}/2 - l'_m}`,
    and that terms with `a_j^{-1} < q` vanish identically so (5.37)'s range restriction is
    never violated.
18. Proposition 5.5 Step 4's geometric split `sum_{n>=J} 2^{-n} q^{nh} <= 2^{1-J} q^{Jh} <=
    2^{1-J} q^{2h(N+1)}` for `J >= 2N+2`.
19. Convergence of the four total moments in Theorem 4.6(v): `M, J` trivially (`U = 0`
    beyond `X_v`), `S` because `E^2 ~ c^2 X^{-1-2h}`, `C_p` likewise, and `I` only after
    subtracting `H_pow` since `H ~ X^{-h}` is not integrable.
20. `Pi(X,eta) = C_p(X,eta) - C_p(inf,eta) = -int_X^inf E^2/(2x) dx`, i.e. (4.25), given
    `C_p(inf,eta) = -Pi_0(eta)`.

---

## 8. Downstream consequences for records already in this Common

- **Supports `ns2026-critique-eq-7-25-missing-weight`.** That critique had to assume,
  without reading Section 4, that (4.27) carries the flat edge weight `zeta` on the
  *derivative* bounds and not only on the lower bound. **(4.27) as printed on p. 33 carries
  `zeta` on both**: `|T_0| >= c zeta` *and* `|d^alpha T_0| <= C_alpha zeta
  delta^{-m_alpha}`. Proposition 5.5's (5.43) has the same shape for the summed stress. So
  the shape the critique argued (7.25) should have is exactly the shape Section 4 uses
  everywhere else - which strengthens the reading that (7.25)'s second bound is a typo
  missing a `zeta`.
- **Partially closes `ns2026-open-question-uniformity-assertions`, external hypothesis (a).**
  See section 3: `kappa` is fixed in profile coordinates before any band index, so
  `kappa > C S_*^{-1/2}` is satisfiable by enlarging `l_0`. Item 7 of that record (uniformity
  in band/label/copy/point) is partially discharged on the profile side by the p. 34 closing
  sentence plus the p. 45 order-of-choices remark.
- **Refutes "the `a_j` are never displayed"** from the Sections 9-10 audit, and
  **confirms "flatness in Lemma 5.4 is conditional"** while showing that the conditionality
  is the ordinary hypothesis/verification split, not a gap.
- **Downgrades Appendix audit flag #1** from high severity: the coverage claim *is* stated as
  a single sentence on p. 40-41 of Section 4, though its support is left to the reader.
- **Confirms Appendix audit flag #2** (derivative-order closure) from the Section 4 side.
- **Resolves Appendix audit flag #8** quantitatively (section 6.6 above).

---

## 9. Artifacts produced

- `commons/ns/artifacts/ns2026_heat_profile_check.py` and
  `ns2026_heat_profile_check_output.txt` - checks `Hcal(0) = 1`, the reduced ODE to `1e-11`
  relative, the full swirl heat equation by centered differences to `5e-5` relative
  (finite-difference truncation dominated), `K_r < 0`, `a = 4/5` on the reference inner
  branch, and the `a <= 2+2h` bound on the heat exterior.
- `commons/ns/artifacts/ns2026_heat_shear_threshold.py` and
  `ns2026_heat_shear_threshold_output.txt` - computes the threshold `Z*` for the asserted
  `-Z Hcal'/Hcal < h/4` of (A.56) and the implied radius bound. (The script's closing prose
  rounds `Z*` to 0.45 and `X` to 4.4; the computed table is the accurate source:
  `Z* = 0.4230, X >= 4.728` at `h = 1e-4`.)

Environment for both: python 3.14.4, numpy 2.5.3, scipy 1.18.1, run as
`./venv/bin/python commons/ns/artifacts/<script>.py` from the scratchpad root.
