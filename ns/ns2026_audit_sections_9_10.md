# Sections 9-10 (pp. 100-126) — "Finite time blowup for Navier-Stokes" (OpenAI, compiled 2026-09-08)

Reviewer scope: pages 100-126 only. Sections 1-8 and the appendices were **not** read;
every claim imported from them is marked `[EXT]` (external dependency, unverified here).
Notation follows the paper; math is in plain text.

---

## 1. Precise statements of the main results

### Section 9 — Residual improvement and the local field

Normalization used throughout (p. 101): physical residual
`R(u,p) = d_t u + (u.grad)u - Laplacian u + grad p`, and `R_* = Q^{2A+1/2} R(u,p)`.
Physical velocity / pressure / residual are recovered from normalized chart coefficients
with factors `Q^{-A}`, `Q^{-2A}`, `Q^{-2A-1/2}`. Throughout Section 9, `kappa_s = 10^{-5}`,
`epsilon = Q^h`, `Q = 2^{-l}`, `q = q(z,t)` (p. 100: the radial integrals hold `(z,t)` fixed,
so **q is a function of z and t only, not of r** — this matters in Section 10).

**Proposition 9.1 (p. 101).** Fix a label, `m` in `Z \ {0}`, `alpha` real. Let `f_m` in `W_alpha`
satisfy the support and smooth-extension hypotheses of Proposition 7.2 `[EXT]`, and assume the
source harmonic `f_m exp(i k m Phi)` extends smoothly by zero outside the fixed local rectangle,
including its pulse endpoints. Let `(t_m, pi_m)` solve `L_m(t_m, pi_m) = -f_m`, where (9.1)

```
L_m(a, pi) = Q^{1+h} N_abs a + K a + eps k^2 m^2 |n_Phi|^2 a + i k m n_Phi . pi.
```

Then the normalized linear residual, after adding the prescribed source `f_m exp(i k m Phi)`,
belongs to `W_{alpha + 1/2 - 3 kappa_s}`, apart from additive pulse-cutoff tails whose Cartesian
space-time derivatives vanish to infinite order as `q -> 0`. The pressure amplitude `pi_m` is in
`W_{alpha + 1/2}`.

The proof expands the linearization as (9.1) plus the remainder (9.2) and reads off the table of
gains (p. 102): slow transport `1 - kappa_s`; phase transport defect `1/2`; base derivatives and
connections `1`; pressure-amplitude gradient `1/2 - kappa_s`; viscous amplitude derivatives
`1 - 2 kappa_s`; mixed derivatives and phase divergence `1/2 - kappa_s`; viscous angular
connection `1/2`. Common lower bound: `alpha + 1/2 - 3 kappa_s`.

**Lemma 9.2 (pp. 102-103).** (i) If `w` is in `W_alpha` and `v` is a mean field with tangential
components in `M_mu` and radial component in `M_{mu+1}`, then the nonzero harmonics of
`(w.grad_*)v + (v.grad_*)w` lie in `W_{alpha + mu - 1/2}`.
(ii) For two curl-generated waves `w, w'` with the **same** label and exponents `alpha, alpha'`,
the nonzero harmonics of `(w.grad_*)w' + (w'.grad_*)w` lie in `W_{alpha + alpha' - kappa_s}`, and
the zero harmonic lies in the mean class `M_{alpha + alpha' - kappa_s}`.
(iii) Distinct labels have zero products on their closed supports (via Lemma 6.1 `[EXT]`).

**Proposition 9.3 (pp. 103-105).** Fix `j` in `N_0`. Let `(u^[j], p^[j])` come from the fixed slow
base and primary waves by a finite sequence of the pulse inverse and curl construction, the signed
amplitude map, and the mean maps of Section 8, with pressure reconstructed by (8.12) after every
update and all products formed from the complete velocity fields. Then there is a decomposition
(9.3) `R(u^[j], p^[j]) = G^[j] + F^[j]` such that, with `G_*^[j] = Q^{2A+1/2} G^[j]` and
`F_*^[j] = Q^{2A+1/2} F^[j]`:

- **(i) Supported sources.** On the coarsest common auxiliary torus, (9.4)
  `G_*^[j] - <G_*^[j]>_theta = sum over gamma in Gamma, m in H_j of f^[j]_{gamma,m} exp(i k_gamma m Phi_gamma)`,
  with `H_j` a subset of `Z \ {0}` of finite cardinality depending on the finite construction but
  **not on the band**; the sum is locally finite;
  `supp(f^[j]_{gamma,m} restricted to E_gamma)` is inside `Omega_gamma`, and
  `f^[j]_{gamma,m}` is in `C^infinity(E_gamma^+; C^3)`, with smooth zero extension across the fixed
  slow and transverse supports and the enlarged local torus rectangle. If `alpha, mu` are lower
  bounds for the wave and mean exponents furnished by Propositions 9.1, 7.2, 6.6, 7.6, Section 8
  and Lemma 9.2 along the sequence, then `f^[j]_{gamma,m}` is in `W_alpha` and
  `<G_*^[j]>_theta` is in `M_mu`.
- **(ii) Flat remainder.** For every fixed `j, m, N`: (9.5) `|F^[j]|_m <= C_{j,m,N} q^N`,
  uniformly over bands, labels and points.
- **(iii) Exact mean equations.** With `(beta^[j], v^[j], gamma^[j], w^[j])` the complete current
  velocity correction as in (8.1) and `E_{B,*} = Q^{2A+1/2} E_B`, (9.6) gives
  `W_ab^[j] = <w_a^[j] w_b^[j]>_theta`, `P^[j] = integral <g_r^[j]>_Y dR`,
  `p_m^[j] = T_0(g_r^[j] - rho P^[j])`, and
  `<G_*^[j] + F_*^[j] - E_{B,*}>_theta = (D_r p_m^[j] - g_r^[j], E_theta^[j], E_z^[j])`.

**Definition 9.4 / Proposition 9.5 (pp. 106-107).** A *state at stage j* has the form (9.7)
`u_* = (b + beta, V + v, G + gamma) + w`, `p_* = p_{B,*} + p_m + p_w`, with
`<w>_theta = <p_w>_theta = 0`, both velocity corrections exactly divergence-free, and (9.8)

```
G_wave^[j] in W_{B_j},   E_theta, E_z in M_{C*_j},   (P, J_theta, J_z) in S_{C*_j},
B_j = 1/2 + sigma_j,     C*_j = 1 + sigma_j,         sigma_j = 1/5 + j/10,
```

together with the *cumulative* bounds (9.9)

```
w in W_{1/2},   w - w_0^tan in W_{0.68},   v, gamma, p_m in M_{0.9},   beta in M_{1.9},
```

and the two exactly preserved moments (9.10)
`integral_0^infinity R^2 <v>_Y dR = 0`, `integral_0^infinity R <gamma>_Y dR = 0`.
Proposition 9.5: the initialization (base plus curls of the primary potentials, then the temporal
mean update (8.20) and the five-equation correction (8.25), with pressure reconstructed by (8.12)
before and after each update) is a state at stage 0 with `B_0 = 0.7`, `C*_0 = 1.2`.

**Proposition 9.6 (pp. 107-111).** Given a state at stage `j`, there is a state at stage `j+1`
with the same base, primary amplitudes, supports and inverse operators, and

```
B_{j+1} = B_j + 1/10,     C*_{j+1} = C*_j + 1/10,
```

obtained by the four *ordered* corrections: (i) solve the inhomogeneous amplitude equation for each
supported nonzero harmonic; (ii) change the wave amplitudes to correct the auxiliary-averaged
tangential residual; (iii) apply the temporal mean inverse to the part with zero auxiliary average;
(iv) apply the five-equation correction to the three current defects. Pressure is reconstructed
after each correction; (9.9) and (9.10) are preserved.

**Lemma 9.7 (pp. 111-112).** There exists `q_big > 0`, independent of the correction stage and of
the derivative order, such that every finite partial sum above (before summation cutoffs) is well
defined on `0 < q < q_big`. For each fixed stage and output amplitude derivative, construction and
estimate require only finitely many input amplitude derivatives. Constants may depend on the stage.

**Lemma 9.8 (pp. 113-114).** With `Z_j = (A_j, B_j, p_j)` and `Delta u_j = curl A_j + B_j` from
(9.15) and the telescoping pressure increment (9.16), for every Cartesian space-time derivative
order `m`: (9.17)

```
|Z_j|_m + |Delta u_j|_m <= C_{j,m} q^{g_j - l_m} (1 + |log q|)^{P_{j,m}},   g_j = (h/10) j -> infinity,
```

with `l_m` independent of `j`; the finite partial sums satisfy
`|(u^[j], p^[j])|_m <= C_{j,m} q^{-K_m} (1 + |log q|)^{P_{j,m}}` with **`K_m` independent of `j`**;
and (9.18)

```
|R(u^[j], p^[j])|_m <= C_{j,m} q^{h sigma_j - K_m} (1 + |log q|)^{P_{j,m}} + E_{j,m},
E_{j,m} <= C_{j,m,N} q^N for every N,
```

again with `K_m` independent of `j`. The per-derivative losses are (9.19): radial physical
derivative `1/2 + h kappa_s`; axial physical derivative `D = 1/2 - h`; time physical derivative
`1 + h`; derivative of the cylindrical frame `1/2`. A sufficient common choice is
`l_m = 2A + (m+1)(1 + 3h/2)`, with `s_x = 1/2 + h/2`, `s_t = 1 + 3h/2`.

**Proposition 9.9 (pp. 114-116).** There are smooth fields `(u_loc, p_loc)` on `Omega_*`, obtained
by summing the finite corrections with shrinking cutoffs, such that `div u_loc = 0` and all
conclusions of Theorem 3.1 hold. Concretely (9.21)

```
A      = A_0 + sum_{j>=1} chi(a_j q) A_j,
B e_theta = B_0 e_theta + sum_{j>=1} chi(a_j q) B_j,
p_loc  = p_0 + sum_{j>=1} chi(a_j q) p_j,      u_loc = curl A + B e_theta,
```

and (9.20) `for all m, N: |R(u_loc, p_loc)|_m = O(q^N)` as `q -> 0`.
Step 3 (p. 115) gives one-sided regularity away from `q = 0`; Step 4 (p. 116) gives the exterior
heat field, with (3.5) `H_ext(s) = 2^A c_infinity H(4s)`, `K(r,tau) = r^{-1-2h} H_ext(tau/r^2)`,
and `sup_{s>=0} |H_ext^{(m)}(s)| <= 2^A c_infinity 4^m (h)_m (1+h)_m < infinity` (from (A.34)
`[EXT]`); Step 5 (p. 116) gives the inner growth
`u_{theta,loc}(sqrt(2 X_in tau), 0, 1 - tau) = tau^{-A} (e_0 + O(tau^{2h}))` with
`e_0 = E_0(X_in, 0) > 0`.

### Section 10 — Compact forcing and whole-space breakdown

Setup (p. 116-117): `u^loc = curl A + B e_theta`, `p^loc`, with `A, B` from (9.21) and `B`
independent of `theta`. Background Stokes streamfunctions (10.1)

```
S_n = q^{1 - A + 2 n h} integral_0^X U_n(X', eta) dX',   S = S_0 + sum_{n>=1} chi(c_n q) S_n,
A_base = (S / r) e_theta,
```

each radial integral vanishing beyond the common support of the axial coefficients (by (4.28) and
(5.10) `[EXT]`, `integral_0^infinity U_n(X, eta) dX = 0` for all `n >= 0`). Near the axis
`S = r^2 a(r^2, z, t)` for a smooth scalar `a`, so
`(S/r) e_theta = a(r^2, z, t) (-x_2, x_1, 0)` — the correct Cartesian smoothness criterion.

**Proposition 10.1 (pp. 117-118).** There are a compact `K` in `R^3`, a smooth vector field `u` and
a smooth scalar `p` on `R^3 x [0,1)`, with `supp u(.,t) union supp p(.,t)` inside `K` for every
`0 <= t < 1`, satisfying (10.2) `div u = 0` and `u(.,t) = p(.,t) = 0` for all sufficiently small
`t >= 0`; they agree with `u^loc, p^loc` on a spatial neighborhood of the origin for all times
sufficiently close to 1.
Key ingredients: (10.3) `q <= C_0 (tau + |z|^{1/D})`, `tau = 1 - t`; a cutoff
`c(x,t) = chi_x(x) chi_t(t)` with `chi_x = 1` for `r <= r_0/2, |z| <= z_0/2`, supported in
`{r < r_0, |z| < z_0}`, and `chi_t = 0` for `1 - t >= tau_0`, `chi_t = 1` for `0 <= 1 - t <= tau_0/2`;
then (10.4) `u = curl(c A) + c B e_theta`, `p = c p^loc`, extended by zero. Product rule:
`u = c u^loc + grad c x A` and `div u = div curl(cA) + r^{-1} d_theta (cB) = 0`.
Force (10.5) `f = R(u,p) = d_t u + (u.grad)u - Laplacian u + grad p` for `0 <= t < 1`, with
`-Laplacian p = sum_{i,j} d_i d_j (u_i u_j) - div f`.

**Lemma 10.2 (pp. 118-120).** For every spatial multi-index `alpha` and integer `j >= 0`,
`d_x^alpha d_t^j f` converges uniformly on `R^3` as `t -> 1`. There are `F_j` in
`C_c^infinity(R^3; R^3)`, all supported in `K`, with (10.6)

```
lim_{t -> 1} d_x^alpha d_t^j f(x,t) = d_x^alpha F_j(x),     d_x^alpha F_j(0) = 0,
```

compatible under spatial and time differentiation, the compatibility being (10.10)
`d_x^alpha d_t^j f(x,t) = d_x^alpha F_j(x) - integral_t^1 d_x^alpha d_t^{j+1} f(x,s) ds`.
Supporting estimates: (10.7) `K(r,tau) = c_infinity s^{-A} H(2 tau / s)`, `s = r^2/2`;
(10.8) `|d_tau^j K(r,tau)| <= C_j r^{-1-2h-2j}` and
`|d_tau^j (K(rho,tau)^2 / rho)| <= C_j rho^{-3-4h-2j}`; and near the origin (10.9)

```
|d_x^alpha d_t^j f(x,t)| <= C_{alpha,j,N} (tau + |z|^{1/D})^N   for every alpha, j, N,
```

on a fixed neighborhood of `(0,1)`.

**Lemma 10.3 (p. 120).** The force (10.5) extends to `f` in
`C_c^infinity(R^3 x (0, infinity); R^3)` with support in `K x [0,2]`. Construction: pick
`chi_0` in `C_c^infinity(R)`, equal to one near zero and zero for arguments at least one; for
`sigma = t - 1 >= 0` set (10.11)

```
f(x, 1 + sigma) = sum_{j>=0} chi_0(b_j sigma) (sigma^j / j!) F_j(x),
```

with (10.12)
`|| d_x^alpha d_sigma^m ( chi_0(b_j sigma) sigma^j F_j / j! ) ||_infinity <= C_{j,m} b_j^{m-j} ||F_j||_{C^{|alpha|}}`,
`b_0 = 1`, and for `j >= 1`

```
b_j = min { b in N : b > b_{j-1},  max over a,m >= 0 with a + m <= floor(j/2) of
            C_{j,m} b^{m-j} ||F_j||_{C^a}  <=  2^{-j} }.
```

Decay: with `K` inside `B(0,R)` and `M_{alpha,m} = || d_x^alpha d_t^m f ||_infinity`,
`|d_x^alpha d_t^m f(x,t)| <= M_{alpha,m} (3+R)^k (1 + |x| + t)^{-k}` for every integer `k >= 0`.

**Lemma 10.4 (p. 121).** Let `F(t) = integral_0^t ||f(s)||_2 ds` for `0 <= t <= 1`. Then
`F(1) < infinity` and (10.13)

```
||u(t)||_2^2 + 2 integral_0^t ||grad u(s)||_2^2 ds <= F(t)^2     (0 <= t < 1).
```

In particular the kinetic energy is uniformly bounded and the total dissipation on `[0,1)` is
finite. Proof route: (10.14) `(1/2) d/dt ||u||_2^2 + ||grad u||_2^2 = <f, u>`; divide by
`(||u||_2^2 + delta^2)^{1/2}`, discard dissipation, Cauchy-Schwarz, integrate, `delta -> 0` to get
`||u(t)||_2 <= F(t)`; then reintegrate (10.14).

**Lemma 10.5 (pp. 121-123).** Fix `T < 1`. If `v, P` is a smooth solution of (1.1) at viscosity one
on `R^3 x [0,T]`, with the force `f` of Lemma 10.3, zero initial velocity, and `v` in
`L^infinity([0,T]; L^2(R^3))`, then `v = u` on that interval, `u` being the localized velocity of
Proposition 10.1. Structure: `w = v - u`, `pi = P - p`,
`g_ij = v_i v_j - u_i u_j = w_i w_j + w_i u_j + u_i w_j`; on `[0,T]`, `||w(t)||_2 <= C_T` and
`sum_{i,j} ||g_ij(t)||_1 <= C_T`; the common force cancels, giving (10.15)

```
d_t w + (v.grad)w + (w.grad)u = Laplacian w - grad pi,     div w = 0.
```

Pressure surrogate (10.16) `pi_* = sum_{i,j} R_i R_j g_ij` (Riesz transforms, multiplier
`-xi_i xi_j / |xi|^2`), shown to lie uniformly in `H^{-s}` for `s > 3/2` and to satisfy
`grad pi = grad pi_*` in the space-time interior. Commutator split (10.18)

```
phi_R^4 pi_* = sum R_i R_j (phi_R^4 g_ij) + sum [phi_R^4, R_i R_j] g_ij,
```

kernel bound `K_R(x-y) = C |x-y|^{-3} min{|x-y|/R, 1}`, `||K_R||_{4/3}^{4/3} <= C R^{-1}`, and with
`E_R = integral chi_R |w|^2`, `A_R = (integral chi_R |grad w|^2)^{1/2}`,
`B_R = ||phi_R^4 w||_6`, `chi_R = phi_R^8`: (10.17) `B_R <= C(A_R + R^{-1} ||w||_2)`; the pressure
flux (10.19)
`|integral pi w . grad chi_R| <= C_T R^{-1} [(B_R + 1) B_R^{1/2} + R^{-3/4} B_R^{3/4}]`;
the difference-energy identity, Young, and
`(1/2) E_R' + (1/2) A_R^2 <= ||grad u||_infinity E_R + C_T / R`; Gronwall
`E_R(t) <= (2 C_T / R) integral_0^t exp(2 integral_s^t ||grad u(a)||_infinity da) ds <= C_T' / R`;
`R -> infinity` gives `w = 0`.

**Theorem 1.1 (proof, p. 124).** The constructed field solves the equation exactly on `[0,1)` by
(10.5); `u` is in `C([0,T]; H^3)` for each `T < 1`; the force is smooth, compactly supported
(Lemma 10.3); the energy is bounded (Lemma 10.4). With the fixed `X_in` of Theorem 3.1, along the
path (10.20) `x_tau = (sqrt(2 X_in tau), 0, 0)`, `t = 1 - tau`, one has `z = 0`, `q = tau`, both
cutoffs equal one for small `tau > 0`, and (10.21)

```
u_theta(x_tau, 1 - tau) = tau^{-A} (e_0 + O(tau^{2h}))  ->  +infinity.
```

`H^3 -> L^infinity` rules out classical `H^3` continuation through `t = 1`; `H^3 -> W^{1,infty}`
plus Gronwall gives uniqueness on shorter intervals, so `[0,1)` is maximal. Any global smooth
solution with the same data and uniformly bounded kinetic energy equals `u` on `[0,1)` by Lemma
10.5, contradicting (10.21). The force is nonzero, since (10.13) would otherwise force `u = 0`.
Viscosity is then arbitrary via the purely spatial rescaling (10.22)

```
u_nu(x,t) = sqrt(nu) u(x / sqrt(nu), t),   p_nu(x,t) = nu p(x / sqrt(nu), t),
f_nu(x,t) = sqrt(nu) f(x / sqrt(nu), t),
```

with `d_x^alpha d_t^m f_nu(x,t) = nu^{(1-|alpha|)/2} (d_y^alpha d_t^m f)(x/sqrt(nu), t)`,
support `K_nu = sqrt(nu) K`, unchanged time, and (10.23)
`||u_nu(t)||_2^2 = nu^{5/2} ||u(t)||_2^2`,
`nu integral_0^1 ||grad u_nu||_2^2 dt = nu^{5/2} integral_0^1 ||grad u||_2^2 dt`. A bounded-energy
competitor for the rescaled data yields the viscosity-one competitor
`v(y,t) = nu^{-1/2} v_nu(sqrt(nu) y, t)`, `P(y,t) = nu^{-1} P_nu(sqrt(nu) y, t)`,
`||v(t)||_2^2 = nu^{-5/2} ||v_nu(t)||_2^2`, already excluded.

**Corollary 10.6 (pp. 125-126).** For every `nu > 0` there is a smooth force on
`T^3 x [0, infinity)`, compactly supported in time, for which the solution `(U, P_per)` of the
periodic Navier-Stokes equations with viscosity `nu` and zero initial velocity is smooth on
`[0,1)` and satisfies `limsup_{t -> 1} ||U(t)||_{L^infinity(T^3)} = infinity`. Velocity and
pressure have support in a fixed compact subset of the interior of
`Q_0 = (-1/2, 1/2)^3` at every time before one; there is no global smooth periodic solution for the
same datum and force. Construction: choose `lambda > 1` with `lambda^{-1} K_nu` compactly inside
`Q_0`, set `t_0 = 1 - lambda^{-2}`, and for `t >= t_0`

```
u~(x,t) = lambda u(lambda x, lambda^2 (t - t_0)),
p~(x,t) = lambda^2 p(lambda x, lambda^2 (t - t_0)),
f~(x,t) = lambda^3 f(lambda x, lambda^2 (t - t_0)),
```

extended by zero for `t < t_0`; then periodize by summing integer translates,
`U = sum_{k in Z^3} u~(x+k, t)`, `P_per = sum_k p~(x+k,t)`, `F_per = sum_k f~(x+k,t)`, the
translates having disjoint supports with a positive gap, so
`(U.grad)U = sum_k (u~(x+k,t).grad) u~(x+k,t)` and

```
d_t U + (U.grad)U - nu Laplacian U + grad P_per = F_per,   div U = 0,   U(.,0) = 0.
```

Growth path `x~_tau = lambda^{-1} sqrt(nu) x_tau`, `t~_tau = 1 - lambda^{-2} tau`, giving
`U_theta(x~_tau, t~_tau) = lambda sqrt(nu) tau^{-A} (e_0 + O(tau^{2h})) -> +infinity`. Periodic
uniqueness on `[0,T]`, `T < 1`, from
`(1/2) d/dt ||w||^2_{L^2(T^3)} + nu ||grad w||^2_{L^2(T^3)} <= ||grad U||_infinity ||w||^2_{L^2(T^3)}`
plus Gronwall, "periodic integration removing both the transport and pressure terms".

---

## 2. Logical dependency chain, and what each stage assumes from Sections 5-8

### 2.1 Skeleton

```
Thm 4.6 / Prop 5.5  [EXT]  fixed slow background (b, V, G), (u_B, p_B), Stokes streamfunctions,
                           expansion (5.1) in powers of q^{2h}, cutoffs chi(c_n q), (5.10), (5.27)
        |
        v
(7.35) primary fields, made divergence-free by Lemma 7.7  [EXT]
        |
        v
Prop 9.5   -----> state at stage 0  (B_0 = 0.7, C*_0 = 1.2)      [uses (8.20), (8.25), (8.12), (8.3)]
        |
        |  induction on j, driven by:
        |     Prop 9.1  (linear residual gain 1/2 - 3 kappa_s)
        |     Lemma 9.2 (nonlinear gains: wave x mean -> +mu-1/2 ; wave x wave -> +alpha'-kappa_s)
        |     Prop 9.3  (the four hypotheses of the next pulse inverse survive a full cycle)
        v
Prop 9.6   -----> state at stage j+1  (B_{j+1} = B_j + 1/10, C*_{j+1} = C*_j + 1/10)
        |
        +-----> Lemma 9.7  (one common domain 0 < q < q_big for ALL stages)
        +-----> Lemma 9.8  (physical Cartesian jet bounds; K_m independent of j)
        |
        v
Lemma 5.4  [EXT, NOT VERIFIED HERE]  summation with shrinking cutoffs chi(a_j q)
        |
        v
Prop 9.9  =  Theorem 3.1 (i)-(iv):  (u_loc, p_loc) on Omega_*, div-free, flat residual (9.20),
                                    exterior heat field (3.5), inner growth (3.6)
        |
        v
Prop 10.1  localization by c = chi_x chi_t  -> compact support, zero datum, (10.2)
        |
        v
(10.5) f := R(u,p)  ---> Lemma 10.2 (uniform derivative limits at t=1)
                    ---> Lemma 10.3 (extension to C_c^infinity, Fefferman decay)
                    ---> Lemma 10.4 (energy bound (10.13))
                    ---> Lemma 10.5 (weak-strong / unconditional-at-infinity comparison)
        |
        v
Theorem 1.1  (viscosity one; then (10.22) for all nu > 0)   and   Corollary 10.6 (periodic)
```

### 2.2 What each stage of the correction cycle assumes

Every stage of Proposition 9.6 consumes, in order:

| Step | Operation | Imported from | What is assumed |
|---|---|---|---|
| 1 | solve the inhomogeneous amplitude equation for each supported nonzero harmonic; increments `Delta w_m = curl_*( i n_Phi x (psi t_m) / (k m |n_Phi|^2) e^{ikm Phi} )`, `Delta p_{w,m} = psi pi_m e^{ikm Phi}` | Prop 7.2 (pulse inverse, `L_m(t_m,pi_m) = -f_m`), Lemma 7.7 (curl correction, `r_m` in `W_{alpha+1/2-kappa_s}`), Prop 9.1, Lemma 9.2 | source in `W_B` with the exact support / smooth-zero-extension hypotheses of Prop 7.2, verified for the *next* cycle by Prop 9.3(i); zero initial data at the entrance to the prescribed path; velocity increment in `W_B`, pressure in `W_{B+1/2}` |
| 2 | change the wave amplitudes to cancel the auxiliary-averaged tangential residual; compactly supported stress `Sigma = Q^{2A}(sigma_2, sigma_1)`; amplitude-correction operator `L_phys^as` | (7.31), (7.33), (7.35), (7.36), Prop 7.5 (primary covariance), Prop 7.6 (stress correction, `W_{alpha-1/2}`), Cor 7.8 (exact covariance expansion), (9.12) | that the total weighted radial moments of `E_theta, E_z` have already been removed (that is (9.12), itself from (8.16) and the preserved constraints (9.10)); a fixed positive lower bound for the primary amplitude `2 sqrt(y_e)`; the moment identity `integral r^e (F_e - b_e M_e) dr = M_e - M_e = 0` |
| 3 | remove the nonconstant auxiliary means by the temporal mean inverse: `Delta v = -c_{i0}^{-1} N_{i0}^{-1} E_theta^o`, `gamma_d = -c_{i0}^{-1} N_{i0}^{-1} E_z^o`, `Psi_* = T_1 gamma_d`, `Delta beta = -eps d_Z Psi_*`, `Delta gamma = (D_r + R^{-1}) Psi_*` | (8.20), Lemma 8.6 (zero-average inverse `N_{i0}^{-1}`, flat axial reconstruction remainder `F_ax`), (8.3), (8.12) | that the part to be inverted has zero auxiliary average (produced by Step 2); `|D^I b| <= C_I S_*^{b_I}` on the shell; (9.9) for the already-accumulated correction |
| 4 | five-equation correction of the three current defects `(P, J_theta, J_z)` | (8.25), (8.26)-(8.27), Lemma 8.7 (fixed linear map), Prop 8.1, Prop 8.3 | that the first two rows preserve (9.10) and the last three remove the linear contributions to the defect changes; fixed support of the test functions; recomputation of pressure from actual fields |

Cross-cutting assumptions used at *every* step: the class calculus of Section 6.4 (`W_alpha`,
`M_alpha`, `S_alpha` at every derivative order); Lemma 6.1 (support separation, hence
distinct-label products vanish); Lemma 6.2 (one common torus for all steps); (6.6) and (6.12)
(chart/graph operators and their derivative costs); (6.23) (logarithmic weight `zeta`, needed for
the boundary integrals on p. 104); (6.28) (enlarged domains `E_gamma^+`); (8.12) (pressure
reconstruction after *every* update — the paper is explicit that all products use the complete
velocity fields); (5.41) (background residual estimates); (5.23) and (5.29) (jet notation and the
tuple `Z_j`).

Closure arithmetic actually displayed (p. 111), with `kappa_s = 10^{-5}` and `B >= 0.7` — I checked
all five numerically:

```
min{1/2 - 3k, 1/2 - k, B - k, 0.4}          = 0.4         >= 0.4          (equality, binding)
min{1/2 - 4k, 0.4 - k, 1/2 - 2k, B - 3k}    = 0.39999     >= 0.4 - k      (equality, binding)
H - B = 1/2 - 2k                            = 0.49998     >  0.1          (margin 0.4)
min{0.17, 1 - 4k}                           = 0.17        >  0.1          (margin 0.07)
0.9 - 4k                                    = 0.89996     >  0.1          (margin 0.8)
```

The *operative* requirement is the `> 0.1` per cycle; the `>= 0.4` lines feed the later steps of the
same cycle. So the induction is **not** tight: the achieved per-cycle gain is `0.4` against a
demand of `0.1`. Consistency checks that pass: `B_j = 1/2 + sigma_j` with `sigma_0 = 1/5` gives
`B_0 = 0.7` and `C*_0 = 1.2` exactly as Prop 9.5 asserts; `B - kappa_s >= 0.69999 > 0.68` matches
`W_{0.68}` in (9.9); "mean-velocity and mean-pressure increments above 0.9" and "radial
mean-velocity increments above 1.9" match `M_{0.9}` and `M_{1.9}`; `beta` in `M_{mu+1}` with
`v, gamma` in `M_mu`, `mu = 0.9`, is exactly the shape Lemma 9.2(i) requires.

---

## 3. Load-bearing estimates, and the steps I judge under-justified

### 3.0 Summary table

| # | Item | Verdict |
|---|---|---|
| A1 | smoothness of `u_loc` for `t < 1` | **Correct, and easier than it looks** (local finiteness) |
| A2 | flat residual (9.20) uniformly | **Not verifiable in Sections 9-10**; rests entirely on Lemma 5.4 `[EXT]` |
| A3 | interplay `C_{j,m}` vs shrinking supports | **Asserted, never displayed**; the quantitative content is off-page |
| A4 | Step 1 exponent-table row label (p. 108) | **Apparent mislabel**; the number 0.4 is right, the label is not |
| A5 | `B_j` notation collision | Readability defect, not an error |
| B | axis and cutoff-transition smoothness for `t < 1` | **Correct**, criterion is the right one |
| C | Lemma 10.2 near the origin given (10.9) | **Correct**, but terse; works with `N = 1` |
| D | Lemma 10.5 | **Correct as written**, and it genuinely closes the growth-at-infinity gap. Best-argued part of my range |
| E | (10.13) and Fefferman's condition | **Correct**; matches Fefferman (D) including the force decay |
| F | Corollary 10.6 pressure normalization | **Essentially correct**; one hypothesis should be stated explicitly in the corollary |

### 3.1 (a) Does the summation over infinitely many stages work?

**Smoothness for `t < 1`: yes, and for a reason that is stated only in passing.** The mechanism is
*local finiteness*, not convergence. Two facts combine:

1. p. 114: "`q >= 1 - t` gives the required local lower bound." So on `{t <= 1 - delta}` the
   coordinate `q` is bounded below by `delta`.
2. `chi(a_j q) = 0` once `a_j q > 1`, and `a_j -> infinity`.

Hence on `{t <= 1 - delta}` only the finitely many `j` with `1/a_j >= delta` contribute to (9.21).
The same argument applies to the background sum `S = S_0 + sum_{n>=1} chi(c_n q) S_n` in (10.1),
since `c_n -> infinity`. So each of `A`, `B e_theta`, `p_loc`, `A_base` is a *finite* sum of smooth
fields on every set `{t <= 1 - delta}` intersected with the local domain, and (9.21) is smooth
there with no estimate required at all. The paper does say "locally finite sum" (p. 115, Step 2) but
does not emphasize that this alone settles interior smoothness; a reader who spends effort trying to
prove `C^infinity` convergence for `t < 1` is doing unnecessary work.

**Flat residual: I could not verify this, and it is the largest single gap in my range.** The
architecture is right and, importantly, is *not* "sum the residuals". Two explicit statements make
that clear:

- p. 105: "The realization step will use a comparison with one fixed finite stage; it will not sum
  these flat residuals directly."
- p. 115: "The tail estimate compares the constructed field with one fixed finite state and its
  flat remainder. The flat residuals separated at each finite stage are estimated through this
  comparison and are not summed."

This is the correct strategy, and the reason it can work is the non-obvious clause in Lemma 9.8:
`K_m` is **independent of `j`**. That is exactly what is needed, because the residual of the sum is
not the sum of residuals — the cross term `(u^[j].grad)(u_loc - u^[j])` pairs the tail against the
partial sum, and the partial sum is only bounded by `C_{j,m} q^{-K_m}(1 + |log q|)^{P_{j,m}}`. If
`K_m` grew with `j`, the tail would have to beat a moving target and the scheme would collapse.
Making `K_m` uniform in `j` is a real design decision and the paper flags it twice (in (9.17) and
again in (9.18)). Credit where due.

What is missing is the actual quantitative step. Section 9.5 Step 1 only *verifies hypotheses*:
"Lemma 9.7 supplies one domain for every finite partial sum before summation cutoffs. The
similarity identities give (5.28), and `q >= 1 - t` gives the required local lower bound. Lemma 9.8
supplies (5.30) with `g_j = h j / 10` and exponent reductions independent of `j`. ... Finally,
(9.18) is (5.33) for the fixed differential polynomial `F(u,p) = R(u,p)`, with
`rho_j = h sigma_j -> infinity`." Every conclusion — the choice of `a_j`, local finiteness, the
`O(q^N)` flatness (9.20) — is then delivered by **Lemma 5.4**, which is in Section 5 and outside my
assignment. `[EXT, UNVERIFIED]`

**The specific thing a reviewer of Section 5 must check** is a quantifier order. The cutoff
parameters `a_j` are chosen once, before any `N` is named; the flatness claim (9.20) is "for all
`m`, `N`". Unwinding the comparison for a given `(m, N)`:

```
|R(u_loc)|_m  <=  |R(u^[j])|_m                              (<= C_{j,m} q^{h sigma_j - K_m} + E_{j,m})
              +  linear terms in the tail  (u_loc - u^[j])
              +  quadratic cross terms     (u^[j] . grad)(u_loc - u^[j])   and   tail x tail
```

Choosing `j = j(m,N)` with `h sigma_j - K_m >= N` handles the first line, and `E_{j,m} <=
C_{j,m,N} q^N` handles its flat part (the `j`- and `N`-dependence of that constant is harmless
because `j` is fixed *after* `N`). But the cross terms require the tail
`sum_{k > j} chi(a_k q) Delta u_k` to be `O(q^{N + K_m})`, and the `a_k` were fixed without
knowing `N` or `m`. So Lemma 5.4's selection must be strong enough to make the tail flat to *all*
orders simultaneously — e.g. `a_k` chosen so that `C_{k,m} q^{g_k - l_m} (1+|log q|)^{P_{k,m}} <=
2^{-k} q^k` for all `m <= k` on `supp chi(a_k q)`. Because `g_k -> infinity` and the supports
shrink, this is achievable (it is the standard Borel / Whitney diagonal device), and I have no
concrete reason to doubt it. But **it is not shown in Sections 9-10, and (9.20) cannot be checked
from these pages.**

**The interplay the caller asked about, made explicit.** Two competing quantities:

- `g_j = h j / 10` (gain), growing *slowly* — `h` is a small parameter (`eps = Q^h`, `D = 1/2 - h`,
  and Lemma 10.2's exponents `r^{-1-2h}` all indicate `h` small).
- `l_m = 2 A + (m+1)(1 + 3h/2)` (per-derivative loss), independent of `j` but growing linearly in
  `m`.

So `g_j > l_m` requires

```
j  >  (10 / h) * [ 2 A + (m+1)(1 + 3h/2) ].
```

For fixed `m` this is a finite but `1/h`-large threshold, and it **grows linearly in `m`**. Read
back into (9.21): for each derivative order `m`, roughly the first `10 l_m / h` stages carry a
*negative* net power of `q` and are controlled purely by the cutoff `chi(a_j q)` restricting them
to `q <= 1/a_j`, i.e. by a *finite* choice made against a stage-dependent constant `C_{j,m}`. That
is precisely why `a_j` must be selected after `C_{j,m}` is known for all `m` up to something that
grows with `j`, and it is why the paper's remark "Constants in the bounds may still depend on the
stage, as in Lemma 9.8" (p. 111) is load-bearing rather than decorative. The paper never writes
down this threshold; I recommend computing it as a sanity check (task V4 below).

**Related flags:**

- **A4 (apparent mislabel, p. 108).** Step 1's table of new nonzero-harmonic exponents reads
  `linear error B + 1/2 - 3 kappa_s`, `cross with old exact waves B + 1/2 - kappa_s`,
  `particular-correction self-interaction 2B - kappa_s`, `cross with total wave correction B + 0.4`.
  The last row cannot be a wave-wave cross term: Lemma 9.2(ii) would give `B + 1/2 - kappa_s`
  (= 0.49999), which is the row above it, and a second wave-wave row would be redundant. The number
  `0.4` is exactly `mu - 1/2` with `mu = 0.9`, i.e. Lemma 9.2(**i**)'s wave-times-mean gain applied
  to `v, gamma` in `M_{0.9}` from (9.9). Step 2's corresponding table (p. 109) confirms this by
  labelling its analogue "mean interaction `B + 0.4 - kappa_s`". So the row should read something
  like "cross with the total *velocity* correction (mean part)". This matters because `0.4` is the
  **binding** entry in the first closure inequality on p. 111 — the one that holds with equality.
  The arithmetic is right; the justification a reader would look up (Lemma 9.2(ii)) is the wrong one.
- **A5 (notation collision).** `B_j = 1/2 + sigma_j` is a residual *exponent* in (9.8) and
  throughout Prop 9.6, while `B_j = b_j e_theta` is an *azimuthal vector field* in (9.15) and is
  used in that sense in (9.21). `B` also denotes `1/2 + sigma_j` inside Prop 9.6's proof, `B_R` is
  a weighted `L^6` norm in Lemma 10.5, `B e_theta` is the summed azimuthal field in Section 10, and
  `b_j` is both the azimuthal coefficient in (9.15) and the cutoff integer in Lemma 10.3. Not an
  error, but (9.21) is genuinely hard to parse and a reader can easily mis-attribute an estimate.
- **A6.** Lemma 9.2(iii) ("distinct labels have zero products on their closed supports") is what
  keeps `H_j` in (9.4) finite and the harmonic range bounded; it is proved in one line by appeal to
  Lemma 6.1 `[EXT]`. Likewise Prop 9.3(i)'s band-independence of `H_j`, supported on p. 113 by "a
  quadratic residual can at most double the current harmonic range, and inverse and curl operations
  preserve each integer". Both are plausible and both are load-bearing for the *uniformity over
  bands* in (9.5); neither is checkable here.
- **A7.** Prop 9.1's "additive pulse-cutoff tails whose Cartesian space-time derivatives vanish to
  infinite order as `q -> 0`" are retained in `F` for the whole iteration, with (9.5)'s constant
  `C_{j,m,N}` depending on `j`. The paper is explicit that these can never be summed, and its
  comparison-with-one-stage device is precisely what avoids needing to. Correct design; but note it
  means the flatness of the *final* field is only as good as Lemma 5.4's tail control, never better.

### 3.2 (b) Smoothness at the axis and at the cutoff transitions for `t < 1`

**Verdict: correct, and the criterion used is the right one.** Three separate places:

1. *Base potential at the axis* (p. 115 Step 2, p. 117 (10.1)). The base Stokes streamfunctions
   have the form `r^2` times a smooth function of `(r^2, z, t)`, so the base vector potential is
   smooth in Cartesian coordinates; the base angular field is `r` times such a scalar times
   `e_theta`. Section 10 makes the criterion explicit: `S = r^2 a(r^2, z, t)` implies
   `(S / r) e_theta = a(r^2, z, t) (-x_2, x_1, 0)`. This is exactly the correct characterization of
   a smooth axisymmetric azimuthal field near the axis (smoothness in `r^2` plus the vanishing
   factor). I verified the identity: `e_theta = (-x_2, x_1, 0)/r`, so
   `(S/r) e_theta = (r^2 a / r^2)(-x_2, x_1, 0) = a(-x_2, x_1, 0)`. Correct.
   The infinite sum in (10.1) is again handled by local finiteness (`c_n -> infinity`, `q >= 1-t`).
2. *Annular corrections at the axis*. p. 114: "All these representatives have smooth zero
   extensions and vanish near the axis for each fixed `t < 1`; their vanishing near the axis gives
   a smooth Cartesian extension there." p. 117 repeats this for the wave and mean potentials of
   (9.21). The neighborhood may shrink with `j` and with `t`, but local finiteness means only
   finitely many terms are present on `{t <= 1 - delta}`, so a common neighborhood exists there.
   Also `curl(Psi_j e_theta) = (-d_z Psi_j, 0, (d_r + r^{-1}) Psi_j)` (p. 114) and
   `d_theta u_loc = 0`, so multiplying `B_j` by a function of `(q, t)` preserves zero divergence —
   I checked that `div(g(q,t) B_j e_theta) = r^{-1} d_theta(g B_j) = 0`. Correct.
3. *Cutoff transition regions.* Two kinds.
   - The *summation* cutoffs `chi(a_j q)`, `chi(c_n q)`: p. 115 Step 3, "Only finitely many terms in
     the expansion in powers of `q^{2h}` and correction stages have nonzero cutoff factors on a
     neighborhood of `c <= q <= c'`, because both cutoff sequences tend to infinity." Correct, same
     mechanism.
   - The *localization* cutoff `c = chi_x chi_t` of Prop 10.1: smoothness of `u = curl(cA) + cBe_theta`
     needs `cA` smooth and extendable by zero, which needs `supp c` to sit at positive distance from
     `{q = q_*}` in the `q` coordinate. That is exactly what (10.3) delivers: on `supp c` we have
     `tau <= tau_0` and `|z| <= z_0`, so `q <= C_0(tau_0 + z_0^{1/D}) < q_*/2` by the choice of
     `tau_0, z_0`. I could not check (10.3) itself, since it cites the coordinates (3.2) `[EXT]`,
     but the two cases quoted are internally consistent with `tau = q(1 - eta^2)` and
     `z = q^D eta`, `|eta| <= 1`: if `1 - eta^2 >= 1/2` then `q = tau/(1-eta^2) <= 2 tau`; otherwise
     `|eta| > 2^{-1/2}` and `q^D = |z| / |eta| < sqrt(2) |z|`, i.e. `q < (sqrt(2)|z|)^{1/D}`. That
     reconstruction reproduces both displayed cases exactly. **Marked: reconstructed, not verified
     against (3.2).**

   One point worth making, because it is easy to get wrong on a first reading: since `q = q(z,t)`
   only, the degeneracy set `{q = 0}` is `{tau = 0, z = 0}` — the *whole plane* `z = 0` at `t = 1`,
   not just the origin. So the cutoff-transition region in `r` (`r` between `r_0/2` and `r_0`, `z`
   near 0) does reach `q = 0` as `t -> 1`, and it is not covered by "endpoint bounds away from
   `q = 0`". The paper handles it by the exterior branch (see (c) below); the reader should notice
   that this branch is doing real work rather than mopping up.

### 3.3 (c) Lemma 10.2's uniform convergence near the origin, given (10.9)

**Verdict: correct, but the write-up compresses the one step where a reader could go wrong.**

The paper's sentence is: "For a prescribed tolerance, we first choose a small ball and then a late
enough time so that this bound is below the tolerance there." Since `C_{alpha,j,N}` in (10.9)
depends on `N`, one cannot let `N -> infinity`; the argument must be run at **fixed `N`**, and
`N = 1` suffices. Made explicit: given a tolerance `eta > 0`, choose `delta > 0` with
`2 C_{alpha,j,1} delta^{1/D} < eta/2` and `delta < min(r_0/2, z_0/2)` (so the cutoffs equal one on
`B(0,delta)`); then for `tau <= delta^{1/D}` and `|x| <= delta`,

```
|d_x^alpha d_t^j f(x,t)|  <=  C_{alpha,j,1} (tau + |z|^{1/D})  <=  C_{alpha,j,1} * 2 delta^{1/D}  <  eta/2,
```

so any two such times give `|d^alpha_x d_t^j f(x,t) - d^alpha_x d_t^j f(x,t')| < eta` on the ball,
and the limit is `0` at the origin, matching `d_x^alpha F_j(0) = 0` in (10.6). Outside the ball,
`K \ B(0,delta)` is compact and avoids the origin, and the two away-from-origin arguments give
uniform Cauchy on a finite cover; off `K` everything vanishes. So uniform convergence on all of
`R^3` follows. The logic is sound.

Two dependencies feed (10.9) and I could only partly check them:

- (10.9) is asserted "Near the origin the cutoffs equal one. We use (3.4) on `X <= X_ext`, and the
  exact zero residual for larger `X`." Where `c = 1` we have `f = R(u^loc, p^loc)` exactly, which is
  flat by (3.4)/(9.20). The conversion from `q`-flatness to `(tau + |z|^{1/D})`-flatness uses (10.3)
  in the **correct direction**: `q <= C_0(tau + |z|^{1/D})` gives
  `q^N <= C_0^N (tau + |z|^{1/D})^N`. Checked. But since (9.20) is itself unverified here (see
  §3.1), so is (10.9). `[EXT chain: (9.20) -> Lemma 5.4]`
- The away-from-origin cases. The paper splits: (i) `q` bounded below and bounded away from `q_*` —
  Theorem 3.1(ii) endpoint bounds `[EXT]` plus the fundamental theorem of calculus
  `sup_{K_0} |d^alpha_x d_t^j G(., t') - d^alpha_x d_t^j G(., t)| <= |t' - t| sup |d^alpha_x d_t^{j+1} G|`;
  (ii) `r` bounded below and `q` small, hence `X >= X_ext`, where `A = 0` and the velocity is the
  explicit exterior azimuthal field `K(r,tau) e_theta` of (3.5), with the pressure given by the
  convergent integral `d_tau^j p^loc(r,tau) = - integral_r^infinity d_tau^j (K(rho,tau)^2/rho) d rho`.
  I verified the exponents in (10.8) against `K(r,tau) = r^{-1-2h} H_ext(tau/r^2)`:
  `d_tau^j` produces `r^{-1-2h} r^{-2j} H_ext^{(j)}`, giving `r^{-1-2h-2j}`; and
  `K^2/rho = rho^{-3-4h} H_ext(tau/rho^2)^2`, giving `rho^{-3-4h-2j}`. Both match. The integrability
  claim ("the second majorant is integrable for `rho >= r_-`") is immediate from
  `rho^{-3-4h-2j}`. Dominated convergence for the differentiated integral is legitimate.
  These two cases do cover the terminal slice minus the origin, by the case split quoted on p. 119
  (`z != 0` gives `q >= |z|^{1/D} > 0`; `z = 0` with `x != 0` gives `r > 0`, hence branch (ii)).

**One internal-consistency check I ran, which passes and is worth recording.** Equating the two
formulas for the same object,

```
(10.7):        K(r,tau) = c_infinity s^{-A} H(2 tau / s),  s = r^2/2   ==>  K = c_infinity 2^A r^{-2A} H(4 tau / r^2)
(3.5), p. 116: K(r,tau) = r^{-1-2h} H_ext(tau / r^2),  H_ext(s) = 2^A c_infinity H(4 s)
                                                                 ==>  K = 2^A c_infinity r^{-1-2h} H(4 tau / r^2)
```

forces `2A = 1 + 2h`, i.e. **`A = 1/2 + h`**. That is consistent with (10.8)'s `r^{-1-2h-2j}` and,
independently, with the growth path: along `r = sqrt(2 X_in tau)` one has
`r^{-1-2h} = (2 X_in)^{-A} tau^{-A}`, reproducing (10.21)'s `tau^{-A}` from the exterior profile.
Three formulas from two different sections agree. A reader should confirm `A = 1/2 + h` against
Section 3's actual definition of `A` `[EXT]`; if it differs, (10.7) and (3.5) are inconsistent.

### 3.4 (d) Lemma 10.5: pressure via Riesz transforms, expanding balls, Gronwall

**Verdict: correct, and it is the most carefully argued item in my range.** I checked essentially
every step and found no error. In particular the growth-at-infinity issue the caller asked about is
real, is explicitly identified by the paper ("These hypotheses leave the growth of spatial
derivatives at infinity unrestricted"; "The pressure term requires separate control because no
spatial growth condition has been imposed on `P`"), and is genuinely resolved rather than finessed.
Details I verified:

1. **Only `L^2` hypotheses are used to get `g` into `L^1`.** `w = v - u` with `v` in `L^infty_t L^2_x`
   and `u` smooth compactly supported gives `||w(t)||_2 <= C_T`. Then
   `g_ij = w_i w_j + w_i u_j + u_i w_j` gives `||g(t)||_1 <= ||w||_2^2 + 2||w||_2||u||_2 <= C_T`.
   No decay of `grad v` or of `P` is needed anywhere. Checked.
2. **`pi_*` is well defined despite Riesz transforms being unbounded on `L^1`.** The multiplier
   `-xi_i xi_j / |xi|^2` is bounded by 1 (its value at `xi = 0` is irrelevant, a null set), and
   `|| g-hat ||_infinity <= ||g||_1`, so
   `|| pi_* ||_{H^{-s}}^2 <= C ||g||_1^2 integral (1+|xi|^2)^{-s} d xi < infinity` for `s > 3/2`.
   Checked; the paper's sentence is exactly right.
3. **`grad pi = grad pi_*`.** For `a` in `C_c^infinity(0,T)`, the conservative form of (10.15) gives
   `integral_0^T a grad pi dt = Laplacian integral_0^T a w dt + integral_0^T a' w dt - div integral_0^T a g dt`.
   I re-derived this: `d_t w + div g = Laplacian w - grad pi` with `div w = div u = div v = 0`;
   multiplying by `a(t)` and integrating in time moves `d_t` onto `a`. The three right-hand terms
   are in `H^{-2}`, `L^2`, `H^{-3}` respectively (the last by the same Fourier estimate with
   `s = 2`), and `pi_*` is in `L_t^infty H_x^{-2}`, so `H_a = integral_0^T a (grad pi - grad pi_*) dt`
   is in `H^{-3}`. Taking the divergence of the difference equation gives
   `Laplacian pi = - sum_{i,j} d_i d_j g_ij = Laplacian pi_*`. I checked the multiplier arithmetic:
   `Laplacian` composed with `R_i R_j` has multiplier `(-|xi|^2)(-xi_i xi_j/|xi|^2) = xi_i xi_j`,
   matching `-sum d_i d_j` whose multiplier is `-(i xi_i)(i xi_j) = xi_i xi_j`. Hence
   `Laplacian H_a = 0`, `H_a-hat` is supported at `{0}` and is a weighted `L^2` function, hence
   zero. This is the standard Liouville step and it is correctly executed. Note that no temperedness
   assumption on `P` is imported: the identity *derives* the distributional class of
   `integral a grad pi dt` from the equation.
4. **Commutator kernel.** `[phi_R^4, R_iR_j] g(x) = integral (phi_R^4(x) - phi_R^4(y)) kappa(x-y) g(y) dy`
   with `|kappa| <= C|x-y|^{-3}`; the delta-multiple in the `R_iR_j` kernel cancels because the
   difference vanishes at `x = y` (the paper says exactly this); and
   `|phi_R^4(x) - phi_R^4(y)| <= C min{|x-y|/R, 1}` since `||grad phi_R^4||_infty <= C/R`. So
   `K_R(x-y) = C|x-y|^{-3} min{|x-y|/R, 1}` is right. I recomputed the `L^{4/3}` norm:
   `||K_R||_{4/3}^{4/3} = C R^{-4/3} integral_0^R r^{-8/3} r^2 dr + C integral_R^infty r^{-4} r^2 dr
   = C R^{-4/3} integral_0^R r^{-2/3} dr + C integral_R^infty r^{-2} dr = 3C R^{-1} + C R^{-1}`,
   i.e. `<= C R^{-1}`, matching the displayed line exactly and giving `||K_R||_{4/3} <= C R^{-3/4}`.
   Young's inequality (`1 + 3/4 = 3/4 + 1`) then bounds the commutator sum in `L^{4/3}` by
   `C_T R^{-3/4}`. All checked.
5. **The `L^{3/2}` sum and the interpolations.** `||phi_R^4 w_i w_j||_{3/2} <= ||phi_R^4 w||_6 ||w||_2
   = B_R ||w||_2` (Hölder `1/(3/2) = 1/6 + 1/2`) and
   `||phi_R^4 w_i u_j||_{3/2} <= ||w||_2 ||u||_6`; Riesz transforms are bounded on `L^{3/2}`. Then
   `||phi_R^3 w||_3 <= ||phi_R^2 w||_3 <= B_R^{1/2} ||w||_2^{1/2}` — I verified via
   `integral |phi_R^2 w|^3 = integral |phi_R^4 w|^{3/2} |w|^{3/2} <= B_R^{3/2} ||w||_2^{3/2}` — and
   `||phi_R^3 w||_4 <= B_R^{3/4} ||w||_2^{1/4}` via
   `integral |phi_R^3 w|^4 = integral (phi_R^4|w|)^3 |w| <= B_R^3 ||w||_2`. Both match the paper.
   With `|grad chi_R| <= C R^{-1} phi_R^7` and the pairing `phi_R^4 pi_*` against `phi_R^3 w`
   (exponents `4 + 3 = 7`), (10.19) follows. Checked.
6. **(10.17).** `chi_R = phi_R^8` is chosen precisely so that
   `||phi_R^4 grad w||_2 = A_R` exactly; then Sobolev plus
   `grad(phi_R^4 w) = phi_R^4 grad w + 4 phi_R^3 (grad phi_R) w` gives
   `B_R <= C ||grad(phi_R^4 w)||_2 <= C(A_R + R^{-1}||w||_2)`. Checked; the eighth power is not
   cosmetic.
7. **The difference-energy identity.** I re-derived it term by term from (10.15) paired with
   `chi_R w` and it is correct as displayed, including the sign and the use of `div v = 0` in the
   transport term and `div w = 0` in `- integral chi_R w . grad pi = integral pi w . grad chi_R`.
8. **Absorption.** With `chi_R = 1` on a neighborhood of `supp u`, `v = w` on `supp grad chi_R`, so
   the transport flux is `<= C R^{-1} integral phi_R^6 |w|^3 <= C_T R^{-1} B_R^{3/2}` (via
   `phi_R^6|w|^3 = (phi_R^4|w|)^{3/2}|w|^{3/2}`), the Laplacian term is `<= C_T R^{-2}`, and by
   (10.17) every power of `A_R` appearing is at most `3/2 < 2`. Young (`R^{-1} A_R^{3/2} <= delta
   A_R^2 + C_delta R^{-4}`) then gives
   `(1/2) E_R' + (1/2) A_R^2 <= ||grad u||_infty E_R + C_T/R`. The `C_T/R` is a crude but valid
   majorant for the collection `R^{-2}, R^{-4}, R^{-4/3}, R^{-7}` at `R >= 1`. Checked.
9. **Gronwall and the limit.** `E_R(0) = 0` (both fields have zero datum), `||grad u||_infty`
   bounded on `[0,T]` because `T < 1` and `u` is smooth with fixed compact support, giving
   `E_R(t) <= C_T'/R`; `chi_R = 1` on each fixed ball for large `R`, so `w = 0`. Checked.

**Does it use only hypotheses stated in Theorem 1.1?** Yes, as far as I can tell from the statement
quoted in my brief ("any smooth bounded-energy solution with the same force and zero datum"). The
lemma uses: `v` smooth on `R^3 x [0,T]`; `v` in `L^infty([0,T]; L^2)`; the same force; zero datum.
It does *not* use decay of `v`, `grad v`, or `P` at infinity, and it does not use `grad v` in `L^2`
globally (only `A_R`, which is a compactly supported quantity, finite by smoothness). The one thing
worth noting is that it is stated at viscosity one; the general-`nu` case is reduced to it on p. 124
by rescaling the *competitor*, using `||v(t)||_2^2 = nu^{-5/2}||v_nu(t)||_2^2`, so the bounded-energy
hypothesis transfers correctly. Checked.

**The two places I would call terse rather than wrong:**
- "These identities hold first for compact smooth cutoffs of `g`, then for `g` by convergence in
  `L^1`, the uniform commutator bound, and convergence of its Riesz transforms in `H^{-s}`. In
  particular `pi_*` is locally integrable in space and time." This is a limiting argument compressed
  into one sentence. It is the right argument (the decomposition (10.18) puts `phi_R^4 pi_*` in
  `L^{3/2} + L^{4/3}`, which is what makes the pairing in (10.19) meaningful), but a referee should
  ask for it in full.
- "for almost every `t` in `(0,T)`" appears twice, inherited from testing the space-time identity
  against `a(t) chi_R w`. Harmless since `E_R` is smooth, but the transition from the a.e.
  differential inequality to Gronwall deserves a clause.

### 3.5 (e) (10.13), "kinetic energy uniformly bounded", and Fefferman's condition

**Verdict: (10.13) is correct — I re-derived it in full — and the statement does match Fefferman's
condition, including the force decay.**

Derivation, checked line by line. (10.14) itself needs `u` smooth with fixed compact support on
`[0,T]`, which Prop 10.1 supplies; then `integral (u.grad)u . u = 0` and
`integral grad p . u = -integral p div u = 0`. Dividing by `(||u||_2^2 + delta^2)^{1/2}`:

```
d/dt (||u||_2^2 + delta^2)^{1/2} = (1/2)(d/dt ||u||_2^2)/(||u||_2^2 + delta^2)^{1/2}
                                <= <f,u> / (||u||_2^2 + delta^2)^{1/2}
                                <= ||f||_2 ||u||_2 / (||u||_2^2 + delta^2)^{1/2}  <=  ||f||_2,
```

so integrating from the zero datum and letting `delta -> 0` gives `||u(t)||_2 <= F(t)`. Returning to
(10.14) and integrating,

```
||u(t)||_2^2 + 2 integral_0^t ||grad u||_2^2 ds = 2 integral_0^t <f,u> ds
                                              <= 2 integral_0^t ||f(s)||_2 ||u(s)||_2 ds
                                              <= 2 integral_0^t F'(s) F(s) ds  =  F(t)^2,
```

using `F' = ||f||_2` and `F(0) = 0`. Exactly (10.13). `F(1) < infinity` because `f` is smooth with
compact space-time support, so `sup_{t<1} ||u(t)||_2^2 <= F(1)^2` — the kinetic energy is uniformly
bounded on `[0,1)`, and monotone convergence gives finite total dissipation. All correct. The
`delta`-regularized division is the standard trick for getting the `F(t)` (rather than `F(t)^2 +
...`) shape and it is used correctly.

**Match with Fefferman.** Two conditions matter.

- *Force decay.* Fefferman's forced alternative requires
  `|d_x^alpha d_t^m f(x,t)| <= C_{alpha,m,K}(1 + |x| + t)^{-K}` for all `alpha, m, K`. Lemma 10.3's
  closing display gives exactly this, in the strongest possible way (compact space-time support):
  `|d_x^alpha d_t^m f(x,t)| <= M_{alpha,m}(3+R)^k (1+|x|+t)^{-k}` for every integer `k >= 0`. The
  paper also notes on p. 120 that "Such a smooth compactly supported force satisfies the decay
  conditions in [13], (5); the zero datum satisfies its initial-data conditions." Match. (I cannot
  check the reference numbering or the erratum text — bibliography is outside my pages `[EXT]`.)
- *Energy.* Fefferman's condition is `integral |u(x,t)|^2 dx < C` for all `t >= 0`. The paper's
  `||u(t)||_2^2 <= F(1)^2` is the same statement up to the factor of `1/2` in "kinetic energy", and
  Lemma 10.5's hypothesis on the competitor is `v` in `L^infty([0,T]; L^2(R^3))` — again the same
  condition, not a stronger one. Match.

Two clarifications worth stating, neither a defect:

- The *logical role* of (10.13) in the negative result is small. What Fefferman's forced alternative
  asks for is the nonexistence of a global smooth bounded-energy solution, and that comes from
  Lemma 10.5 plus (10.21). (10.13) is needed for the *positive* half of Theorem 1.1 (the constructed
  `u` has bounded energy) and for the one-line remark that `f` is nonzero — "since (10.13) would
  otherwise give `u = 0`". That remark is correct and cute: `f = 0` forces `F = 0` forces `u = 0`.
- (10.13) is stated only on `[0,1)`, which is all that exists; there is no claim of an energy bound
  at or past `t = 1`, and none is needed.

**One independent consistency check I ran on the energy, which passes.** With `A = 1/2 + h` and
`D = 1/2 - h` (from (10.19)'s axial loss), the terminal profile on `z = 0` is `~ r^{-1-2h}` by
(3.5), which is *not* locally square-integrable in 3D by itself
(`integral_0 r^{-2-4h} r dr = integral_0 r^{-1-4h} dr` diverges). The construction is saved by the
inner cutoff: the exterior form holds only where `X >= X_ext`, i.e. `r^2 >~ q ~ |z|^{1/D}`. Then

```
integral |u|^2 dx  ~  integral_{|z| <= z_0} integral_{r >~ |z|^{1/(2D)}}^{r_0} r^{-1-4h} dr dz
                   ~  C integral_0^{z_0} |z|^{-2h/D} dz,
```

which converges precisely when `2h / D < 1`, i.e. roughly `4h < 1` for small `h`. So the energy is
finite at the terminal time, consistent with (10.13) — but only barely, and the margin is governed
by `h`. Similarly `integral_0^1 ||grad u||_2^2 dt ~ integral q^{-1/2 - 3h} dq < infinity`,
consistent with "the total dissipation on `[0,1)` is finite". Both checks pass. This is a good sign
and also identifies `h` small as doing structural work well beyond the exponent bookkeeping.

### 3.6 (f) Corollary 10.6 and the pressure normalization

**Verdict: essentially correct, with one hypothesis that should be stated explicitly in the
corollary and one harmless omission.**

What is correct, and checked:

- **The scaling preserves `nu`.** For `u~(x,t) = lambda u(lambda x, lambda^2(t-t_0))`,
  `p~ = lambda^2 p(...)`: `d_t u~ = lambda^3 (d_t u)`, `(u~.grad)u~ = lambda^3 (u.grad)u`,
  `Laplacian u~ = lambda^3 Laplacian u`, `grad p~ = lambda^3 grad p`. Every term scales by
  `lambda^3`, so with `f~ = lambda^3 f(...)` the viscosity is unchanged — exactly as the paper
  states. Checked.
- **Times line up.** `lambda^2(t - t_0) = 1` at `t = t_0 + lambda^{-2} = 1`. Force time support:
  `f` supported in `K x [0,2]` gives `t <= t_0 + 2 lambda^{-2} = 1 + lambda^{-2}`, matching the
  claimed `[t_0, 1 + lambda^{-2}]`. Checked.
- **Smooth gluing at `t_0`.** Extension by zero for `t < t_0` is smooth because the original fields
  and force vanish on an interval after the initial time — which is (10.2) plus `f = R(u,p) = 0`
  where `u = p = 0`. Checked.
- **Periodization is exact.** Distinct translates have disjoint supports with a positive gap, so at
  every point at most one translate is nonzero; derivatives commute with the locally finite sums,
  and the quadratic term periodizes term by term, `(U.grad)U = sum_k (u~(x+k).grad)u~(x+k)`. Checked.
- **Growth path.** `x~_tau = lambda^{-1} sqrt(nu) x_tau` gives `lambda x~_tau = sqrt(nu) x_tau`, and
  `lambda^2(t~_tau - t_0) = lambda^2(1 - lambda^{-2} tau - 1 + lambda^{-2}) = 1 - tau`. Combined
  with `u_nu = sqrt(nu) u(x/sqrt(nu), t)` from (10.22) this gives
  `U_theta = lambda sqrt(nu) u_theta(x_tau, 1-tau) = lambda sqrt(nu) tau^{-A}(e_0 + O(tau^{2h}))`.
  Matches the displayed formula. Checked.
- **Torus decay.** Compact support in time trivially gives
  `||d_x^alpha d_t^m F_per(t)||_infty <= C_{alpha,m,N}(1+t)^{-N}`. Checked.
- **The mean-velocity constraint is automatically consistent.** Integrating the momentum equation
  over `T^3` gives `d/dt <U> = <F_per>`; and indeed
  `integral f dx = integral (d_t u + div(u tensor u) - Laplacian u + grad p) dx = d/dt integral u dx`,
  since every divergence/gradient/Laplacian of a compactly supported field integrates to zero. So no
  hidden obstruction. Checked.

**The pressure point, which is where the corollary is doing something nontrivial.** The paper flags
it: "The pressure is periodic as well, as required in the erratum to the problem statement [13]."
This is the right thing to flag, and it is load-bearing in exactly one place — the final uniqueness
paragraph on p. 126, where "Periodic integration removes both the transport and pressure terms."
The pressure term is removed via `integral_{T^3} grad pi . w = - integral_{T^3} pi div w = 0`, which
requires `pi` to be a genuine periodic function on the torus (no boundary terms). Constructed:
`P_per = sum_k p~(x+k,t)` is periodic by construction. Correct.

Two remarks:

- **The corollary's last sentence should say "with periodic pressure".** As written, "There is no
  global smooth periodic solution for the same datum and force" leaves the competitor's pressure
  class implicit. It matters: if a competitor is allowed `P = P_per + c(t).x` (periodic *gradient*
  but non-periodic pressure — a spatially uniform pressure gradient), then
  `integral grad pi . w` need not vanish and uniqueness genuinely fails, since one can absorb a
  uniform force into a spatially uniform time-dependent shift of `U`. This is precisely why
  Fefferman's erratum imposes periodic `p`, and the paper cites it — so the intent is clear and the
  mathematics is right. But the corollary statement, read on its own, is under-specified, and a
  careless reader could take it as a stronger claim than is proved.
- **Normalization of `P_per` is unfixed and unmentioned.** `integral_{T^3} P_per dx =
  integral_{Q_0} p~ dx` is generally nonzero, whereas the usual periodic convention is mean-zero
  pressure. This is harmless (subtract the mean, which changes nothing in the equation), but it is
  the one thing the phrase "pressure normalization" would most naturally refer to, and the paper
  does not address it.

Also worth noting: the corollary claims the solution's support stays in a fixed compact subset of
the *interior* of `Q_0`, which is what makes the disjointness of translates work and is arranged by
choosing `lambda` with `lambda^{-1} K_nu` compactly inside `Q_0`. Consistent.

---

## 4. How the local construction becomes a genuine whole-space solution with a force

Plain language, in the order the machinery actually runs.

**Start with a shape, not a solution.** Sections 4-5 produce a fixed *background* — an explicitly
designed, self-similar-ish swirling profile living in a curved coordinate system whose scale is a
parameter `q` that shrinks to zero as `t -> 1` at the origin. This background is chosen to have the
growth one wants: along the path `x_tau = (sqrt(2 X_in tau), 0, 0)` at time `1 - tau`, the azimuthal
velocity is about `tau^{-A} e_0` with `e_0 > 0`, so it blows up. Crucially, the growth is a property
of the background *alone*. Everything in Sections 6-9 exists to fix up the fact that the background
does not solve Navier-Stokes; none of it is needed to produce the blowup. In fact the corrections are
supported in an annulus and *vanish identically* in the inner region where the growth is read off
(p. 116, Step 5). That separation is the strategic core of the paper.

**Kill the error, one power at a time.** Plug the background into Navier-Stokes and you get a
leftover — the *residual*. Section 9's job is to make that residual not merely small but *flat*:
vanishing faster than every power of `q` as you approach the singular point. It does this with a
repeating four-move cycle (Prop 9.6): cancel the oscillating (nonzero-harmonic) part of the error by
solving a linear amplitude equation and taking a curl (so the fix is automatically
divergence-free); then adjust the wave amplitudes to cancel the averaged tangential error via a
compactly supported stress; then invert a time-averaging operator to kill the remaining
non-constant averages; then run a five-equation finite-dimensional correction to restore three
integral compatibility conditions that compact support demands. Each pass buys a fixed `1/10` of a
power (`sigma_j = 1/5 + j/10`) while re-establishing every hypothesis the next pass needs — that
last part is Prop 9.3, and it is why the induction is stated as "a state at stage `j`" rather than
as a bare estimate. Each pass also creates new nonlinear errors, which is why the gain has to be
*fixed and positive* rather than merely nonzero: Prop 9.1 and Lemma 9.2 are the ledger that tracks
what every new interaction costs.

**Add up infinitely many fixes without breaking anything.** After `j` passes the error is
`O(q^{h sigma_j})` but the fix itself has grown a stage-dependent constant. So the fixes are summed
with cutoffs `chi(a_j q)` that switch stage `j` on only very close to the singularity — close enough
that its (large) constant is beaten by its (small) power of `q`. Two things make this work. First,
for any time before `t = 1` the coordinate `q` is bounded below, so only finitely many cutoffs are
switched on: the infinite sum is a *finite* sum on every slab `t <= 1 - delta`, and smoothness there
is free. Second, flatness at the terminal time is not obtained by summing the errors — that would
fail, because there are infinitely many of them and each carries its own constant. Instead the
finished field is compared against *one* finite stage, chosen as late as the desired order of
flatness requires, and the difference is controlled because the partial sums lose only a *fixed*
number of powers of `q` per derivative (the `K_m` independent of `j` in Lemma 9.8). The result is
Prop 9.9 / Theorem 3.1: a genuinely smooth, exactly divergence-free field on a shrinking
neighborhood of the singular point, whose Navier-Stokes residual is flatter than any power.

**Cut it out and call the leftover a force.** The field so far lives on a curved local domain, not
on `R^3`. Section 10 multiplies the *vector potential*, the azimuthal coefficient and the pressure —
not the velocity — by a cutoff `c(x,t) = chi_x(x) chi_t(t)`, and only then takes the curl (10.4).
Cutting the potential rather than the velocity is the trick that preserves incompressibility
exactly: `u = curl(cA) + c B e_theta` is divergence-free by construction, for any cutoff. The time
cutoff kills the field at early times, so the initial data is exactly zero. The price is that `u` is
no longer a solution near the edge of the cutoff, where `u = c u^loc + grad c x A` picks up an
error. The paper simply *defines* the force to be whatever is left over: `f := R(u,p)` (10.5). This
is the standard "solve backwards" move — with a force at your disposal, any smooth divergence-free
field is a solution of something — and it is legitimate here only because the resulting `f` turns
out to be genuinely nice. That is the content of the next three lemmas.

**Make the force legal.** A force must be smooth (including at and after the blowup time) and decay.
Two obstacles. (i) `f` is defined only for `t < 1`, and `t = 1` is the singular time. Lemma 10.2
shows every space-time derivative of `f` has a uniform limit as `t -> 1`. Near the origin this is
exactly where the flatness bought in Section 9 pays off: the residual is smaller than any power, so
all derivatives go to zero. Away from the origin it is the explicit exterior heat field — an
elementary closed-form solution — that supplies the bounds. (ii) Having the limits, one still needs
an actual smooth continuation past `t = 1`. Lemma 10.3 builds it as a Borel-type series
`sum_j chi_0(b_j sigma) sigma^j F_j(x) / j!`, with each term switched off after a time `1/b_j` that
shrinks fast enough (chosen recursively) that the series and all its derivatives converge. The
result is `C_c^infinity` in space and time, hence trivially satisfies every decay condition the
official problem statement asks of a force.

**Check the energy, then close the trap.** Lemma 10.4 pairs the equation with `u` and gets
`||u(t)||_2^2 + 2 integral ||grad u||_2^2 <= (integral_0^t ||f||_2)^2`, so the kinetic energy stays
bounded right up to the singular time — the constructed solution is legal, not a large-energy
cheat. The final step is a uniqueness argument, and it has to be unusually careful. The whole point
is to say "no *other* global smooth solution with this force and this datum can exist", so the
competitor must be granted as little as possible: smooth, bounded `L^2` norm, and nothing else — in
particular no control on how fast its derivatives or its *pressure* grow at infinity. Lemma 10.5
does this by (i) reconstructing the pressure difference from the equation via Riesz transforms,
proving it agrees with the unknown pressure difference by a Liouville argument (a harmonic
distribution that is also a weighted `L^2` function must vanish), and (ii) running the energy
estimate on balls of radius `R` and showing every boundary flux — transport and pressure alike — is
`O(1/R)`, with the pressure flux tamed by a commutator estimate. Gronwall then gives `w = 0` on each
ball, and `R -> infinity` gives `v = u`. So a hypothetical global smooth bounded-energy solution
*is* the constructed field, and therefore blows up at `t = 1`. Contradiction. Rescaling `x` by
`sqrt(nu)` (which changes viscosity without touching time) gives every `nu > 0`; summing integer
translates of a shrunken copy gives the periodic version.

---

## 5. Five concrete verification tasks for a newcomer

**V1 (10-30 minutes) — Re-derive (10.13).** From (10.14), reproduce the `delta`-regularized
division step: divide by `(||u||_2^2 + delta^2)^{1/2}`, discard the dissipation, apply
Cauchy-Schwarz, integrate from the zero datum, let `delta -> 0` to get `||u(t)||_2 <= F(t)`; then
return to (10.14) and integrate to get `||u(t)||_2^2 + 2 integral_0^t ||grad u||_2^2 <= F(t)^2`.
Then answer three questions: (i) where exactly is `u(0) = 0` used, and what does the bound become
for nonzero datum? (ii) why is the `delta` needed at all — what goes wrong if you divide by
`||u||_2` directly? (iii) confirm that `F(1) < infinity` uses only compact space-time support of
`f`, so the bound is uniform on `[0,1)`. *Expected outcome:* full agreement; this derivation is
correct as printed.

**V2 (10-30 minutes) — Pin down `A` and cross-check the growth rate three ways.** Equate (10.7)
(`K = c_infinity s^{-A} H(2 tau / s)`, `s = r^2/2`) with the p. 116 form
(`K = r^{-1-2h} H_ext(tau/r^2)`, `H_ext(s) = 2^A c_infinity H(4s)`) and derive `2A = 1 + 2h`, i.e.
`A = 1/2 + h`. Then (i) verify (10.8)'s exponents `r^{-1-2h-2j}` and `rho^{-3-4h-2j}` by
differentiating; (ii) substitute `r = sqrt(2 X_in tau)` into `r^{-1-2h}` and check you recover
`tau^{-A}` as in (10.21); (iii) look up Section 3's definition of `A` and confirm it equals
`1/2 + h` — if it does not, (10.7) and (3.5) are mutually inconsistent and something is wrong.
Finally note that `A = 1/2 + h > 1/2` makes the blowup strictly *faster* than the self-similar rate
`(1-t)^{-1/2}`, hence consistent with Leray's lower bound `||u(t)||_infty >= c(T-t)^{-1/2}` and with
the exclusion of exactly-self-similar blowup.

**V3 (1-3 hours) — Audit the commutator and flux estimates in Lemma 10.5 (pp. 122-123).** Verify, in
order: (i) `|phi_R^4(x) - phi_R^4(y)| <= C min{|x-y|/R, 1}` and that the delta-multiple in the
`R_i R_j` kernel cancels in the commutator, giving
`K_R(x-y) = C|x-y|^{-3} min{|x-y|/R, 1}`; (ii) `||K_R||_{4/3}^{4/3} <= C[R^{-4/3} integral_0^R
r^{-2/3} dr + integral_R^infty r^{-2} dr] <= C R^{-1}`, hence `||K_R||_{4/3} <= C R^{-3/4}`, and
that Young's inequality applies with `g` only in `L^1`; (iii) the two interpolations
`||phi_R^3 w||_3 <= B_R^{1/2}||w||_2^{1/2}` and `||phi_R^3 w||_4 <= B_R^{3/4}||w||_2^{1/4}`;
(iv) (10.17), noting where the eighth power in `chi_R = phi_R^8` is used; (v) assemble (10.19) and
confirm that after substituting (10.17) no power of `A_R` exceeds `3/2`, so Young's inequality
closes against `A_R^2`. Then re-derive the difference-energy identity on p. 123 from (10.15) paired
with `chi_R w` and check every sign. *Expected outcome:* full agreement; I found no error here.
Stretch goal: write out in full the limiting argument compressed into "These identities hold first
for compact smooth cutoffs of `g`, then for `g` by convergence in `L^1`..." — this is the one step in
the lemma that is asserted rather than shown.

**V4 (1-3 hours) — Audit Proposition 9.6's exponent bookkeeping and the summation thresholds.**
(i) For each of the four steps (pp. 108-111), match every table row against the lemma it cites, and
resolve the row "cross with total wave correction `B + 0.4`" on p. 108: show that Lemma 9.2(ii)
would give `B + 1/2 - kappa_s` (already the row above), whereas Lemma 9.2(i) with `mu = 0.9` from
(9.9) gives exactly `B + 0.9 - 1/2 = B + 0.4`, and compare with Step 2's explicitly labelled "mean
interaction". (ii) Evaluate all five closure inequalities on p. 111 at `kappa_s = 10^{-5}`,
`B >= 0.7` and identify which are binding (answer: the first two, with equality) and which have
large margin; confirm the operative per-cycle demand is `0.1`, not `0.4`, so the induction is not
tight. (iii) Confirm the arithmetic of the schedule: `sigma_j = 1/5 + j/10` gives
`B_0 = 0.7, C*_0 = 1.2` as in Prop 9.5, `B_j - kappa_s >= 0.69999 > 0.68` matching `W_{0.68}` in
(9.9), and `H - B = 1/2 - 2 kappa_s`. (iv) Then, using `l_m = 2A + (m+1)(1 + 3h/2)` and
`g_j = h j / 10` from Lemma 9.8, compute the threshold `j > (10/h)[2A + (m+1)(1 + 3h/2)]` above which
stage `j` carries a positive net power of `q` at derivative order `m`, and write down explicitly how
many stages must be controlled by the cutoff alone at order `m`. This last number is the
quantitative content of "shrinking supports" and is never displayed in the paper.

**V5 (multi-session) — Numerically test the growth path (10.20)-(10.21) and the necessary
consistency conditions.** With `A = 1/2 + h`, `D = 1/2 - h` and a small `h` (try `h = 0.01, 0.05,
0.1`):
1. *Growth path.* Implement `x_tau = (sqrt(2 X_in tau), 0, 0)`, `t = 1 - tau`, and
   `u_theta = tau^{-A}(e_0 + O(tau^{2h}))`. Confirm divergence, confirm `|x_tau| ~ tau^{1/2}` so
   that `|u| ~ |x|^{-1-2h}` — the same power as the exterior profile `r^{-1-2h}` in (3.5), which is
   a nontrivial internal consistency check — and confirm the rate beats Leray's
   `c(1-t)^{-1/2}`. Note how slowly `tau^{-A}` diverges for small `h`: this is a *barely*
   supercritical construction, and any numerics must reach very small `tau` to see anything.
2. *Energy.* Model the terminal field as `~ r^{-1-2h}` on `{r^2 >~ |z|^{1/D}}` and zero inside, and
   evaluate `integral |u|^2 dx` over `{|z| <= z_0, r <= r_0}`. You should find it reduces to
   `C integral_0^{z_0} |z|^{-2h/D} dz`, finite iff `2h/D < 1`. Check numerically that this is
   consistent with (10.13) and identify how small `h` must be. Then evaluate
   `integral_0^1 ||grad u||_2^2 dt ~ integral q^{-1/2 - 3h} dq` and confirm it converges, matching
   "the total dissipation on `[0,1)` is finite".
3. *Caffarelli-Kohn-Nirenberg.* A genuine singularity must violate `epsilon`-regularity. Using the
   anisotropic singular region (radial extent `~ sqrt(q)`, axial extent `~ q^D`, time extent `~ q`),
   compute the scale-invariant quantities `r^{-2} integral integral_{Q_r} |u|^3` and
   `r^{-1} integral integral_{Q_r} |grad u|^2` at `r = sqrt(q)`. You should get `~ q^{-4h}` and
   `~ q^{-3h}` respectively — divergent as `q -> 0`, hence compatible with a singularity, but
   diverging only at rate `h`. Verify they actually exceed the CKN threshold `epsilon_0` for the
   relevant range of `q`, and note that the singular set is the single point `(0,1)`, compatible
   with CKN's vanishing one-dimensional parabolic Hausdorff measure. This is the strongest
   independent necessary condition a newcomer can test without reading Sections 4-8, and the fact
   that all three checks pass with matching powers of `h` is meaningful evidence that the
   construction's exponents are internally coherent.

---

## 6. Overall assessment, with explicit uncertainty

### What I checked and believe is correct

- **Lemma 10.4 / (10.13).** Re-derived in full. Correct.
- **Lemma 10.5.** Re-derived essentially every step: the `L^1` bound on `g` from `L^2` hypotheses
  alone; `pi_*` in `H^{-s}` for `s > 3/2` via `||g-hat||_infty <= ||g||_1`; the conservative-form
  identity and its `H^{-2}/L^2/H^{-3}` classification; `Laplacian H_a = 0` and the Liouville step
  (including the multiplier arithmetic); the commutator kernel and `||K_R||_{4/3} <= C R^{-3/4}`;
  the `L^{3/2}` Hölder estimates; both interpolations; (10.17) and the role of `phi_R^8`; the
  difference-energy identity term by term; the Young absorption with maximal `A_R` power `3/2`; and
  the Gronwall conclusion. **Correct.** This lemma is careful, correctly identifies that the
  competitor's pressure and derivative growth at infinity are unconstrained, and resolves that
  properly rather than by an implicit decay assumption. It is the strongest-written part of my range.
- **Prop 10.1's algebra.** `div u = 0` for `u = curl(cA) + c B e_theta` with axisymmetric `c` and
  `theta`-independent `B`; `u = c u^loc + grad c x A`. Correct, and cutting the potential rather
  than the velocity is the right design.
- **Lemma 10.3's construction.** The `b_j` recursion is well posed (`m <= floor(j/2) < j` so
  `m - j < 0` and `b^{m-j} -> 0`); (10.12)'s monomial-degree bookkeeping is right; the tail
  converges in every `C^{alpha,m}`; the `sigma = 0` derivative computation gives `F_m`; support and
  decay claims follow. Correct, modulo the Whitney/Borel gluing theorem being left implicit.
- **Corollary 10.6's scalings.** Verified that all terms scale by `lambda^3` so `nu` is unchanged;
  the time alignment `lambda^2(t-t_0) = 1` at `t = 1`; the force time support
  `[t_0, 1 + lambda^{-2}]`; the growth path composition; and that `d/dt <U> = <F_per>` is
  automatically consistent. Correct.
- **Section 9's exponent arithmetic where checkable.** All five closure inequalities on p. 111 at
  `kappa_s = 10^{-5}`; the consistency of `sigma_j = 1/5 + j/10` with `B_0 = 0.7`, `C*_0 = 1.2`,
  `W_{0.68}`, `M_{0.9}`, `M_{1.9}`; the `mu`/`mu+1` shape matching Lemma 9.2(i).
- **Interior smoothness of the summed field.** Correct, by local finiteness rather than by
  convergence — this is worth stating loudly since it removes a whole class of imagined difficulties.
- **Three independent consistency checks I ran that pass:** `A = 1/2 + h` forced by equating (10.7)
  and (3.5), then confirmed by (10.8) and by the growth path; finiteness of the terminal energy under
  the inner cutoff, requiring `2h/D < 1`; finiteness of `integral_0^1 ||grad u||_2^2 dt`. The CKN
  quantities diverge like `q^{-4h}` and `q^{-3h}`, as a singularity requires.

### What I could not verify, in decreasing order of importance

1. **Lemma 5.4, the summation lemma `[EXT]`.** Every conclusion of Prop 9.9 — the choice of the
   cutoff parameters `a_j`, the local finiteness, and above all the flat residual (9.20) — is
   delegated to it. Section 9.5 only verifies its hypotheses (5.28), (5.30), (5.33). **The
   paper's headline local result (Theorem 3.1) is therefore not checkable from pages 100-126.** The
   specific thing to check in Section 5 is the quantifier order described in §3.1: the `a_j` are
   fixed once, before any flatness order `N` is named, so the tail must be flat to *all* orders
   simultaneously, strongly enough to beat the partial sums' fixed loss `q^{-K_m}` in the quadratic
   cross terms. This is achievable in principle (`g_j -> infinity`, shrinking supports, standard
   Borel/Whitney diagonalization) and the architecture — comparison against one fixed finite stage
   rather than summation of residuals, with `K_m` independent of `j` — is exactly right. But it is
   asserted, not shown, in my range.
2. **Everything imported from Sections 3-8.** Prop 7.2 (the pulse inverse), Lemma 7.7 (the curl
   correction), Props 7.5/7.6 and Cor 7.8 (covariance and stress), (8.20)/(8.25) (the mean update
   and five-equation system), Lemma 8.6 (zero-average inverse), Lemma 8.7, (8.12) (pressure
   reconstruction), Lemma 6.1 (support separation, hence the finiteness of `H_j`), Lemma 6.2, the
   class calculus of Section 6.4, Theorem 4.6 and Prop 5.5 (the background and its growth), (3.2)
   (the coordinates behind (10.3)), (3.4)-(3.6), (A.32)/(A.34). I checked internal consistency
   wherever two of these surface in my pages, and found no contradiction; I checked none of them
   directly.
3. **Prop 9.3(ii)'s uniformity over bands and labels**, and the band-independence of `H_j` in
   (9.4). Load-bearing for (9.5); supported by short arguments (p. 104, p. 113) I cannot audit.
4. **Bibliography and the Fefferman erratum text.** The claim that Lemma 10.3's force satisfies
   "the decay conditions in [13], (5)" and that Corollary 10.6 gives "the periodic breakdown
   alternative (D) in [13]" matches my own recollection of the official problem statement (force
   decay `(1+|x|+t)^{-K}` for all `K`; bounded `integral |u|^2 dx`; periodic pressure per the
   erratum), and Lemma 10.3's closing display and Lemma 10.5's hypothesis do match those. But I did
   not read the references.

### Defects I would report to the authors

- **D1 (substantive, presentational).** p. 108, Step 1's table row "cross with total wave correction
  `B + 0.4`" appears mislabelled: the exponent is Lemma 9.2(**i**)'s wave-times-*mean* gain
  (`0.4 = 0.9 - 1/2` with `mu = 0.9` from (9.9)), not a wave-wave cross term (Lemma 9.2(ii) gives
  `B + 1/2 - kappa_s`, which is the preceding row). Step 2's analogous table labels it "mean
  interaction". This row is the **binding** entry in the first closure inequality on p. 111, so a
  reader chasing the wrong lemma will not be able to confirm the tightest step in the induction.
- **D2 (substantive, statement).** Corollary 10.6's final sentence should state that the competitor's
  *pressure* is periodic. The Gronwall step on p. 126 needs it ("periodic integration removes both
  the transport and pressure terms"), and it fails for a competitor with periodic pressure gradient
  but non-periodic pressure (`P = P_per + c(t).x`). The paper cites the erratum, so the intent is
  clear, but the statement as written is stronger than what is proved.
- **D3 (minor).** `B_j` denotes both a residual exponent (`1/2 + sigma_j`, (9.8)) and an azimuthal
  vector field (`b_j e_theta`, (9.15), used in (9.21)); `b_j` denotes both that azimuthal
  coefficient and Lemma 10.3's cutoff integer; `B_R` is a norm in Lemma 10.5. (9.21) is materially
  harder to read because of this.
- **D4 (minor).** The `N = 1` suffices point in Lemma 10.2's near-origin argument (§3.3) should be
  said, since `C_{alpha,j,N}` depends on `N` and a reader may try to send `N -> infinity`.
- **D5 (minor).** Corollary 10.6 does not fix the mean of `P_per`; harmless, but it is the natural
  reading of "pressure normalization" and deserves a clause.
- **D6 (minor).** The Whitney-type gluing at `t = 1` in Lemma 10.3 (two functions smooth on closed
  half-neighborhoods with all derivatives matching on the interface glue to `C^infinity`) is used
  without being named.

### Bottom line

Within pages 100-126 I found **no mathematical error**, one apparent mislabel that affects the
binding step of Section 9's induction (D1), one statement that is under-specified relative to its
proof (D2), and four presentational defects. The Section 10 material is, on my reading, correct and
in places notably careful — Lemma 10.5 in particular anticipates and closes the growth-at-infinity
loophole that a weaker argument would have papered over, and Lemmas 10.2-10.4 are correct as
printed. The energy bound genuinely matches Fefferman's condition, and the force genuinely satisfies
his decay requirement (indeed trivially, being compactly supported in space-time).

**However, my confidence in Sections 9-10 as a whole is limited by one dependency and it is the
decisive one.** Section 10 is entirely conditional on Theorem 3.1, and Theorem 3.1 is Proposition
9.9, and Proposition 9.9's proof consists of verifying the hypotheses of **Lemma 5.4** and invoking
it. The single most important quantitative claim in my range — that the infinite sum with shrinking
cutoffs yields a residual flatter than every power of `q` (9.20) — is therefore *not established on
these pages*. The design that would make it work is visibly correct in outline (compare against one
finite stage; never sum the flat remainders; keep `K_m` independent of `j`), and I found no sign of
trouble, but I want to be explicit: **a reader who verifies only pages 100-126 has verified that
Section 10 correctly converts Theorem 3.1 into Theorem 1.1 and Corollary 10.6, and has verified the
internal arithmetic of the correction cycle — not that the local field of Theorem 3.1 exists.** The
weight of the paper's claim rests on Section 5's summation lemma and on Sections 6-8's estimate
calculus, which I did not read.

Secondary uncertainty: Section 9's estimates are stated in a bespoke class calculus (`W_alpha`,
`M_alpha`, `S_alpha`, `S_*` polynomial factors, band/label/copy/point uniformity) defined in Section
6.4. I took the class arithmetic at face value. Several of the tightest steps — the `0.4` row, the
`1.49` tangential-mean exponents, (9.12)'s "improve the weighted moments by one" — are one-sentence
appeals to that calculus, and a full referee's report would need to expand each. My checks there
were consistency checks against the stated schedule, not verifications from first principles.
