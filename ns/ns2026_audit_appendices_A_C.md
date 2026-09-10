# Appendices A, B, C of "Finite time blowup for Navier-Stokes" (OpenAI, compiled 2026-09-08)

Reader's report on pages 126-165 of `/Users/benjamincox/dsm:haidaa/navier-stokes.pdf`.
All math is written in plain text. Page numbers are the printed page numbers, which
coincide with the PDF page indices here (page 126 = PDF page 126).

For the dependency mapping I also consulted Section 4 (pages 24-39) for the exact
statements of Theorem 4.6, Lemmas 4.3-4.5, 4.7-4.9, 4.11, Proposition 4.10 and the
definitions (4.1)-(4.36). Everything else in this report is from pages 126-165.

---

## 0. Notation recap (from Section 4, needed to read the appendices)

- Similarity variables: `tau = 1 - t`, `A = 1/2 + h`, `D = 1/2 - h`, `0 < h < 1/2`,
  `z = q^D eta`, `tau = q(1 - eta^2)`, `d = 1 - eta^2`, `L = 1 - 2 h eta^2`,
  `s = r^2/2`, `X = s/q` (equation (4.1), p. 24).
- Leading field: `u_theta = q^{-A} E`, `u_z = q^{-A} U`, `r u_r = V_0`,
  `p = q^{-2A} Pi`, with `E = C^{-1} sqrt(2X) phi` (equation (4.3), p. 25). Axis
  regularity is imposed on `F = phi/C`, `U`, `v_0 = V_0/X`, `Pi` (equations (4.4)-(4.5),
  p. 25).
- `D_X = X d/dX` is the logarithmic radial derivative. `A_X(f) = (1/X) int_0^X f dx`
  is the radial average (4.6).
- `H = sqrt(2X) E` (angular momentum profile), `l = D_X log H`,
  `W = 1 - 2 D eta A_X(U) - d d_eta A_X(U)`, `H_c = D eta + d U` (4.8).
- `Q_s`, `N_s` are the smooth-at-axis solutions of
  `D_X Q_s + (1 + l) Q_s = S_q`, `D_X N_s + N_s = S_n` (4.9), with the explicit
  moment formulas (4.16) (Lemma 4.3, p. 28).
- Five cumulative radial integrals (4.15), p. 28:
  `M = int_0^X U dx`, `I = int_0^X H dx`, `J = int_0^X U H dx`,
  `S = int_0^X (U^2 - E^2/2) dx`, `C_p = int_0^X E^2/(2x) dx`, `Pi = Pi(0,eta) + C_p`.
- Leading stress: `T_0 = F (p_s - s)` with
  `p_s = (X Q_s / L, X N_s / (L E))`, `s = (a, -b_s)`,
  `a = 1 - 2 D_X log E = 2 - 2 l`, `b_s = 2 D_X U / E` (4.11), p. 27.
- Cone coordinates (4.20): `t_s = -b_s/a`, `v_s = a(1 + t_s^2)`,
  `P_c = p_{s,1} + t_s p_{s,2}`, `J_c = p_{s,2} - t_s p_{s,1}`.
- Relaxed cone condition (4.21): `P_c > 2` and `v_s < U(P_c, J_c)` where
  `U(P_c,J_c) = P_c + J_c^2/4 - |J_c| sqrt((P_c-2)/2 + J_c^2/16)`.
  **Admissible** stress cone = relaxed cone **plus** `v_s > 2`.
- Lemma 4.5 (p. 31) equivalent test: for `v_s > 2`, admissible is equivalent to
  `P_c > v_s` and `(v_s - 2) J_c^2 < 2 (P_c - v_s)^2`. Sufficient condition used in
  Appendix A: with `w = p_{s,2}/p_{s,1}`, if `a - b_s w > 0` and
  `2 b_s w + b_s^2/a + (a-2) w^2 < 2` on a compact set, then `p_{s,1}` large enough
  gives the relaxed cone.
- Cone-gap vector (4.35): `Psi(a,b,p) = (a, v-2, c-v, 2(c-v)^2 - (v-2) j^2)` with
  `v = a + b^2/a`, `c = p_1 - (b/a) p_2`, `j = p_2 + (b/a) p_1`. Positivity of all
  four components is exactly `a > 0`, `v_s > 2`, and (4.22).

---

## 1. Precise statements of the main results, with page numbers

### Appendix A - Matching radial moments and constructing the heat exterior (pp. 126-144)

**Lemma A.1 (p. 126).** Let `alpha_1, ..., alpha_m` be distinct reals. For each `j` let
`beta_j` be a nonnegative, nonzero smooth function supported in the interior of a
compact interval `I_j` in `(0, infinity)`, with `I_1 < ... < I_m`. Then the matrix
`B_ij = int_0^infinity x^{alpha_i} beta_j(x) dx` is invertible. The same holds for
exponential weights in a logarithmic coordinate. If the data vary smoothly over a
compact parameter set and the separations persist, `B^{-1}` and each fixed parameter
derivative of it are bounded.
*Proof (p. 127):* a nonzero combination of `m` distinct powers has at most `m-1`
positive zeros (Rolle plus induction), so `det[x_j^{alpha_i}]` has constant sign on
`0 < x_1 < ... < x_m`, and `det B = int_{I_1 x ... x I_m} det[x_j^{alpha_i}] prod_j beta_j dx`
is nonzero. The paper explicitly warns: *"This assertion does not give a uniform
inverse when two exponents approach one another."*

**Quantitative degeneration (A.1), p. 127.** For two identical translated bumps
`beta(y - c_1)`, `beta(y - c_2)` with `Delta = c_2 - c_1 > 0` and
`b(s) = int e^{st} beta(t) dt`,
`det B = b(s_1) b(s_2) (e^{s_2 Delta} - e^{s_1 Delta})`, so `||B^{-1}|| = O(lambda^{-1})`
when `s_1 - s_2 = lambda`. The same `O(lambda^{-1})` holds for weights `1, x^{-lambda}`.

**Lemma A.2 (p. 127) - quadratic moment solve.** Let `K = [-1,1]`, and suppose the
correction coefficients `c` in `R^m` change the required moments by
`F_eta(c) = B(eta) c + Q_eta(c,c) = d(eta)` (A.2), with `B`, `Q` smooth, `B` invertible
on `K`, `beta_0 = sup_K ||B^{-1}||`, `kappa_0 = sup_K ||Q||`, `d_0 = ||d||_{C^0(K)}`. If
`8 beta_0^2 kappa_0 d_0 <= 1`, there is a unique solution in `||c||_{C^0} <= 2 beta_0 d_0`,
smooth in `eta`, selected by iteration from zero of
`c -> B(eta)^{-1}(d(eta) - Q_eta(c,c))`. If `Q = 0` there is no smallness requirement.
For a preselected finite `k`, if also `8 beta_k^2 kappa_k ||d||_{C^k} <= 1` then
`||c||_{C^k} <= 2 beta_k ||d||_{C^k}` (A.3).

**Corollary A.3 (p. 128) - the five-moment Jacobian.** For `(M,I,J,S,C_p)` of (4.15)
under `X = rho x`, the normalizing factors are `(rho, rho^{3/2}, rho^{3/2}, rho, 1)` and
`H = rho^{1/2} sqrt(2x) E` (A.4). On a correction interval where `U_0 = u_c(eta)` and
`E_0 = e_* f_*(eta) x^alpha` with `e_* > 0`, `f_* > 0` smooth, two additive `U` bumps
and three additive `E` bumps give an invertible Jacobian of the normalized five-moment
map **provided `alpha` is not in `{-1/2, 1/2, 3/2}`**. The moment changes have the exact
quadratic form (A.2). At `U_0 = 0`, three `E` corrections can prescribe changes in
`(I, S, C_p)` while preserving `M, J`.
*Proof:* row-reduce `J -> J - u_c I`, `S -> S - 2 u_c M`; the resulting matrix splits
into a `U` block with weights `1, x^{alpha + 1/2}` and an `E` block with weights
`x^{1/2}, x^alpha, x^{alpha - 1}`. Lemma A.1 applies since these are distinct.
Specializations: axis correction `alpha = 1/10` gives blocks `(0, 3/5)` and
`(1/2, 1/10, -9/10)`; intermediate power-law interval `alpha = -1/2 - lambda` gives
`(0, -lambda)` and `(1/2, -1/2-lambda, -3/2-lambda)`. The first of the latter blocks
"can introduce an inverse bound of `lambda^{-1}`".

**Section A.2 (pp. 129-130) - the outer radial profile.** Fixed smooth step (A.5):
`sigma(y) = e^{-1/y^2} / (e^{-1/y^2} + e^{-1/(1-y)^2})` on `0 < y < 1`, `0` for
`y <= 0`, `1` for `y >= 1`. Parameter order (A.6):
`M_d`, `T_d = e^{M_d} + 10`, `P_* > e^{T_d}`, `0 < lambda << 1`,
`0 < h << min{lambda, e^{-T_d}}`; `X_R` fixed afterward. Inner reference data (A.7),
with `x = X/X_R` and `y = log x <= 0`:
`U = 4 eta`, `E = P_* f(eta) e^{y/10} = P_* f(eta) x^{1/10}`,
`f = (1+eta^2)^{-1} = e^{-J_0}`, `J_0 = log(1+eta^2)`.
Exterior target power `E_pow = c_infinity X^{-A}`, `c_infinity > 0` independent of `eta`,
and the two moment identities (A.8):
`int_0^infinity (U^2 - E^2/2) dX = 0`, `int_0^infinity (H - H_pow) dX = 0`,
`H_pow = sqrt(2X) E_pow`.
Stages: axial reduction with `l = (3/5)(1 - sigma(y))` then `l = 0`,
`U = k(y) eta` with `k(y) = 4[1 - sigma(log(1+y)/M_d)]`; intermediate power law
`l = -lambda` on logarithmic length `T_w = 60 log(1/lambda)` (A.9), inside which four
disjoint reserved sub-intervals `(T_w - 25, T_w - 20)`, `(T_w - 20, T_w - 15)`,
`(T_w - 14, T_w - 9)`, `(T_w - 8, T_w - 3)` are set aside for the cone-realization
correction, heat compensation, higher-order background moments, and phase-averaged
corrections. Axial pulse (p. 130): `R_b = Amp(eta) R_0(xi_b) + c_1 beta_1(y) + c_2 beta_2(y)`
with `xi_b = lambda y`, `phi_b(xi_b) = int_0^{xi_b} sigma(v/.02) dv`,
`R_0 = phi_b(xi_b)[1 - sigma(xi_b - 10)]`. Profile interpolation (A.10):
`log E = log e_end - (1/2 + lambda) y - vartheta_f(y) J_0 - (1 - vartheta_f(y)) log 2`,
`vartheta_f(y) = 1 - sigma(y/T_f)`, with the two conditions (A.11)
`I = X H/(1 - lambda)` at the end and `int_{two bumps} (E^2 - E_unedited^2) dy = 0`.
Terminal transition (A.12)-(A.13): `psi_o(y) = 1 - sigma((y-1)/2)`,
`f_o(y) = 1 - rho_o psi_o(y)`, `rho_o = c_o h`, `l = -h + f_o'/f_o`, and
`Q_p = int_0^3 exp(int_0^v (1 + l(s)) ds) (f_o'(v)/f_o(v)) dv`.

**Proposition A.4 (pp. 130-131).** There are finite choices, in the order (A.6), for
which the profile of Section A.2 defines smooth `U`, `E` for `X > 0`,
`eta in [-1,1]`, with `E > 0`, satisfying (A.8). All stage boundaries divided by `X_R`
are independent of `X_R` and `eta`. After the compact axial pulse, `U = M = J = 0`;
after the terminal interval `E = E_pow`. There are four disjoint unaltered patches
(those after (A.9)) on which `U = 0` and `E` is a positive multiple of
`f X^{-1/2 - lambda}`. For each fixed `x_- in (0,1]`, a sufficiently large `X_R`
(allowed to depend on `x_-`) makes the reference profile satisfy the **strict relaxed**
cone condition from `x = x_-` through the first half-unit of its terminal tail, and the
**admissible** cone throughout the intermediate power-law interval, axial pulse, profile
interpolation, exterior transition, and as soon as `l < 0` on the transition into the
constant-slope interval.

**Lemma A.5 (p. 133) - the axis pressure datum.** With `E_id,sched` the reference inner
profile (A.7) followed by the outer profile of Section A.2, with the azimuthal
corrections that preserve the total pressure increment omitted,
`Pi_0(eta) = -(1/2) int_{-infinity}^{infinity} E_id,sched(y,eta)^2 dy` (A.21)
is independent of `X_R` and analytic on a complex neighborhood of `[-1,1]`. It is even,
`Pi_0'(eta)` has the sign of `eta` for `eta != 0`, and
`Pi_0 <= -(5/2) P_*^2 f(eta)^2` (A.22). For the complete reference radial profile,
`Pi(X,eta) = -(1/2) int_{log(X/X_R)}^{infinity} E(v,eta)^2 dv` (A.23).

**Section A.5 (pp. 134-137) - uniform cone on the outer interval.** With
`w = N_s/(E Q_s)` the sufficient test of Lemma 4.5 is (A.24):
`a - b_s w > 0`, `2 b_s w + b_s^2/a + (a-2) w^2 < 2`.
Key explicit lower bound on the reference inner interval and axial transition (A.25):
`Q_s = [(3/5)(4L - 1) - h(1 - 8 eta^2) + (D + 4d) eta J_0'] / (8/5) >= c > 0`.
Axial transition source (A.26): `S_q = -h(1 - 2 k eta^2) + (D + dk) eta J_0'`,
`Q_s >= c(eta^2 + e^{-y})` on `0 <= y <= T_d`. Intermediate interval (A.27) and the
pulse analysis (A.28)-(A.30) culminate in
`w = 2 R_b - C_d d_{xi_b} R_b + o(1)`, `C_d = lambda(D + c_eta)(1-lambda) / (beta^2(lambda - h + c_eta)) <= 2 + o(1)`,
`b_s w <= sup_{R >= 0}(-2R^2 + 2.41 R) + o(1) < .74`,
`2 b_s w + b_s^2/a + (a-2) w^2 <= sup_{R >= 0}(-3.5 R^2 + 4.82 R) + o(1) < 1.68` (p. 136).

**Section A.6 (pp. 137-139) - the heat exterior.** With `a_K = 1 + h`, the heat factor
is
```
(A.32)   Hcal(Z) = (1/Gamma(a_K)) int_0^infinity e^{-v} v^{a_K - 1} (1 + Z v)^{-h} dv,   Z >= 0.
```
`Hcal(0) = 1`; `E_heat(X,eta) = E_pow(X) Hcal(2d/X)`.

**Lemma A.6 (p. 138).** `Hcal` is positive and smooth up to `Z = 0` from the right,
with every fixed derivative bounded on bounded nonnegative intervals. The field
```
(A.33)   K(r,t) = c_infinity s^{-A} Hcal(2 tau / s),   s = r^2/2,
```
is independent of `z`, satisfies `d_t K = (d_rr + r^{-1} d_r - r^{-2}) K`, and has
`K_r < 0`. Its profile is `E_pow(X) Hcal(2d/X)`, smooth with all `eta`-derivatives
extending continuously to `eta = +-1` from the interior.
*Proof ingredients:*
```
(A.34)   Hcal^{(m)}(Z) = ((-1)^m (h)_m / Gamma(1+h)) int_0^infinity e^{-v} v^{h+m} (1 + Zv)^{-h-m} dv
(A.35)   Hcal^{(m)}(0) = (-1)^m (h)_m (1+h)_m
(A.36)   Hcal(Z) = 1 - h(1+h) Z + O_h(Z^2),   |Hcal(Z) - 1| + |Z Hcal'(Z)| <= C h Z
(A.37)   Z^2 Hcal'' + (1 + 2 a_K Z) Hcal' + a_K(a_K - 1) Hcal = 0
```
and the two-sided bounds `0 <= -Z Hcal'/Hcal < h`, `r K_r/(2K) = -A - Z Hcal'/Hcal < -1/2`.
Large-`X` expansion:
```
(A.38)   E_heat(X,eta)/E_pow(X) = Hcal(2d/X) = 1 - 2 h (1+h) d / X + O_{h,j,m}(X^{-2}).
```

**Proposition A.7 (p. 139) - heat replacement with moment compensation.** Let
`X_K = e^2 X_tail` and let `chi_K` be a smooth step, zero for `y <= .2` and one for
`y >= .5`. The terminal replacement
```
(A.39)   E_cl  |-->  E_cl [1 + chi_K(y)(Hcal(2d/X) - 1)]
```
can be compensated by three additive `E`-bumps in the second reserved patch `I_2` of the
intermediate power-law interval. The complete edit preserves `M`, `J`, the totals
`S(infinity)`, `C_p(infinity)` of (4.15), and the renormalized angular moment
`int_0^infinity (H - H_pow) dX`, pointwise on `[-1,1]`. For sufficiently large `X_R`
it preserves `E > 0` and the strict admissible stress cone from that patch through
`y = .5`. The original axis pressure datum and the profiles before the compensation
patch are unchanged. Discrepancy bounds:
```
(A.40)   |d_eta^m Delta int_0^infinity E^2/(2X) dX| <= C_m (e_K^2/X_K) int_1^infinity x^{-2A-2} dx
(A.41)   |d_eta^m Delta S(infinity)| <= C_m e_K^2 int_1^infinity x^{-2A-1} dx
(A.42)   |d_eta^m Delta int_0^infinity (H - H_pow) dX| <= C_m e_K X_K^{1/2} int_1^infinity x^{-A-1/2} dx <= C_{m,h} e_K X_K^{1/2}
```
normalized by `(e_*^2, X_* e_*^2, X_*^{3/2} e_*)` (A.43); the three weights are
`f x_*^{-3/2 - lambda}`, `-f x_*^{-1/2 - lambda}`, `sqrt 2 x_*^{1/2}`.

**Lemma A.8 (p. 140) - backward stress representation.** Suppose the leading
axisymmetric field is smooth at the axis, has `U = 0` and
`E = E_pow f_o [1 + chi_K(Hcal(2d/X) - 1)]` for `y >= 0`, and satisfies the exact moment
conditions
```
M(infinity) = J(infinity) = S(infinity) = 0,   int_0^infinity (H - H_pow) dX = 0,
Pi(X) = -(1/2) int_X^infinity E(x)^2 / x dx.
```
Then the leading tangential stresses vanish for `X >= X_b`; on `y >= 1/2`, where
`chi_K = 1`, they are given by the weighted integrals of the leading residual from `r`
to infinity, with positive integration sign, in (A.46):
`T_theta(r) = r^{-2} int_r^infinity r'^2 Rcal_theta^{(0)}(r') dr'`,
`T_z(r) = r^{-1} int_r^infinity r' Rcal_z^{(0)}(r') dr'`.
Tail representations used: `I(X) = X H_pow/(1-h) - int_X^infinity (H - H_pow) dx`,
`S(X) = (1/2) int_X^infinity E^2 dx`,
`Q_s = -1 + ((1-h) I - D eta I_eta)/(X H)`,
`N_s = (4 h eta S - d S_eta)/X + 4 A eta Pi - d Pi_eta`.
Decay lemma inputs (A.44): `K - K_pow = O_h(tau r^{-3-2h})`,
`d_t K = O_h(r^{-3-2h})`, `r^2 K_r - r K = O_h(r^{-2h})`. Key identities:
`int_0^infinity r^2 (u_theta - K_pow) dr = q^{3/2 - A} int_0^infinity (H - H_pow) dX = 0` (A.45),
`int_0^infinity r(u_z^2 + p) dr = q^{1 - 2A} int_0^infinity (U^2 - E^2/2) dX = 0`.

**Lemma A.9 (p. 141) - flat endpoint factorization.** For `c > 0`, fixed integer
`j >= 0`, and `b(u,eta)` smooth on `0 <= u <= delta_0`, `eta` in a compact `K`, there is
a coefficient `B` with
```
(A.47)   int_0^delta e^{-c/u^2} u^{-j} b(u,eta) du = (1/2) e^{-c/delta^2} delta^{3-j} B(delta,eta)
```
for `0 < delta <= delta_0`, `B` smooth on `[0,delta_0] x K`, `B(0,eta) = b(0,eta)/c`.
*Proof:* substitute `u = delta/(1 + delta^2 v)^{1/2}`.

**Proposition A.10 (pp. 142-144) - stress direction at the outer edge.** Under the
hypotheses of Lemma A.8, the profile on the terminal interval satisfies the admissible
stress cone on `.5 <= y < 3`. Its stress direction extends smoothly to the outer
endpoint, and `2 - (a-2)(T_z/T_theta)^2`, interpreted through this extension, has a
uniform positive lower bound. For `delta = 3 - y` sufficiently small,
```
(A.48)   T_{0,theta} = e^{-4/delta^2} delta^{-3} b_theta(delta,eta),   b_theta(0,eta) > 0
(A.49)   T_{0,z}     = e^{-4/delta^2} delta^{3} b_z(delta,eta)
(A.50)   T_{0,z}/T_{0,theta} = delta^6 b_z/b_theta -> 0
(A.51)   |d^I T_0| <= C_I e^{-4/delta^2} delta^{-N_I},   |T_0| >= c e^{-4/delta^2} delta^{-3}
```
Constants may depend on all fixed profile choices and on `I`, but not on the physical
scale `q`. Exact formulas used:
```
(A.52)   Rcal_theta^{(0)} = K f'/(q L) - K(f_rr + r^{-1} f_r) - 2 K_r f_r
(A.53)   p_z(r) = (2 eta/(q^D L)) int_y^3 K(y')^2 f(y') f'(y') dy'
(A.54)   T_theta = K f_r + r^{-2} (1/(q L)) int_r^infinity (r' K - r'^2 K_{r'}) f_{r'} dr' + (1/(q L r^2)) int_r^infinity r'^2 K f' dr'
(A.55)   |T_z/T_theta| <= C q^{1-D} K = C E_pow Hcal(2d/X)
(A.56)   a = 2 + 2h + 2 Z Hcal'/Hcal - 2 f'/f > 2 + h,   Z = 2d/X
```
and `b_theta(0,eta) = 16 rho_o E_pow(X_b) Hcal(2d/X_b) g(0) / sqrt(2 X_b) > 0`.

---

### Appendix B - Analytic profiles near the axis and their continuation (pp. 144-157)

**Section B.1 - choice of axis data (p. 144).** `Pi_0` is taken from Lemma A.5:
real analytic near `[-1,1]`, even, `Pi_0 <= -c P_*^2 f^2`, `eta Pi_0' > 0` for
`eta != 0`. All schedule parameters, including `0 < h <= 10^{-2}`, are already fixed.
Choose `0 < j_0 <= .05` and set (B.1):
```
U_* = 4 eta + j_0,        H_* = D eta + d U_*,
W_* = 1 - d U_{*,eta} - 2 D eta U_*,
Z_* = -A(1 - 2 eta U_*) U_* - H_* U_{*,eta} - d Pi_{0,eta} + 4 A eta Pi_0.
```
`H_*` has exactly one zero `eta_0` in `(-1,0)`; `H_*/d = D eta/d + 4 eta + j_0` is
strictly increasing from `-infinity` to `+infinity` on `(-1,1)`, and `|eta_0| ~ j_0`.
Moreover `-W_* = 3 - 8 h eta^2 + (1 - 2h) j_0 eta > 2.8` and
`Z_*(eta_0) >= c j_0 P_*^2 > 0`. Choose `delta_* > 0` so that `{|Z_*| <= delta_*}`
avoids a neighborhood of `eta_0`, then `sigma_* > 0` so small that
```
(B.2)   chi = H_*^2 / (H_*^2 + sigma_*^2) > .99   where |Z_*| <= delta_*.
```
Azimuthal axis datum (B.3), for `Lambda >= 1`:
```
zeta_* = -L H_* / (H_*^2 + sigma_*^2),   xi_0 = Lambda zeta_*,
phi_* = exp(Lambda int_0^eta zeta_*(w) dw).
```

**Section B.2 - the analytic coefficient space (pp. 144-146).** `Y = Lambda X`; for
`F(Y,eta) = sum_{alpha >= 0} F_alpha(eta) Y^alpha`,
```
(B.4)   a_{alpha beta} = 20^{-alpha} rho^{-beta} beta! binom(alpha+beta, beta) / ((alpha+1)^2 (beta+1)^2),
        ||F||_rho = sup_{alpha,beta} sup_{eta in I} |d_eta^beta F_alpha(eta)| / a_{alpha beta}.
```
`Bcal_rho` = sequences of smooth coefficients with finite norm; it is complete.
`I F = int_0^Y F(Y',eta) dY'`; `A_X(F) = Y^{-1} I F`.

**Lemma B.1 (p. 145).** Multiplication is bounded on `Bcal_rho`. For `nu = 1, 2` let
`J_nu F` be the solution `G` of `Y G_YY + nu G_Y = F` extending smoothly to `Y = 0` with
`G(0) = 0`; its coefficients are
```
(B.5)   (J_nu F)_{alpha+1} = F_alpha / ((alpha+1)(alpha+nu)),   (J_nu F)_0 = 0.
```
`J_nu` is bounded on `Bcal_rho`, as are radial averaging, multiplication by `Y`, and
`I`, `d_eta I`. In addition
```
(B.6)   || J_nu[(d_eta F)(D_X G)] ||_rho <= C_rho ||F||_rho ||G||_rho,
```
also with either derivative omitted or `F` replaced by `A_X(F)`, and extra
undifferentiated factors admissible at the cost of one algebra constant each.
*Proof tools:* `sum_{i=0}^N (N+1)^2 / ((i+1)^2 (N-i+1)^2) <= 8 sum i^{-2} =: C_sq` (B.7),
the binomial inequality (B.8), and the weight ratios (B.9)-(B.10)
`a_{alpha beta}/a_{alpha+1,beta} = 20 ((alpha+2)^2/(alpha+1)^2)(alpha+1)/(alpha+beta+1) <= 80`,
`a_{i,k+1}/a_{i+1,k} = (20/rho)(i+1)((i+2)^2/(i+1)^2)((k+1)^2/(k+2)^2) <= (80/rho)(i+1)`.

**Section B.3 - the nonlinear analytic axis profile (pp. 146-148).** Scalar comparison
function `f_0(z) = sum_{alpha >= 0} (-z/2)^alpha / (alpha! (alpha+1)!)` (a Bessel-type
series). For `0 <= z <= 4.1`, with `t = z/2`,
```
(B.11)   f_0(z) >= 1 - t/2 + t^2/12 - t^3/144 >= 305719/1152000 > .265,   f_0(z) <= 1.
```

**Proposition B.2 (p. 146) - the analytic axis profile.** With the axis data of Section
B.1, there are `Lambda_0` and, for each `Lambda >= Lambda_0`, an amplitude threshold
`C_0(Lambda)` such that every `C >= C_0(Lambda)` admits an analytic axis profile with
**vanishing leading residual stress** on `0 <= Y = Lambda X <= 4.1`. It has the form
```
(B.12)   phi = phi_* Phi,   U = U_* + Lambda^{-1} u,   Pi = Pi_0 + Lambda^{-1} I(g^2 Phi^2),   g = phi_*/C,
```
with `Phi(0,eta) = 1`, `u(0,eta) = 0`. The profile is analytic in `Y` and in `eta` on a
common neighborhood of `[0,4.1] x I`, and for every fixed `r, s >= 0`,
```
(B.13)   sup_{[0,4.1] x I} | d_Y^r d_eta^s ( Phi - f_0(Y chi), u + Y Z_*/(2L) ) | <= C_{r,s} / Lambda.
```
Constants and a smaller common parameter neighborhood are independent of sufficiently
large `Lambda` and of `C >= C_0(Lambda)`. In unscaled radial derivatives the bound is
`C_{r,s} Lambda^{r-1}`.
*Proof structure (pp. 147-148):* put (B.14)
`B = -2 D eta A_X(u) - d d_eta A_X(u)`, `W = W_* + Lambda^{-1} B`,
`H_c = H_* + Lambda^{-1} d u`, `p = I(g^2 Phi^2)`; the stress-free equations (B.15) are
`-2L(X phi_XX + 2 phi_X)/phi = S_q`, `-2L(X U_XX + U_X) = S_n`, `Pi_X = phi^2/C^2`.
Substitution gives
```
2(Y Phi_YY + 2 Phi_Y) = -chi Phi + Lambda^{-1} R_1,
2(Y u_YY + u_Y)       = -Z_*/L + Lambda^{-1} R_2,
```
with explicit `L R_1`, `L R_2`. Unperturbed solution `Phi^0 = (1+T)^{-1} 1 = f_0(Y chi)`,
`u^0 = -Y Z_*/(2L)`, where `T = J_2 chi / 2` and multiplication precedes integration;
```
||J_2 F||_rho <= 80 ||F||_rho / ((b+1)(b+2)),   ||T^k|| <= (40 M_chi)^k / (k!(k+1)!)
```
so `sum_{k >= 0} (-T)^k` converges absolutely in operator norm and inverts `1 + T`.
The solution map is
`(Phi,u) |-> (Phi^0 + (1+T)^{-1} J_2 R_1/(2 Lambda), u^0 + J_1 R_2/(2 Lambda))`,
a strict contraction on a fixed ball for `Lambda` large. Requirement (B.16):
`C >= sup_{eta in Omega} |phi_*(eta)|`, so `|g| <= 1` on `Omega` and Cauchy's
inequality bounds its derivatives uniformly. Finally the coefficient norm gives a
positive radius of convergence via
`sum_{alpha >= 0} binom(alpha+beta,beta)(R/20)^alpha = (1 - R/20)^{-beta-1}` for `R < 20`,
with `4.1 < R < 20`. The scalar bound (B.11) gives `Phi > 0` for large `Lambda`, hence
`phi > 0`, so division by `phi` in (B.15) is valid.

**Proposition B.3 (p. 148) - positive sources and the outer-endpoint alternative.**
Increasing `Lambda` and then `C_0(Lambda)` if necessary, the normalized azimuthal
profile satisfies `Phi >= c_0 > 0` on `[0,4.1] x [-1,1]`, and on that rectangle
```
(B.17)   S_q >= 2.5 + .95 Lambda L chi.
```
With `p_1 = p_{s,1} = X Q_s/L`, `p_2 = p_{s,2} = X N_s/(L E)`, `n_s = N_s/L`, on
`0 < X <= 4.1/Lambda`,
```
(B.18)   p_1/X >= c_1 > 0,   0 < p_1 <= C_1,   n_s = -2 U_X = Z_*/L + O(Lambda^{-1}).
```
At `X_0 = 4/Lambda` there is a fixed `c_ex > 0` with
```
(B.19)   p_1 + p_2^2/p_1 > 2 + c_ex.
```
`Phi`, `log Phi`, `u`, `p_1`, `n_s` have bounded derivatives in every fixed parameter
order, uniformly in `C >= C_0(Lambda)`.
*Proof highlights (pp. 149-150):*
```
(B.20)   |D_X log Phi| + |H_* d_eta log Phi| <= C chi + C/Lambda
(B.21)   -H_c xi_0 = Lambda L chi + O_{sigma_*}(sqrt(chi)),
         -H_c (log phi)_eta >= .95 Lambda L chi - C_{sigma_*}/Lambda
```
`Q_s(X,eta) = int_0^X X' Phi(Lambda X',eta) S_q(X',eta) dX' / (X^2 Phi(Lambda X,eta))`
is a positive integral representation. Endpoint inequality: at `Y = 4`, in the region
`chi > .99` one has `1.98 < t = 2 chi <= 2` and
`f_0(z) + z f_0'(z) <= 1 - t + t^2/4 - t^3/36 + t^4/576 < -.18`, giving
`-2 z f_0'/f_0 = 2 - 2(f_0 + z f_0')/f_0 > 2.36` and hence `p_1 > 2.3` at the endpoint.
In the complementary region `chi <= .99`, (B.2) gives `|Z_*| > delta_*`, so
`|n_s| >= delta_*/2` and `|p_2| = C sqrt(2/Lambda) |n_s|/(phi_* Phi)` becomes large with
`C`, making `p_2^2/p_1 > 2.3`. Both alternatives give (B.19), e.g. with `c_ex = .2`.

**Section B.5 - reference continuation (pp. 150-151).** `X_0 = 4/Lambda`, `X_b = 100`,
`X_i = 110`; fixed `bar t > 0` with `4 e^{2 bar t} < 4.1`; `y = log(X/X_0)`. For
`t_1 < y < 2 t_1` prescribe (B.22)
```
d_y log phi_r = [1 - sigma((y-t_1)/t_1)] D_X log phi_nat(X,eta),
d_y U_r       = [1 - sigma((y-t_1)/t_1)] D_X U_nat(X,eta),
```
then extend both constantly in `log X`.

**Lemma B.4 (p. 150).** After the outer and axis parameters are fixed, one can take
`Lambda`, then `C`, sufficiently large and then `t_1` sufficiently small so that on
`[X_0, X_i] x [-1,1]`
```
(B.23)   S_{q,r} >= .94 L Lambda chi + 2.4,   l_r <= 1,   v_r := p_{1,r} + p_{2,r}^2/p_{1,r} > 2 + c.
```
Here `p_{1,r} > 0`, and `p_{1,r} > 3` for `X_b <= X <= X_i`. For each fixed `k`, the
parameter norms of `C E_r`, its reciprocal, `p_{1,r}`, `n_{s,r}` on `[X_0,X_i]` have
bounds independent of sufficiently large `C` and of `0 < t_1 <= bar t`.
Key relations: (B.24) `Pi_r(X) - Pi_0 = C^{-2} int_0^X phi_r(X',eta)^2 dX'`,
`sup_{X <= X_i} ||Pi_r - Pi_0||_{C_eta^k} <= C_k C^{-2}`;
(B.25) `D_X p_{1,r} = X S_{q,r}/L - l_r p_{1,r}`, `D_X n_{s,r} + n_{s,r} = S_{n,r}/L`;
comparisons `p_{1,r}(y) >= p_{1,r}(0) e^{-y} + 1.88 chi (e^y - e^{-y})` and
`p_{1,r}(X) >= p_{1,r}(X_0) X_0/X + 1.2(X - X_0^2/X)`.

**Proposition B.5 (p. 151) - continuation to nonzero leading stress.** The analytic
axis profile can be continued to `X_i` with positive `E`, unchanged for `X <= X_0`.
Its stress has the factorization `T_0 = e_a B_0`, where the scalar `e_a` vanishes to
infinite order at `X_0` and the smooth vector coefficient `B_0` is nonzero there. The
**admissible** stress cone holds on a first collar, and the **strict relaxed** cone
condition holds on the rest of `(X_0, X_i)`. At `X_i`: `a = .8`, `D_X U = 0`, `p_1 > 2`.
*Construction (B.26), p. 152:* with `0 < kappa_0 < 1/2`,
```
e_a = (1 - kappa_0) sigma(y/t_1),   kappa = 1 - e_a,   a = kappa p_{1,r},
D_X U = -kappa X n_{s,r}/2,        d_y log phi = -kappa p_{1,r}/2,
(B.27)   d_y(log phi - log phi_r) = e_a p_{1,r}/2,   d_y(U - U_r) = e_a X n_{s,r}/2,
(B.28)   p_s - p_{s,r} = O(y e_a),   E_r/E = 1 + O(y e_a),
(B.29)   t_s = (p_{2,r}/p_{1,r})(E_r/E),   v_s = kappa v_r + O(y e_a),
         P_c = v_r + O(y e_a),   J_c = O(y e_a),
         (v_s - 2)_+ J_c^2 < 2(P_c - v_s)^2   for 0 < y <= t_1.
(B.30)   T_0 = e_a B_0(y,eta),   B_0(0,eta) = F(X_0,eta) p_{s,r}(X_0,eta) != 0,   B_0 in C^infinity
(B.31)   |T_0| >= c e^{-t_1^2/y^2},   |d^I T_0| <= C_I e^{-t_1^2/y^2} y^{-N_I}
```
`e_a = e^{-t_1^2/y^2} g_a(y)` near zero with `g_a` smooth and positive; Lemma A.9 with
`c = t_1^2`, `j = 0` supplies `(1/e_a(y)) int_0^y e_a(u) b(u,eta) du = y^3 B(y,eta)`,
`B` smooth, which is what makes division by the flat factor smooth.
Completion of the continuation: keep (B.26) with `kappa = kappa_0` through the reference
cutoff and constant-slope interval; at `X_b` multiply `D_X U` by a smooth factor `beta`
decreasing from one to zero, then interpolate `a` convexly from `kappa_0 p_{1,r}` to `.8`.
On the final constant-slope interval `a = .8`, `l = .6`, and the gradient absorption of
Lemma B.4 gives `S_q > 1` (constant contribution `.6 (-W_*) > 1.68`), so at a
hypothetical downward crossing `p_1 = 2` one would have
`D_X p_1 = X S_q/L - .6 p_1 > X - 1.2 > 0`, a contradiction.

**Corollary B.6 (p. 153) - the reserved analytic collar.** There is `t_c > 0` with
`t_c < t_1` and `4 e^{t_c} < 4.1` such that `E/sqrt(2X)`, `U`, `V_0/X`, `Pi` on
`[0, X_0 e^{t_c}] x [-1,1]` are smooth in `(X,eta)` through `X = 0` and analytic in
`eta` on one complex neighborhood; every fixed radial derivative has the same parameter
neighborhood and finite bounds on a smaller one. On `X_0 < X <= X_0 e^{t_c}` the
admissible cone condition and the factorization (B.30) hold. All subsequent profile
modifications are supported strictly to the right of `X_0 e^{t_c}`.

**Lemma B.7 (p. 154) - transition length before the amplitude.** Fix the outer- and
axis-profile data, `Lambda`, and a finite parameter order `k`. With
`ell_i = log(C E(X_i,eta))`, `G_i = U(X_i,eta)`, and `omega_fin` the sum of the
logarithmic radial lengths of the final axial cutoff and of the interpolation of `a`
to `.8` in Proposition B.5,
```
(B.32)   ||ell_i||_{C_eta^k} <= Bcal_k,
         ||G_i - U_*||_{C_eta^k} <= C_k^nat / Lambda + C_k^join(Lambda)(t_1 + kappa_0 + omega_fin).
```
`C_k^nat` is independent of both `Lambda` and sufficiently large `C`; `Bcal_k` and
`C_k^join(Lambda)` may depend on `Lambda` but not on `C` or on smaller transition widths.
Transition choice (B.33)-(B.34):
```
(B.33)   T_sh >= 20 ||sigma'||_infinity (Bcal_0 + ||log f||_infinity)
(B.34)   log E = -log C + y_i/10 + (1 - sigma(y_i/T_sh)) ell_i + sigma(y_i/T_sh) log f,   U = G_i,
```
with `y_i = log(X/X_i)`. Then `l = .6 + T_sh^{-1} sigma'(y_i/T_sh)(log f - ell_i)` lies
in `[.55, .65]`, so `a` is in `[.7,.9]` and `b_s = 0`.

**Section B.8 - matching the five radial moment functions (pp. 155-157).**
`X_R = X_i (C P_*)^{10}`, `x = X/X_R`, `X_sep = X_i e^{T_sh}`.

**Proposition B.8 (p. 155).** For sufficiently large amplitude `C`, followed by
sufficiently small activation transitions, the profile of Proposition B.5 continued by
(B.34) can be continued to `log x = -5`, preserving the strict relaxed cone condition,
so that its fields, pressure `Pi`, and all five radial moment functions `M, I, J, S, C_p`
agree exactly with those of the reference inner profile (A.7). It can then follow the
reference outer radial profile without changing its prescribed pressure datum or the
subsequent values of `Q_s`, `N_s`.
*Normalization (p. 155):* `Hhat = sqrt(2x) E`,
`M = X_R Mhat`, `I = X_R^{3/2} Ihat`, `J = X_R^{3/2} Jhat`, `S = X_R Shat`,
`C_p = int_0^x E^2/(2 xi) dxi`, `Pi = Pi_0 + C_p`, and (B.35) rewrites `W`, `Q_s`, `N_s`
in terms of the hatted moments with `X_R` eliminated. Tolerance (B.36) on
`R = [e^{-8}, e^{-5}]`:
```
Q_s >= Q_min/2,   E >= e_*/2,   |a - .8| <= .1,   |N_s/(E Q_s)| <= w_*,   |b_s| <= .1/(1+w_*)
```
whence `v_s <= .9 + .01/.7 < 1` and (B.37) `G = Q_s - b_s N_s/(a E) >= (6/7) Q_s >= 3 Q_min/7`,
`P_c = X_R x G/L`, so `P_c > 2` follows by increasing `X_R`. Endpoint coordinate (B.38):
`x_sep = X_sep/X_R = e^{T_sh}/(C P_*)^{10}`. Inner-interval contributions (B.39):
```
||Mhat||_{C_eta^k} + ||Shat||_{C_eta^k} <= C_k x_sep,
||Ihat||_{C_eta^k} + ||Jhat||_{C_eta^k} <= C_k C^{-1} x_sep^{3/2},
||C_p||_{C_eta^k} <= C_k C^{-2}   (x = x_sep).
```
Restoration: on `-8 < log x < -7` restore `G_i` smoothly to `4 eta` keeping `E` ideal;
on `-6 < log x < -5` add two fixed bumps to `U` and three to `E`. Replacing the `J` row
by `J - 4 eta I` and the `S` row by `S - 8 eta M`, the normalized moment map has blocks
```
U weights: 1, f x^{3/5};      E weights: x^{1/2}, f x^{1/10}, f x^{-9/10}
```
which is Corollary A.3 with `alpha = 1/10`.

**Remark B.9 / (B.40), p. 157 - order of choices.**
```
M_d, T_d, P_*, lambda, h  ->  tolerance for matching moments,  j_0  ->  delta_*, sigma_*, Lambda
  ->  (Bcal_k), T_sh  ->  C, X_R  ->  kappa_0, t_1, widths of final transitions.
```

**Corollary B.10 (p. 157).** The analytic axis profile of Proposition B.2 admits a
smooth continuation with `E > 0` for `X > 0`, matching the reference outer radial
profile at `log(X/X_R) = -5` with exact fields, pressure, and radial moment functions.
It has vanishing leading residual stress through `X_0`, followed by an inner collar on
which the profile remains analytic in `eta` and satisfies the admissible stress cone
condition with the flat factor (B.30). It satisfies the strict relaxed cone inequalities
thereafter up to the matching point. Its finite parameters can be chosen in the order
(B.40).

---

### Appendix C - Realizing the admissible stress cone (pp. 157-165)

**Setup (pp. 157-158).** The joined profile of Corollary B.10 + Proposition A.4 +
Proposition A.7 has nonzero leading stress `T_0` exactly on `(X_a, X_b)` with
`eta`-independent endpoints. It satisfies the **strict relaxed** cone throughout, and
the **admissible** cone on an inner collar and from the intermediate power-law interval
onward. `I = [X_-, X_+]` is a compact subinterval of `(X_a, X_b)` containing every point
where the admissible condition may fail, starting inside the first collar and ending in
the intermediate power-law interval before its first reserved patch. `X_an` denotes the
reserved analytic-rectangle endpoint from Corollary B.6, and `X_-` is chosen in the
shorter inner collar.

**Lemma C.1 (pp. 158-160) - a loop with the prescribed mean.** Let `a, b_s, p_s` be
smooth on `I x [-1,1]` with `a > 0`, and set `t_s = -b_s/a`, `v_s = a(1 + t_s^2)`.
Suppose the strict inequalities of the relaxed cone condition (4.21) hold everywhere,
and the admissible stress cone holds in neighborhoods of the two radial boundary
components. There is a smooth loop `(a_L, -b_L)(X,eta,phi)`, of period one in `phi`,
with `a_L > 0`, such that
```
(C.1)   int_0^1 (a_L, -b_L) dphi = (a, -b_s).
```
With `Psi` the vector of cone gaps in (4.35), there are constants `kappa_L, delta_d > 0`
depending on the fixed data `a, b_s, p_s, I` such that
```
(C.2)   Psi_j(a_L(X,eta,phi), b_L(X,eta,phi), p_s(X,eta)) >= kappa_L,
        (X,eta,phi) in I x [-1,1] x (R/Z),   j = 1,2,3,4
(C.3)   (a_L, -b_L)(X,eta,phi) = (a, -b_s)(X,eta)   for dist(X, dI) < delta_d.
```
*Construction:* write the shear as a positive multiple of `(1,t)` and vary the ratio `t`
with mean `t_s`; seek admissible vectors `v(1,t)/(1+t^2)` with a scalar `v > 2`, so the
**variance** of `t` determines `v`. With `<f>_{theta'} = (2 pi)^{-1} int_0^{2 pi} f dtheta'`,
`P_c(t) = p_{s,1} + p_{s,2} t`, `J_c(t) = p_{s,2} - p_{s,1} t`,
`0 < d_0 < (1/2) min(P_c(t_s) - 2)`,
```
(C.4)   M_e(z) = <e^{z sin theta'}>_{theta'}
(C.5)   t(theta'; mu) = t_s + d_0 (e^{mu p_{s,2} sin theta'}/M_e(mu p_{s,2}) - 1)/p_{s,2},   mu >= 0
```
(removable singularity at `p_{s,2} = 0` via `G(p)/p = int_0^1 G'(mu p) du`), giving
`<t>_{theta'} = t_s`, `P_c(t) >= P_c(t_s) - d_0 > 2`, and
`V(mu, p_{s,2}) = <(t - t_s)^2>_{theta'} = (d_0^2/p_{s,2}^2)[M_e(2 mu p_{s,2})/M_e(mu p_{s,2})^2 - 1]`,
equal to `d_0^2 mu^2/2` at `p_{s,2} = 0`. Monotonicity: with `g = log M_e`,
`g'' > 0` (variance of `sin theta'` under a positive tilted density) and
```
(C.6)   d_mu{g(2 mu p) - 2 g(mu p)} = 2 p {g'(2 mu p) - g'(mu p)} > 0
```
so `V` increases strictly in `mu`; `M_e(z) ~ e^{|z|}/sqrt(|z|)` for `|z| >= 1` gives
`V(mu,p) -> infinity`, hence a finite `mu_max` with
```
(C.7)   V(mu_max, p_{s,2}) > 3/min a   throughout I x [-1,1].
```
Upper admissible bound `Ucal(P_c, J_c)` (the smaller root in `v` of
`2(P_c - v)^2 = (v-2) J_c^2`) is continuous and `> 2` on the compact family, so there is
`0 < delta_L < 1` with
```
(C.8)   Ucal(P_c(t), J_c(t)) > 2 + delta_L   for 0 <= mu <= mu_max.
```
Cutoff (C.9): `zeta_L(v_s)` smooth, in `[0,1]`, `= 1` for `v_s <= 2 + delta_L/8`, `= 0`
for `v_s >= 2 + delta_L/4`;
```
(C.9)   v_* = 2 + delta_L/2,   rho = zeta_L(v_s)^2 (v_* - v_s),   v = v_s + rho
(C.10)  V(mu, p_{s,2}) = rho/a,   0 <= mu < mu_max
```
solved uniquely by (C.7) and strict monotonicity; smooth through `rho = 0` because the
signed square root of `V` is smooth and odd near `mu = 0` with derivative `d_0/sqrt 2`.
Then `v = v_*` where `v_s <= 2 + delta_L/8`, `2 < v_s <= v <= v_*` in the transition, and
`v = v_s > 2` where `v_s >= 2 + delta_L/4`; where the cutoff vanishes, `mu = 0` and the
loop equals the given data. Finally reparametrize `theta'` by a lifted angle `phi`:
```
dphi/dtheta' = a(1 + t^2)/(2 pi v),   (a_L, -b_L) = v(1,t)/(1 + t^2)
```
which has period one because `a<1+t^2>_{theta'} = v_s + a V = v`, and
`int_0^1 a_L dphi = a`, `int_0^1 (-b_L) dphi = a<t>_{theta'} = a t_s = -b_s`.

**Proposition C.2 (pp. 161-163) - radial modulation with phase `N log X`.** For the
fixed input profile `E, U, Pi`, a sufficiently large finite integer `N` and a correction
supported in the **first** reserved patch of the intermediate power-law interval give
smooth profiles `Etil, Util, Pitil` with `Etil > 0` for `X > 0`, satisfying the
**admissible** stress cone throughout `(X_a, X_b)`. Both endpoint collars, the rectangle
near the axis, and the two patches reserved for higher-order and mean corrections are
unchanged. With `mtil, Qtil_s, Ntil_s` defined from these profiles by (4.15) and (4.16),
```
(mtil, Pitil, Qtil_s, Ntil_s) = (m, Pi, Q_s, N_s)   beyond the first correction patch.
```
Before that point, profile values and each fixed number of `eta`-derivatives differ from
those of the input profile by `O(N^{-1})`.
*Construction:*
```
(C.11)   d_phi Acal = -(1/2)(a_L - a),   d_phi Bcal = (1/2) E (b_L - b_s)
```
(the unique zero-mean periodic antiderivatives, which exist by (C.1), are smooth, and
vanish identically on the boundary neighborhoods of `I`; extend by zero radially)
```
(C.12)   E_N = E exp( Acal(X,eta, N log X)/N ),   U_N = U + Bcal(X,eta, N log X)/N
(C.13)   a_N = a_L - 2 D_X Acal / N,   b_N = e^{-Acal/N} ( b_L + 2 D_X Bcal/(N E) )
(C.14)   sup_X sum_{j=0}^m ( |d_eta^j(E_N - E)| + |d_eta^j(U_N - U)| ) <= C_m N^{-1}
(C.15)   sup_X sum_{j=0}^m |d_eta^j(m_N - m)| + sup_X sum_{j=0}^m |d_eta^j(Pi_N - Pi)| <= C_m N^{-1}
(C.16)   sup_X sum_{j=0}^m |d_eta^j(p_{s,N} - p_s)| <= C_m N^{-1}
```
The exact shear identities (C.13) come from `X d_X` at fixed `phi` acting as
`D_X + N d_phi`. The **shear changes by order one** (`a_L` replaces `a`) while the
profiles, moments, and `p_s` change by `O(N^{-1})`. Moment restoration is done on the
first reserved patch, where `U = 0` and `E = K(eta) X^{-1/2 - lambda}` with `K` smooth
and positively bounded below; two `U`-bumps and three `E`-bumps, ordered as
`(M, J; I, S, C_p)`, give the block-diagonal weight table
```
rows    | weights in dX      | powers of X
M, J    | 1, H               | 0, -lambda
I,S,C_p | sqrt(2X), -E, E/X  | 1/2, -1/2 - lambda, -3/2 - lambda
```
with distinct exponents for the already-fixed `lambda > 0`; Lemma A.1 / Corollary A.3
give an invertible matrix with uniform inverse on `[-1,1]` (the bound may depend on
`lambda`; **the two `U`-weights coalesce as `lambda -> 0`, so no uniformity in that
limit is asserted**). Lemma A.2 then solves the quadratic system, with coefficients
`O_m(N^{-1})`; `N` is chosen after `lambda`, the patch scales, and all choices entering
the fixed input profile. At the right endpoint `X_rep`,
```
(C.17)   mtil(X_rep,eta) = m(X_rep,eta),   (Etil, Util) = (E, U) for X >= X_rep
```
and Lemma 4.4 propagates equality of pressure, `Q_s`, `N_s` and their parameter
derivatives for all `X >= X_rep`.

**Proposition C.3 (pp. 163-165) - the two edges and the completed profile.** The
profiles of Proposition C.2 satisfy **all** conclusions of Theorem 4.6, with the weight
```
(C.18)   y_a = log(X/X_a),   y_b = log(X_b/X),   delta = min{1, y_a, y_b},
         zeta = exp( -t_1^2/y_a^2 - 4/y_b^2 )   (X_a < X < X_b),
```
extended by zero for `X <= X_a` and `X >= X_b`. (`t_1` is the fixed activation width of
Proposition B.5 / (B.26).)
*Proof content:* stress support (part (ii)) - zero up to `X_a` by Proposition B.5 and
beyond `X_b` by Lemma A.8; nonzero in between because a vanishing point would force
`p_s = (a, -b_s)`, hence `P_c = v_s`, contradicting the strict quadratic cone test
`(v_s - 2) J_c^2 < 2(P_c - v_s)^2`. Direction and margins (part (iii)) - `n = T_0/|T_0|`
on `X_a < X < X_b`; at the inner endpoint `a = p_{1,r} > 0` and `v_s > 2 + c`, at the
outer endpoint (A.56) with `b_s = 0` and `Hcal > 0`. Inner edge factorization
`T_0 = e_a B_0` from (B.30) gives limiting direction parallel to `(a, -b_s) = a(1,t_s)`,
so `n_z - t_s n_theta = 0`, `n_theta + t_s n_z > 0`; outer edge from Proposition A.10:
```
n = (b_theta, y_b^6 b_z) / sqrt(b_theta^2 + y_b^{12} b_z^2),   b_theta(0,eta) > 0
```
so the outer limiting direction is `(1,0)` and `t_s = 0` there. Interior identities from
(4.11):
```
n_theta + t_s n_z = (P_c - v_s)/|p_s - (a,-b_s)| > 0,
n_z - t_s n_theta = J_c/|p_s - (a,-b_s)|
```
so `2 - (v_s - 2)(n_z - t_s n_theta)^2/(n_theta + t_s n_z)^2 > 0`, with value two at the
outer edge; compact positive minima give
`n_theta + t_s n_z >= kappa`, `(v_s - 2)(n_z - t_s n_theta)^2 <= (2 - kappa)(n_theta + t_s n_z)^2`
for a single `0 < kappa < 2`. Weighted estimates (part (iv)):
```
(C.19)   |T_0| >= c zeta,   |d^I T_0| <= C_I zeta delta^{-m_I}   (X_a < X < X_b, -1 <= eta <= 1)
```
proved by comparing `zeta` with `e^{-t_1^2/y_a^2}` on an inner collar ((B.31)) and with
`e^{-4/y_b^2}` on an outer collar ((A.48)-(A.51)), and by compactness in between.
Part (i): Proposition B.2 supplies smooth axis profiles with `F > 0`; Corollary B.6
preserves the analytic rectangle `[0, X_an] x [-1,1]`; pressure normalization (4.25) from
Lemma A.5 and Propositions A.7, B.8, retained by (C.17). Part (v): exact exterior
moments (4.28) from Propositions A.7 and B.8, retained by (C.17); heat exterior (4.29)
from Lemma A.6; `X_v = X_end` in `(X_a, X_b)` is the pulse-end radius before (A.10), and
beyond it `U = M = 0` (Proposition A.4 and (C.17)), so `A_X(U) = M/X = 0` and (4.7)
gives `V_0 = 0`. Part (vi): `I_pos`, `I_mean` are the third and fourth reserved intervals
listed after (A.9); they lie in `(X_a, X_v)` with `sup I_pos < inf I_mean`, are left
unchanged by Proposition C.2, and there `U = 0`,
`E = K(eta) X^{-1/2 - lambda}` with `K(eta) = c_patch (1 + eta^2)^{-1}` (Proposition
A.4), i.e. (4.30). The finite frequency `N` is fixed in Proposition C.2 before the
physical scale `q` and all later bands and correction stages, which gives the
`q`-independence asserted at the end of Theorem 4.6.

---

## 2. Dependency chain and the map into Section 4

### 2.1 Intra-appendix dependency graph

```
Lemma A.1 (moment matrix invertible)
  |-- (A.1) quantitative O(lambda^{-1}) degeneration
  |-- Corollary A.3 (five-moment Jacobian, alpha not in {-1/2,1/2,3/2})
Lemma A.2 (quadratic solve, smallness 8 beta^2 kappa d <= 1)
  |
  +--> Proposition A.4  (outer profile: A.5-A.13 shapes; A.14-A.20 moment closure;
  |        A.24-A.31 cone; needs A.1, A.2, A.3, Lemma 4.5)
  |        `-- Lemma A.5 (axis pressure datum Pi_0, A.21-A.23)  [uses only E, not Amp]
  |
  +--> Lemma A.6 (heat factor Hcal, A.32-A.38)
  |        `--> Proposition A.7 (heat replacement + 3-bump compensation on I_2,
  |                A.39-A.43; needs A.1, A.2, A.4, A.6)
  |
  +--> Lemma A.8 (backward stress representation; needs A.4, A.6, A.7 moments
  |        + regular axis from Prop B.2/Cor B.10)
  |        `-- Lemma A.9 (flat-endpoint factorization, A.47)
  |             `--> Proposition A.10 (edge factorizations A.48-A.51; needs A.6, A.8, A.9)
  |
  +--> Lemma B.1 (analytic coefficient algebra Bcal_rho, B.5-B.10)
           `--> Proposition B.2 (analytic axis profile, B.11-B.16; needs Lemma A.5
                    for Pi_0, Section B.1 axis data, Lemma B.1)
                 `--> Proposition B.3 (positive sources, exit inequality B.19)
                       `--> Lemma B.4 (reference continuation bounds B.23-B.25)
                             `--> Proposition B.5 (nonzero stress, factorization B.30-B.31;
                                     uses Lemma A.9 for the flat-factor division)
                                   |-- Corollary B.6 (reserved analytic collar)
                                   `-- Lemma B.7 (T_sh before amplitude, B.32-B.34)
                                         `--> Proposition B.8 (exact five-moment match,
                                                 B.35-B.39; needs A.1, A.2, A.3, Lemma 4.4)
                                               `--> Corollary B.10 (= B.5 + B.8 + B.6)
Lemma C.1 (loop with prescribed mean, C.1-C.10; needs Lemma 4.5 only through (4.21))
  `--> Proposition C.2 (modulation at phase N log X, C.11-C.17; needs C.1, A.1/A.3, A.2,
          Lemma 4.4, and the fixed input profile A.4 + A.7 + B.10)
        `--> Proposition C.3 (= Theorem 4.6; needs C.2, A.8, A.10, A.5, A.6, A.7,
                B.2, B.5, B.6, B.8, and Prop 4.2 for the residual identities)
```

Note the **cross-appendix cycle that is not actually circular but must be read in the
right order**: Lemma A.5 (outer) supplies `Pi_0` to Proposition B.2 (axis); Proposition
B.2 / Corollary B.10 supplies the regular axis to Lemma A.8 (outer). The paper resolves
this by fixing `Pi_0` first from `E` alone (Lemma A.5 explicitly omits the azimuthal
corrections that would depend on the amplitude) and only afterwards using the axis
construction. The order-of-choices remark (B.40), p. 157, is the audit trail for this.

### 2.2 Which Section 4 statements each appendix result supplies

Section 4 packages the appendices as three lemmas/propositions (all on pp. 34-38):

| Section 4 statement | Proved by |
|---|---|
| Lemma 4.7 (moment matrix + quadratic solve) | Lemmas A.1 and A.2 (stated verbatim on p. 34-35) |
| Lemma 4.8 (outer family, parts i-iii) | Proposition A.4 + Lemma A.5 (part i, incl. the `Pi_0` bounds); the four intervals after (A.9) (part ii); Lemma A.6 + Proposition A.7 (part iii, heat tail and its compensation on `I_2`) |
| Lemma 4.9 (backward stress + outer edge) | Lemma A.8 and Proposition A.10 |
| Proposition 4.10 (inner extension, parts i-iii) | Propositions B.2, B.3 (through `Y = 4.1` and the exit inequality `p_{s,1} + p_{s,2}^2/p_{s,1} > 2 + c_ex`); Proposition B.5 + Corollary B.6 (part i regularity, part ii first-collar factorization (4.33)); Proposition B.8 (part iii, exact moment match, with (4.34) the reference integrals) |
| Lemma 4.11 (periodic shear family) | Lemma C.1, via (C.8), (C.9), (C.10) |
| **Theorem 4.6 in full** | Proposition C.3 |

Theorem 4.6 part by part, as verified explicitly on pp. 163-165:

- **(i)** smoothness of `F = phi/C`, `U`, `Pi`, `V_0/X` on `[0,R] x [-1,1]`, `phi > 0`,
  `E = sqrt(2X) F > 0`, existence of `X_an` with `eta`-analyticity on
  `[0, X_an] x [-1,1]`, inner-collar directional margin, pressure normalization (4.25).
  *Supplied by:* Proposition B.2 (analytic axis, `Phi >= c_0 > 0`), Corollary B.6
  (analytic rectangle `[0, X_an] x [-1,1]`), Proposition C.2 (later modifications smooth
  and positivity-preserving), Lemma A.6 (`eta = +-1` smoothness of the heat profile - no
  analyticity of `Hcal` at `Z = 0` is needed), Lemma A.5 + Propositions A.7 and B.8 for
  (4.25), retained by (C.17). Stated on p. 164, last two paragraphs.
- **(ii)** exact radial pressure balance `Pi_X = E^2/(2X)` with smooth extension
  `Pi_X = F^2` at `X = 0`; exact tangential residual identities; `T_0 = 0` for
  `X <= X_a` and `X >= X_b`, nonzero on `(X_a, X_b)`; (4.13) solved for `X <= X_a`.
  *Supplied by:* Proposition B.2/B.5 (`T_0 = 0` and (4.13) up to `X_a`), Lemma A.8
  (`T_0 = 0` beyond `X_b`), Proposition C.2 (preservation), the strict quadratic cone
  test (nonvanishing in between), Proposition 4.2 (residual identities). p. 163, p. 164.
- **(iii)** positive lower bounds for `F`, `a`, `v_s - 2` on the closed annulus; smooth
  extension of `n = T_0/|T_0|` to the closed annulus; `n(X_a,eta)` parallel to
  `(a, -b_s)(X_a,eta)`; `n(X_b,eta) = (1,0)`, `b_s(X_b,eta) = 0`; the two margins (4.26).
  *Supplied by:* Proposition C.2 (interior admissible cone via Lemma C.1's `kappa_L`),
  Proposition B.5 / (B.30) (inner edge direction), Proposition A.10 / (A.48)-(A.50)
  (outer edge direction and `t_s = 0`), (A.56) (`b_s = 0`, `a > 2 + h` at the outer
  edge). pp. 163-164.
- **(iv)** the flat weight `zeta` and the bounds (4.27) `|T_0| >= c zeta`,
  `|d^alpha T_0| <= C_alpha zeta delta^{-m_alpha}`.
  *Supplied by:* (C.18) definition, (B.31) inner collar, (A.48)-(A.51) outer collar,
  compactness in between: this is (C.19). p. 164.
- **(v)** existence of `M(infinity)`, `J(infinity)`, `S(infinity)`; the four identities
  (4.28); `U = V_0 = 0` and `E = c_infinity X^{-A} Hcal(2d/X)` for `X >= X_b`; the exact
  radial heat solution `c_infinity s^{-A} Hcal(2 tau/s)`; positivity and one-sided
  smoothness of `Hcal` at `Z = 0`; existence of `X_v` with `U = V_0 = 0` on
  `[X_v, infinity)`.
  *Supplied by:* Proposition A.4 (A.8), Proposition A.7 (moment-preserving heat
  replacement), Proposition B.8 (exact inner match), Lemma A.6 (heat evolution and
  smoothness), (C.17) (retention); `X_v = X_end` the pulse-end radius, with
  `A_X(U) = M/X = 0` beyond it so (4.7) gives `V_0 = 0`. p. 164-165.
- **(vi)** two reserved intervals `I_pos`, `I_mean` in `(X_a, X_v)`,
  `sup I_pos < inf I_mean`, with `U = 0` and
  `E = c_patch (1 + eta^2)^{-1} X^{-1/2 - lambda}` (4.30), left unchanged.
  *Supplied by:* the third and fourth intervals reserved after (A.9), Proposition A.4
  (`U = 0`, `E` a positive multiple of `f X^{-1/2 - lambda}`), Proposition C.2 (both
  patches untouched). p. 165.

Downstream users of these outputs (for orientation, from the Section 4/5/7 text I
sampled): Theorem 4.6(iii) feeds Proposition 7.5, which converts the admissible cone into
positive squared wave amplitudes; (iv) feeds the weighted comparisons in Section 5;
(vi) feeds Lemma 5.2 and Lemma 8.7.

---

## 3. Load-bearing estimates, and what I could / could not verify

I hand-checked every computation below that is marked **[verified]**. Items marked
**[unverified]** I could not confirm from the text alone; items marked **[FLAG]** are
places where I judge the write-up under-justified even if I believe the claim.

### 3.1 Existence and regularity of the axis profiles (Proposition B.2)

**The mechanism.** After `Y = Lambda X`, the two stress-free profile equations (B.15)
become `2(Y Phi_YY + 2 Phi_Y) = -chi Phi + Lambda^{-1} R_1` and
`2(Y u_YY + u_Y) = -Z_*/L + Lambda^{-1} R_2`. Everything nonlinear carries `Lambda^{-1}`,
**except the order-one angular term `-chi Phi`**. The paper does not assume `chi` small.
Instead it inverts `1 + T`, `T = J_2 chi / 2`, exactly, using the fact that `J_2` raises
the minimum nonzero coefficient degree by one, whence
`||J_2 F||_rho <= 80 ||F||_rho / ((b+1)(b+2))` and
`||T^k|| <= (40 M_chi)^k / (k!(k+1)!)`. The Neumann series is therefore entire, not
merely convergent for small data. **This is the right architecture and I regard it as
the strongest single idea in Appendix B.**

- **[verified] Sign of the leading angular term.** p. 147 asserts
  `H_* xi_0/(Lambda L) = -chi`. From (B.3), `xi_0 = Lambda zeta_*` and
  `zeta_* = -L H_*/(H_*^2 + sigma_*^2)`, so
  `H_* xi_0/(Lambda L) = -H_*^2/(H_*^2 + sigma_*^2) = -chi`. Exact.
- **[verified] `-W_* > 2.8`.** With `U_* = 4 eta + j_0`, `d = 1 - eta^2`, `D = 1/2 - h`,
  `W_* = 1 - d U_{*,eta} - 2 D eta U_* = -3 + 8 h eta^2 - (1-2h) j_0 eta`, so
  `-W_* = 3 - 8 h eta^2 + (1-2h) j_0 eta >= 3 - .08 - .05 = 2.87` for `h <= 10^{-2}`,
  `j_0 <= .05`. Matches the paper exactly. This `2.8` is reused at least four times
  (B.21, B.23, p. 153, p. 155) and is genuinely load-bearing.
- **[verified] `H_*` has exactly one zero, `|eta_0| ~ j_0`, `eta_0 < 0`.**
  `H_*/d = D eta/(1-eta^2) + 4 eta + j_0` is strictly increasing on `(-1,1)` with range
  `R`, so exactly one zero; near zero `(D+4) eta_0 = -j_0 + O(j_0^2)`, giving
  `eta_0 ~ -j_0/4.5`.
- **[verified] `Z_*(eta_0) >= c j_0 P_*^2 > 0`, the pivotal positivity.** At `eta_0`,
  `H_* = 0` kills the `-H_* U_{*,eta}` term; `4 A eta_0 Pi_0(eta_0) > 0` because
  `eta_0 < 0` and `Pi_0 <= -(5/2) P_*^2 f^2 < 0` (Lemma A.5), with size
  `>= 4 A |eta_0| (5/2) P_*^2 f(eta_0)^2 ~ c j_0 P_*^2`; `-d Pi_{0,eta}(eta_0) >= 0`
  because `eta Pi_0'(eta) > 0` for `eta != 0` (Lemma A.5) and `eta_0 < 0`; the remaining
  `-A(1 - 2 eta_0 U_*) U_*` is `O(j_0)` with an absolute constant. Since `P_* > e^{T_d}`
  is fixed **before** `j_0` in (B.40), and both the good and the bad term are
  proportional to `j_0`, the comparison is `j_0`-free and holds. **This is exactly the
  place where a positivity assumption could have been smuggled in, and it is not: it is
  proved, and it is proved from Lemma A.5's two sign properties of `Pi_0`.**
- **[FLAG] The whole positivity of the axis construction rests on `Pi_0 < 0` with
  `eta Pi_0' > 0`, which Lemma A.5 proves only for `E_id,sched`, i.e. for the schedule
  with the azimuthal (E-) bumps of (A.11) **omitted**.** The reinstatement argument -
  (A.11) makes the two bumps change the total pressure increment by exactly zero, so
  `Pi_0` is unchanged - is correct in outline, but the paper's proof of the sign
  property (p. 134: "every part of `E` has the form `c(y) f(eta)^{vartheta(y)}` with `c`
  and the transition lengths independent of `eta`") is precisely the property the bumps
  can violate. The paper never verifies that the reinstated bumps do not destroy the
  *strict* inequality `eta Pi_0' > 0` pointwise; it only argues the *total integral* is
  unchanged. Since only the total is used to define `Pi_0`, I believe the claim, but
  the reader has to notice that `Pi_0` is defined from `E_id,sched` and that (A.11)
  makes the two definitions coincide. **Terse but, on my reading, correct.**
- **[verified] Why `E/sqrt(2X)` is smooth.** It is smooth *by construction, not by
  regularity theory*: `E = C^{-1} sqrt(2X) phi` (4.3) with `phi = phi_* Phi`, `phi_*`
  a function of `eta` alone (B.3) and `Phi` a convergent power series in `Y = Lambda X`
  with `Phi(0,eta) = 1` (B.12). Hence `F = phi/C = E/sqrt(2X)` is analytic in `X` near
  zero and positive for large `Lambda` by (B.11). The `sqrt(2X)` factor in `E` is the
  exact linear vanishing of `u_theta` at the axis, and (4.5) then gives a smooth
  Cartesian field. Consistency check I performed: `Pi = Pi_0 + Lambda^{-1} I(g^2 Phi^2)`
  gives `Pi_X = Lambda Pi_Y = g^2 Phi^2 = phi^2/C^2 = E^2/(2X)`, matching (4.7). The
  remaining part of (4.4), smoothness of `V_0/X`, follows from (4.7) because `U` and
  `A_X(U)` are analytic in `X` - **the paper does not spell this out** (p. 148 mentions
  only `E`, `phi`, and Cartesian smoothness), but it is immediate.
- **[FLAG] The `Y`-radius / `eta`-radius bookkeeping (p. 148).** The step
  `sum_alpha binom(alpha+beta,beta)(R/20)^alpha = (1 - R/20)^{-beta-1}`, valid for
  `R < 20`, converts the coefficient-norm bound into an honest analytic function on
  `|Y| <= R` with `4.1 < R < 20`, and the resulting `eta`-analyticity radius is
  "any sufficiently small number below `rho(1 - R/20)`". Three shrinkings are chained
  (`Omega` -> smaller neighborhood -> `rho` strictly below its boundary distance ->
  `rho(1-R/20)`) and the paper never names the final neighborhood. Theorem 4.6(i)
  demands "one common complex neighborhood of `[-1,1]`" for the scalars **and all their
  fixed radial derivatives**; Cauchy estimates on a smaller `Y` disk supply the radial
  derivatives, and the paper says so, but the reader must assemble the final
  neighborhood himself. **Compressed, not wrong; the constant `20` in the weight (B.4)
  and the target radius `4.1` are the two magic numbers here and their relation
  (`4.1 << 20`) is where all the slack lives.**
- **[verified] The choice of axis data does not require smallness of `chi`.** (B.2) makes
  `chi > .99` only on `{|Z_*| <= delta_*}`, which is a *largeness* statement, and it is
  achieved by shrinking `sigma_*` after `delta_*` - legitimate, because `eta_0` is the
  unique zero of `H_*` and `Z_*(eta_0) > 0`. The two-region argument of Proposition B.3
  (p. 149-150) then uses `chi > .99` on one region and `|Z_*| > delta_*` on the other:
  these are complementary **by (B.2)**, which is why (B.2) is stated as an implication
  and not as a global bound. Clean.
- **[verified] (B.11) and the endpoint numerics of Proposition B.3.** Since
  `f_0(z) = sum_{alpha}(-z/2)^alpha/(alpha!(alpha+1)!)`, with `t = z/2` the truncation
  `1 - t/2 + t^2/12 - t^3/144` is exact through cubic order (`1!2! = 2`, `2!3! = 12`,
  `3!4! = 144`), the tail alternates with decreasing magnitude for `t <= 2.05`, and the
  cubic is decreasing there (its derivative `-1/2 + t/6 - t^2/48` has negative
  discriminant). Its value at `t = 2.05` is `305719/1152000 = .26538...`, matching the
  paper digit for digit. Likewise
  `f_0(z) + z f_0'(z) = sum_alpha (-t)^alpha/(alpha!)^2 = 1 - t + t^2/4 - t^3/36 + t^4/576 - ...`
  because `(1+alpha)/(alpha!(alpha+1)!) = 1/(alpha!)^2`; at `t = 1.98` the quartic
  truncation equals `-75535511/400000000 = -.18883878`, again matching exactly, and the
  next term `-t^5/14400` is negative so the truncation is a genuine upper bound. Then
  `-2 z f_0'/f_0 = 2 - 2(f_0 + z f_0')/f_0 > 2 + 2(.18) = 2.36` using `0 < f_0 <= 1`.
  With `p_1 = a = -2 Y Phi_Y/Phi` and (B.13) this gives `p_1 > 2.3`. **Every number in
  this argument checks out.** This is the best-verified page in the appendices.
- **[verified] The complementary alternative.** On `chi <= .99`, (B.2) forces
  `|Z_*| > delta_*`, hence `|n_s| >= delta_*/2` by (B.18); then
  `|p_2| = C sqrt(2/Lambda) |n_s|/(phi_* Phi)` grows linearly in `C` with `Lambda`
  fixed, while `p_1 <= C_1` is `C`-independent, so `p_2^2/p_1 > 2.3`. Both branches give
  (B.19) with `c_ex = .2`. **Sound; note this is where the amplitude `C`, not `Lambda`,
  does the work, which is why `C_0(Lambda)` must be increased after `Lambda`.**
- **[FLAG] `C_0(Lambda)` is exponentially large in `Lambda`.** p. 148 states plainly:
  "It is essential here that the bound on `g` holds on a complex neighborhood, since
  `phi_*` itself can grow exponentially with `Lambda`." Combined with
  `X_R = X_i (C P_*)^{10}` and `P_* > e^{T_d} = e^{e^{M_d}+10}`, the construction is
  entirely non-quantitative in absolute terms. This is not an error - the theorem only
  claims existence of fixed finite parameters - but any attempt to make the paper
  numerical or computer-assisted at this stage would fail, and no reader can form an
  intuition for the sizes involved.

### 3.2 Joining and moment matching (Lemmas A.1, A.2, Corollary A.3, Propositions A.7, B.8, C.2)

- **[verified] Lemma A.1.** The Rolle/Descartes induction is correct: dividing
  `sum c_i x^{alpha_i}` by `x^{alpha_1}` and differentiating gives a combination of
  `m-1` distinct powers, so `m` positive zeros would give `m-1` for the derivative,
  contradicting the induction hypothesis (at most `m-2`). Multilinearity in the columns
  gives `det B = int det[x_j^{alpha_i}] prod beta_j`, whose integrand has one sign and
  is nonzero on a positive-measure set because the `beta_j` are nonzero and supported in
  *ordered* intervals. **Correct and standard.**
- **[verified] Corollary A.3's exponent hypothesis is genuinely necessary and is checked
  at every application.** The row operations `J -> J - u_c I`, `S -> S - 2 u_c M` split
  the `5x5` Jacobian into a `2x2` `U` block with weights `1, x^{alpha+1/2}` and a `3x3`
  `E` block with weights `x^{1/2}, x^alpha, x^{alpha-1}`. Distinctness within blocks is
  exactly `alpha != -1/2, 1/2, 3/2`. Applications:
  `alpha = 1/10` (axis correction, Prop B.8) - fine;
  `alpha = -1/2 - lambda` (intermediate interval, Prop C.2) - fine for `lambda > 0`;
  Prop A.7 uses only three `E` bumps with weights `x^{-3/2-lambda}, x^{-1/2-lambda},
  x^{1/2}` - distinct for `lambda > 0`. **All checked.**
- **[FLAG - the one real fragility in the moment machinery] The `O(lambda^{-1})` loss.**
  (A.1) and p. 128 record that in the intermediate power-law interval the two `U`-block
  exponents are `0` and `-lambda`, which coalesce as `lambda -> 0`, giving
  `||B^{-1}|| = O(lambda^{-1})`. Proposition C.2 (p. 162) states honestly that "the two
  `U`-weights coalesce when `lambda -> 0`, so no uniformity in that limit is asserted."
  The paper's answer is that `lambda` is fixed very early in (B.40), so `lambda^{-1}` is
  a fixed constant, and the discrepancies are made small **afterwards**: `O(X_K^{-1})`
  by large `X_R` in Prop A.7, `O(N^{-1})` by large `N` in Prop C.2, and
  `C_pre e^{-c/lambda}(1+|Amp|)` (A.15) - exponentially small in `1/lambda` - for the
  pulse corrections. That last one is why the pulse must end at `11/lambda`, at least
  `2/lambda - 3` before the first bump center at `13/lambda - 3`. **The logic is sound.
  But the derivation of `O(lambda^{-1})` for the `1, x^{-lambda}` pair is one sentence
  ("expanding `x^{-lambda} = 1 - lambda log x + O(lambda^2)`: the two bump averages of
  `log x` have strictly separated values") with no display of the remainder control, and
  it is the hinge on which three separate correction stages hang.** I believe it (the
  bump supports are fixed in the logarithmic coordinate, so `log x` ranges over a fixed
  compact set), but it deserves a proof, not a clause.
- **[verified] Lemma A.2's quadratic solve.** `c -> B^{-1}(d - Q(c,c))` maps the closed
  ball of radius `r = 2 beta_0 d_0` to itself (`r/2 + beta_0 kappa_0 r^2 <= r` iff
  `2 beta_0 kappa_0 r <= 1`, i.e. `4 beta_0^2 kappa_0 d_0 <= 1`, implied by the stated
  `8 beta_0^2 kappa_0 d_0 <= 1`) with Lipschitz constant `2 beta_0 kappa_0 r <= 1/2`.
  The same inequality makes `B + D_c Q(c,c)` invertible, giving smoothness by the
  pointwise implicit function theorem. **Correct, and the factor 8 rather than 4 buys the
  derivative bound.**
- **[FLAG] The "finitely many derivative orders need smallness" device.** Lemma A.2's
  `C^k` conclusion needs `8 beta_k^2 kappa_k ||d||_{C^k} <= 1` for the *specific* `k`
  chosen in advance; higher orders are only asserted finite. This device is used at
  least four times (p. 133 "Only finitely many orders of `eta`-derivatives must satisfy
  prescribed smallness bounds"; p. 140 "Only the finitely many derivative bounds used to
  preserve the cone must meet prescribed smallness thresholds"; p. 157 Remark B.9;
  p. 162 "no simultaneous smallness condition for all derivative orders is needed").
  Remark B.9 (p. 157) is the paper's answer: "Fix any finite list of `eta`-derivative
  orders whose bounds are needed to preserve the profile inequalities and solve the
  finite moment equations, including the extra parameter derivative in (B.35)." **This
  is the correct fix and it is stated. What is NOT done anywhere is an audit that the
  list is finite and closed: every use of Lemma 4.4(ii) (4.17) costs one extra
  `eta`-derivative (`||Delta(Q_s,N_s,p_s)||_k <= C_k(||Delta(U,E)||_{k+1} + ||Delta m||_{k+1})`),
  and the appendices chain several such comparisons (Prop A.7 -> Lemma A.8 ->
  Prop C.2 -> Prop C.3). If the required order grows by one at each link, the list is
  still finite, but nobody counts the links.** I judge this the most likely place for a
  hidden circularity, and I could not close the audit from the text.
- **[verified] Proposition A.7's normalization arithmetic.** The three discrepancies
  (A.40)-(A.42), divided by the natural scales (A.43), all come out `O_m(X_K^{-1})`:
  `(e_K^2/X_K)/e_K^2 = X_K^{-1}`; `e_K^2/(X_K e_K^2) = X_K^{-1}`;
  `e_K X_K^{1/2}/(X_K^{3/2} e_K) = X_K^{-1}`. The convergence exponents also check:
  `int_1^inf x^{-2A-2} = 1/(2+2h)`, `int_1^inf x^{-2A-1} = 1/(1+2h)`,
  `int_1^inf x^{-A-1/2} = 1/h` - the last requires `h > 0` **fixed before `X_R`**, and
  the paper says exactly that. The underlying profile estimate
  `|(X d_X)^j d_eta^m (E - E_cl)| <= C_{j,m} e_K X_K^{-1} x^{-A-1}` follows from (A.38)
  `Hcal(2d/X) = 1 - 2h(1+h)d/X + O(X^{-2})`. **Fully consistent.**
- **[verified] Proposition B.8's scale elimination.** (B.35) rewrites `W`, `Q_s`, `N_s`
  entirely in the normalized moments, so `X_R` disappears from those formulas; then
  `P_c = X_R x G / L` with `G >= 3 Q_min/7` bounded below by outer-profile data alone, so
  `P_c > 2` is obtained by increasing `X_R` **after** the tolerance (B.36) is fixed.
  `x_sep = e^{T_sh}/(C P_*)^{10} < e^{-8}` is then arranged by large `C` after `T_sh`.
  The estimate `v_s <= .9 + .01/.7 < 1` from (B.36) is correct arithmetic
  (`v_s = a(1+t_s^2) <= .9(1 + (.1/.9)^2)`, and the paper's `.9 + .01/.7` is a valid
  crude bound). **The order of choices here is the crux and Lemma B.7 is exactly the
  lemma that makes it non-circular: `T_sh` is bounded using `Bcal_k`, which bounds
  `ell_i = log(C E(X_i,eta))` - and `C E = sqrt(2X) phi` is `C`-free, so `Bcal_k` is
  legitimately `C`-independent.** [verified]
- **[verified] Proposition C.2's exact shear identities (C.13).** Because
  `X d_X` applied to `Acal(X,eta,N log X)` equals `D_X Acal + N d_phi Acal`, and
  `d_phi Acal = -(1/2)(a_L - a)`,
  `a_N = 1 - 2 D_X log E_N = a + (a_L - a) - 2 D_X Acal/N = a_L - 2 D_X Acal/N`,
  and similarly `b_N = e^{-Acal/N}(b_L + 2 D_X Bcal/(N E))`. **Both identities are exact
  and I reproduced them. This is the heart of Proposition C.2 and it is airtight.**
  The `O(N^{-1})` profile estimate (C.14) uses that the phase `N log X` carries no `eta`,
  so `eta`-differentiation never produces a factor of `N`; and (C.15)-(C.16) then follow
  because (4.16) expresses `Q_s`, `N_s` through radial *integrals* and *parameter*
  derivatives, never radial derivatives.
- **[FLAG] Radial derivatives of the modulated profile grow like `N^{r-1}`.** p. 161
  admits this: "The derivative `X d_X(E_N - E)` can have order one. In general, for fixed
  `r >= 1, m >= 0`, radial differentiation of (C.12) gives a bound `C_{r,m} N^{r-1}`."
  p. 163 then fixes one `N` and allows all later constants to depend on it. **Logically
  fine, and Theorem 4.6's closing paragraph does permit dependence on "the fixed profile
  choices". But every downstream estimate in Sections 5-8 that uses a radial derivative
  of the leading profile now carries a factor `N^{r-1}` in its constant. That
  bookkeeping lives outside Appendices A-C and I could not verify it. If any later
  argument needs `N -> infinity`, the construction breaks. p. 163 pre-empts this ("This
  choice is made before any limit in the physical scale `q` or any later dyadic band or
  correction stage"), so the intent is clear, but the reader of Section 5 onward must
  check it.** This is my single largest residual concern about Appendix C.

### 3.3 The heat-exterior replacement preserving the axis pressure (Lemmas A.6, A.8, Proposition A.7)

This is the part I was able to verify most completely, and it holds up.

- **[verified] (A.34), (A.35) and the ODE (A.37) by direct computation.** Writing
  `Hcal(Z) = (1/Gamma(1+h)) int_0^inf e^{-v} v^h (1+Zv)^{-h} dv` (note
  `a_K - 1 = h`, so (A.32)'s `v^{a_K-1} = v^h`), dominated differentiation gives (A.34);
  setting `Z = 0` gives `Hcal^{(m)}(0) = (-1)^m (h)_m Gamma(1+h+m)/Gamma(1+h) = (-1)^m (h)_m (1+h)_m`,
  which is (A.35). For (A.37) I integrated `h d_v[e^{-v} v^{1+h}(1+Zv)^{-1-h}]` over
  `(0,infinity)` (both boundary terms vanish: `v^{1+h} -> 0` at zero, `e^{-v} -> 0` at
  infinity) and reduced the three resulting integrals using the two algebraic identities
  `Z v^{h+1}(1+Zv)^{-1-h} = v^h(1+Zv)^{-h} - v^h(1+Zv)^{-1-h}` and
  `Z v^{h+2}(1+Zv)^{-2-h} = v^{h+1}(1+Zv)^{-1-h} - v^{h+1}(1+Zv)^{-2-h}`. The result is
  exactly `Z^2 Hcal'' + [1 + 2(1+h)Z] Hcal' + h(1+h) Hcal = 0`, i.e. (A.37) with
  `a_K = 1+h`, `a_K(a_K-1) = h(1+h)`. **Confirmed symbolically.**
- **[verified] The exponent identities that make (A.33) a heat solution.** With
  `K = c_infinity s^{-A} Hcal(2 tau/s)`, `s = r^2/2`, `Z = 2 tau/s`:
  `d_t K = -2 c_infinity s^{-A-1} Hcal'` and
  `(d_rr + r^{-1} d_r - r^{-2})K = 2 c_infinity s^{-A-1}[Z^2 Hcal'' + (2A+1) Z Hcal' + (A^2 - 1/4) Hcal]`.
  Since `A = 1/2 + h`: `2A+1 = 2(1+h) = 2 a_K` and
  `A^2 - 1/4 = h + h^2 = h(1+h) = a_K(a_K - 1)`. Equating gives precisely (A.37).
  **The choice `a_K = 1 + h` is forced by `A = 1/2 + h`, and it works exactly. This is
  the cleanest identity in the appendices.**
- **[verified] `-Z Hcal'/Hcal` is a probabilistic average, hence in `[0,h)`.**
  `-Z Hcal' = (1/Gamma) int e^{-v} v^h (1+Zv)^{-h} [h Z v/(1+Zv)] dv`, so
  `-Z Hcal'/Hcal = <h Z v/(1+Zv)>` under the positive density proportional to (A.32)'s
  integrand: nonnegative, and `< h` because `Zv/(1+Zv) < 1`. Consequently
  `r K_r/(2K) = -A - Z Hcal'/Hcal` lies in `[-1/2 - h, -1/2)`, giving both `K_r < 0` and
  the strict bound `< -1/2` of Lemma A.6. **All verified.** I also verified
  `r K_r = 2K(-A - Z Hcal'/Hcal)` directly from `K = c s^{-A}Hcal(2tau/s)`.
- **[verified] (A.56) `a > 2 + h` at the outer edge.** `a = 2 + 2h + 2 Z Hcal'/Hcal - 2 f_o'/f_o`.
  The paper needs both `f_o'/f_o < h/4` and `-Z Hcal'/Hcal < h/4`. The first is
  *arranged by construction*: `f_o = 1 - c_o h psi_o` with `c_o` chosen so that
  `0 <= f_o'/f_o < h/4`, and `c_o` is `h`-independent because
  `f_o'/f_o = -c_o h psi_o'/(1 - c_o h psi_o)` and `psi_o' <= 0`, so it suffices that
  `c_o ||psi_o'||_infinity < 1/4`. The second needs `Z = 2d/X` small, which holds for
  large `X_R`: `-Z Hcal'/Hcal = <h Z v/(1+Zv)> <= h Z <v> = h Z (1+h) + o(Z) -> 0`.
  **The paper asserts the second bound without this line; it is immediate but not shown.**
  [FLAG - minor]
- **[verified] The four backward representations in Lemma A.8.** In the tail `U = 0`
  gives `M = J = 0` and `A_X(U) = 0`, so `W = 1`, and (4.16) collapses exactly to
  `Q_s = -1 + ((1-h)I - D eta I_eta)/(XH)` and
  `N_s = (4h eta S - d S_eta)/X + 4A eta Pi - d Pi_eta`, as displayed. Using
  `int_0^infinity (H - H_pow) dX = 0` and `H_pow = sqrt 2 c_infinity X^{-h}`:
  `int_0^X H_pow dx = sqrt 2 c_infinity X^{1-h}/(1-h) = X H_pow(X)/(1-h)` (this is where
  `h < 1` enters, as the paper notes), so
  `I(X) = X H_pow/(1-h) - int_X^infinity (H - H_pow) dx`. Using `S(infinity) = 0` and
  `U = 0`: `S(X) = (1/2) int_X^infinity E^2 dx`. **All four verified.** The paper spends
  one sentence on this ("The others follow directly from the exact moments and the
  integral identities for `Q_s`, `N_s` (4.16)") - correct, but it is the algebraic
  keystone of the appendix and deserved a display.
- **[verified] The three total-integral identities and the boundary terms.**
  `int_0^infinity r^2 (u_theta - K_pow) dr = q^{3/2-A} int_0^infinity (H - H_pow) dX` -
  I checked the change of variables (`dX = r dr/q`, `H = sqrt(2X)E`, `u_theta = q^{-A}E`).
  `int_0^infinity r(u_z^2 + p) dr = q^{1-2A} int_0^infinity (U^2 - E^2/2) dX` - this
  uses `int_0^infinity Pi dX = -int_0^infinity E^2/2 dX`, which is integration by parts
  in `Pi_X = E^2/(2X)` and needs `X Pi -> 0`; since `E ~ c_infinity X^{-A}` gives
  `Pi ~ -c X^{-1-2h}` and `X Pi ~ -c X^{-2h} -> 0`, this holds **because `h > 0`**. The
  radial-viscosity boundary term is
  `int_0^infinity r^2 (d_rr + r^{-1}d_r - r^{-2}) u_theta dr = [r^2 (u_theta)_r - r u_theta]_0^infinity`
  (I verified `d_r(r^2 u' - r u) = r^2 u'' + r u' - u`), which vanishes at infinity by
  (A.44)'s third estimate `r^2 K_r - r K = O_h(r^{-2h})` and at zero by axis regularity.
  The exact pressure identity
  `int_0^infinity r p dr = -(1/2) int_0^infinity r' u_theta(r')^2 dr'` follows from
  Fubini. **All verified.**
- **[verified] The exterior `N_s` bound on p. 137.** With `l <= -3h/4` beyond the
  exterior transition, `E(v)^2 <= E(y)^2 e^{-(1+3h/2)(v-y)}`, so
  `h int_y^infinity e^{v-y} E(v)^2 dv <= h E(y)^2 int_0^infinity e^{-3hs/2} ds = (2/3) E(y)^2`
  and `int_y^infinity E(v)^2 dv <= E(y)^2/(1+3h/2)`; with `A = 1/2+h` these give
  `N_s = 2 eta int_y^infinity (h e^{v-y} - A) E(v)^2 dv = O(E(y)^2)` **uniformly in `h`**,
  which is what the paper claims. Note that the `h`-uniformity comes from the exact
  cancellation `h * (2/(3h)) = 2/3`; this is a nice and non-obvious point.
- **[verified] Proposition A.10's factorizations via Lemma A.9.** Lemma A.9's
  substitution `u = delta/(1+delta^2 v)^{1/2}` gives Jacobian
  `-(1/2) delta^3 (1 + delta^2 v)^{-3/2}` and `c/u^2 = c/delta^2 + c v`, hence (A.47) with
  the displayed `B`; setting `delta = 0` gives `B(0,eta) = b(0,eta)/c`. **Correct.**
  Applying it with `j = 3, c = 4` to the two integral terms of (A.54) makes them vanish
  to order `delta^3` relative to the boundary term `K f_r`, which is what produces the
  asymmetry `delta^{-3}` vs `delta^{+3}` in (A.48)-(A.49) and hence
  `T_{0,z}/T_{0,theta} = delta^6 b_z/b_theta -> 0`. The three-term nonnegativity in
  (A.54) (`K > 0`, `K_r < 0`, `f_o' >= 0`) is exactly Lemma A.6 plus the construction of
  `f_o`. **Sound.**
  **[FLAG]** The claimed exponent `q^{A+1/2} q^{1/2 - D - 2A} = q^{1 - D - A} = 1` on
  p. 143 uses `A + D = 1`, which is (4.1); I verified it. But the reader must supply
  `A + D = 1` from `A = 1/2 + h`, `D = 1/2 - h` - it is used silently at least three
  times in Appendix A.

### 3.4 The cone realization (Lemma C.1, Propositions C.2, C.3)

- **[verified] Lemma C.1's variance identity and monotonicity.** With
  `t - t_s = d_0 (e^{mu p_2 sin theta'}/M_e(mu p_2) - 1)/p_2`, expanding the square gives
  `<(t-t_s)^2> = (d_0^2/p_2^2)[M_e(2 mu p_2)/M_e(mu p_2)^2 - 1]`, exactly the displayed
  `V`. `<t> = t_s` is immediate. Monotonicity: `g = log M_e` has
  `g''(z) = Var_tilted(sin theta') > 0`, and
  `d_mu[g(2 mu p) - 2 g(mu p)] = 2 p[g'(2 mu p) - g'(mu p)] > 0` for `p != 0, mu > 0`
  (both factors flip sign together when `p < 0`). **Verified.**
- **[verified] `V -> infinity` and the constant `3`.** `M_e(z) ≍ e^{|z|}/sqrt{|z|}` gives
  `M_e(2z)/M_e(z)^2 ≍ sqrt{|z|/2}`, so `V ≍ (d_0^2/p^2) sqrt{mu|p|/2} -> infinity`; at
  `p = 0`, `V = d_0^2 mu^2/2 -> infinity`. The target in (C.7) is `3/min a` because
  `rho = zeta_L(v_s)^2 (v_* - v_s) < v_* = 2 + delta_L/2 < 3`, using `delta_L < 1`, and
  (C.10) demands `V = rho/a <= 3/min a`. **Verified - and the `delta_L < 1` in the
  statement is doing real work here, which the paper does not point out.**
- **[verified] `P_c(t) > 2` along the whole loop.** `p_2 (t - t_s) = d_0(e^{...}/M_e - 1) >= -d_0`
  since the exponential ratio is nonnegative, so
  `P_c(t) = P_c(t_s) + p_2(t-t_s) >= P_c(t_s) - d_0 > 2` by the choice
  `d_0 < (1/2) min(P_c(t_s) - 2)`. **Verified; note this bound is one-sided by design
  and uses only positivity of an exponential.**
- **[verified] The lift has period one and the prescribed means.**
  `int_0^{2pi} (dphi/dtheta') dtheta' = (a/v)<1 + t^2> = (a(1+t_s^2) + aV)/v = (v_s + rho)/v = 1`
  by (C.10); and
  `int_0^1 (a_L,-b_L) dphi = (a/2pi) int_0^{2pi} (1,t) dtheta' = (a, a t_s) = (a,-b_s)`.
  **Both verified. This is a genuinely elegant construction: the whole point is that the
  cone requires `v_s > 2`, `v_s` is determined by the *variance* of the shear ratio, and
  a mean-preserving spread raises the variance without moving the mean.**
- **[FLAG - unstated but true step] Smoothness of the solution of (C.10) through `rho = 0`.**
  The paper says "the expansion of (C.4) gives, uniformly for bounded `p`,
  `V(mu,p) = (1/2) d_0^2 mu^2 + O(mu^4 p^2)`. The signed square root of `V` is smooth and
  odd near `mu = 0`, with derivative `d_0/sqrt 2`." This requires the **cubic** term of
  `g(2z) - 2g(z)` to vanish, which it does because `<sin^3 theta'> = 0` makes
  `g(z) = z^2/4 + O(z^4)` with no `z^3` term, whence
  `g(2z) - 2g(z) = z^2/2 + O(z^4)`. **The paper never mentions the vanishing of the
  cubic, and without it the signed square root would not be smooth at `mu = 0`.** I
  verified it; but this is exactly the kind of step that should be displayed.
- **[verified] (C.2) follows from Lemma 4.5.** Along the loop `c_L = P_c(t) > 2`,
  `v_L = v`, `j_L = J_c(t)`, and `2 < v < Ucal(P_c(t),J_c(t))` in all three cutoff
  regimes (`v = v_*` where `v_s <= 2 + delta_L/8`, `2 < v_s <= v <= v_*` in the
  transition, `v = v_s > 2` where the cutoff vanishes; and `v_* = 2 + delta_L/2 < 2 + delta_L < Ucal`
  by (C.8)). Lemma 4.5 then gives `P_c > v` and the quadratic inequality, and compactness
  gives `kappa_L`. **Verified.**
- **[FLAG - coverage claim not displayed] `I = [X_-, X_+]` must contain every point where
  the admissible cone can fail.** p. 158 asserts this. Its justification is distributed
  across four separate results: Corollary B.6 (admissible on the first collar),
  Proposition A.4 (admissible throughout the intermediate power-law interval, pulse,
  interpolation, exterior transition, and as soon as `l < 0`), Proposition A.7 (preserves
  it from the compensation patch through `y = .5`), and Proposition A.10 (admissible on
  `.5 <= y < 3`). **The paper nowhere states, as a single claim, that these regions cover
  `(X_a, X_b) \ I`, and there is one visible seam: between the end of the axis collar and
  the start of the region where Proposition A.4 asserts the admissible (rather than
  merely relaxed) condition. Proposition A.4's phrase "as soon as `l < 0` on the
  transition into that constant-slope interval" is the boundary of that seam, and I could
  not confirm from the text that `I`'s right endpoint `X_+` is inside that region while
  its left endpoint `X_-` is inside the collar. This is the coverage gap I would raise
  first with the authors.** [unverified]
- **[FLAG - ambiguous phrasing] Ordering of `X_an` and `X_-`.** p. 158 says the analytic
  rectangle's "radial endpoint" is chosen "in the shorter inner collar", denoted `X_an`,
  while `I` starts "inside the first collar". Proposition C.2 then claims "the rectangle
  near the axis described above" is unchanged, which requires `X_an < X_-` (with the
  `delta_d` margin of (C.3) as slack). Theorem 4.6(i) separately calls `[X_a, X_an]` "the
  inner collar" with the directional margin. **The two uses of "collar" are not the same
  interval and the ordering `X_a < X_an < X_-` is never displayed.** I believe it is
  intended and consistent, but it costs the reader a re-read.
- **[verified] Proposition C.3's nonvanishing argument.** If `T_0 = F(p_s - s)` vanished
  at an interior point with `F > 0`, then `p_s = (a, -b_s)`, hence
  `P_c = p_{s,1} + t_s p_{s,2} = a + t_s(-(-b_s))`... concretely `p_s = s` gives
  `P_c = a(1 + t_s^2) = v_s` and `J_c = 0`, so
  `2(P_c - v_s)^2 - (v_s - 2)J_c^2 = 0`, contradicting the strict quadratic test of
  Lemma 4.5 / (C.2). **Verified.**
- **[verified] The interior cone identities used for (4.26).**
  `n_theta + t_s n_z = (P_c - v_s)/|p_s - (a,-b_s)|` and
  `n_z - t_s n_theta = J_c/|p_s - (a,-b_s)|`: with `n = (p_s - s)/|p_s - s|` and
  `s = a(1,t_s)`, the projections onto `(1,t_s)` and `(-t_s,1)` give
  `(p_{s,1} - a) + t_s(p_{s,2} - a t_s) = P_c - a(1+t_s^2) = P_c - v_s` and
  `-t_s(p_{s,1}-a) + (p_{s,2} - a t_s) = J_c - a t_s + a t_s = J_c`. **Verified.** Hence
  `2 - (v_s-2)(n_z - t_s n_theta)^2/(n_theta + t_s n_z)^2 > 0` is exactly (4.22) and has
  value `2` at the outer edge, where `J_c`-type numerator vanishes to order `delta^6`.
- **[verified] The weight (C.18) matches the two edge factorizations.** With
  `y_b = log(X_b/X)` and `X_b = e^3 X_tail`, Proposition A.10's `delta = 3 - y` (where
  `y = log(X/X_tail)`) equals `log(X_b/X) = y_b`. So (A.51)'s
  `|T_0| >= c e^{-4/delta^2} delta^{-3}` dominates `zeta = e^{-t_1^2/y_a^2 - 4/y_b^2}` on
  the outer collar (the extra `delta^{-3} >= 1` and the bounded `e^{-t_1^2/y_a^2}` are
  harmless), and (B.31)'s `|T_0| >= c e^{-t_1^2/y_a^2}` does the same on the inner collar.
  The upper bounds match `C_I zeta delta^{-m_I}` in the same way. **Verified.**

### 3.5 ODE existence claims: does any rely on an unproved smallness or positivity?

My conclusion: **no ODE *existence* claim in Appendices A-C rests on an unproved
smallness or positivity hypothesis.** In detail:

1. **The only nonlinear existence result is Proposition B.2**, and it is a contraction
   argument in `Bcal_rho` whose hypothesis is *largeness* (`Lambda >= Lambda_0`,
   `C >= C_0(Lambda)`), not smallness of data. Crucially, the order-one term `-chi Phi`
   is **not** treated perturbatively: `1 + T` is inverted by a Neumann series that
   converges for every `M_chi` because of the `1/(k!(k+1)!)` gain. The positivity needed
   to divide by `phi` in (B.15) is a *conclusion* ((B.11) plus (B.13)), not an
   assumption. **[verified]**
2. **All other ODEs are first-order linear with explicit integrating factors** ((4.9),
   (4.10), (A.16), (B.25), (A.9)), so existence is never at issue. What is at issue is
   positivity of `Q_s` (needed to define `w = N_s/(E Q_s)` and to run the cone test), and
   every instance is backed by an explicit source bound or a barrier:
   - (A.25) on the reference inner interval. **[verified by hand]** I reproduced this
     formula from scratch: on (A.7), `l = 3/5` is constant so
     `Q_s = S_q/(1 + 3/5) = S_q (5/8)`; and `S_q = -Wl - h(1 - 2 eta U) - H_c (log E)_eta`
     with `A_X(U) = 4 eta`, `W = -3 + 8 h eta^2 = -(4L - 1)`, `H_c = (D + 4d) eta`,
     `(log E)_eta = -J_0'`, giving exactly
     `Q_s = [(3/5)(4L-1) - h(1 - 8 eta^2) + (D + 4d) eta J_0']/(8/5)`. Positivity:
     `(3/5)(4L-1) >= (3/5)(3 - 8h)`, `|h(1-8eta^2)| <= 7h`,
     `(D+4d) eta J_0' = 2 eta^2 (D + 4d)/(1+eta^2) >= 0`, so the numerator is
     `>= 9/5 - 24h/5 - 7h > 0` for `h <= 10^{-2}`. **The identity `4L - 1 = -W` on the
     reference profile is a small miracle that the paper does not comment on.**
   - (A.26) on the axial transition: `Q_s >= c(eta^2 + e^{-y})`. Here the source is
     `S_q >= (1/2 - h)eta^2 - h`, which is **negative near `eta = 0`**, so the positive
     lower bound comes from the homogeneous decay `e^{-y}` off a positive initial value
     over an interval of length `T_d`. **This does require `h << e^{-T_d}`, which is
     exactly the smallness built into the parameter order (A.6).** So the assumption is
     stated, not hidden. **[verified in structure; the constant `c` is not exhibited]**
   - (A.27) on the intermediate interval, (A.29) on the pulse, `Q_s >= c lambda` during
     profile interpolation, `Q_s ≍ h` and (A.16) on the terminal interval: each is
     asserted with a variation-of-constants sketch. **[unverified in detail]**
   - (B.23)/(B.25): `S_{q,r} >= .94 L Lambda chi + 2.4` and the comparison
     `p_{1,r}(y) >= p_{1,r}(0)e^{-y} + 1.88 chi(e^y - e^{-y})`. **[unverified]**
   - The barrier at the end of Proposition B.5 (p. 153): at a hypothetical downward
     crossing `p_1 = 2` on the final constant-slope interval (`a = .8`, `l = .6`),
     `D_X p_1 = X S_q/L - .6 p_1 > X - 1.2 > 0` since `S_q > 1`, `L <= 1`, and
     `X >= X_b = 100`. **[verified]** The input `S_q > 1` comes from
     `-W l > .6 * 2.8 = 1.68` plus the nonnegative gradient term (B.21) minus `O(h)`.
     **[verified]** The closeness of `W` to `W_*` uses Lemma B.4's neat observation that
     radial averaging is a contraction toward a constant. **[verified]**
3. **The one construction that is "run the ODE until it hits a target"** is the terminal
   transition of Section A.2: "retain `l = -h` until the solution of `Q' + (1+l)Q = -l-h`,
   initially `Q = (lambda-h)/(1-lambda)`, has decreased to `Q_p`". With `l = -h` the
   source is `-l - h = 0`, so `Q` decays like `e^{-(1-h)y}` from a value comparable to
   `log(1/h)` (built up during the `l = -1` interval, where `Q_s' = 1-h`) down to
   `Q_p ≍ rho_o ≍ h`. Since `h << log(1/h)`, the target is attained in finite
   logarithmic length, and the paper says so (p. 132). **[verified in structure]** I also
   checked that `c_o` in `rho_o = c_o h` is `h`-independent, so `Q_p ≍ h` is legitimate.
4. **The scalar root-finding for the pulse amplitude (A.19)-(A.20)** is not an ODE but is
   the other place where existence could fail. It is proved by an explicit sign change
   plus a derivative lower bound, and **all three numbers check out [verified]**:
   `K_b <= int_0^infinity xi^2 e^{-2xi} dxi = 1/4` and `K_b > .20` (I computed
   `e^{-.04} int_0^{8.98} e^{-2u}u^2 du ~ .240`); at `Amp = .9`,
   `.81 K_b - (1-e^{-26})/4 <= .81/4 - .25 = -.0475 <= -.047`; at `Amp = 1.2`,
   `1.44 K_b - .25 >= 1.44(.20) - .25 = .038`; the derivative
   `2 Amp K_b >= 2(.9)(.20) = .36`. The perturbation `||Ecal||_{C^1} <= C_pre lambda(1+log(1/lambda))`
   is `o(.36)` for small `lambda`, so the root exists, is unique in the bracket, and is
   smooth with `|d_eta Amp| <= C_pre lambda(1+log(1/lambda))` (A.20). **Well done.**
5. **The pulse cone numbers `< .74` and `< 1.68` (p. 136)** [verified]: with
   `b_s = -R_b + O(lambda)`, `w = 2 R_b - C_d d_{xi_b}R_b + o(1)`, `a = 2 + 2 lambda`,
   `d_{xi_b}R_b <= 1.2`, `C_d <= 2 + o(1)`,
   `b_s w <= -2R^2 + 2.41 R`, whose max is `2.41^2/8 = .7260 < .74`; and
   `2 b_s w + b_s^2/a + (a-2)w^2 <= -4R^2 + 4.82R + R^2/2 = -3.5R^2 + 4.82R`, whose max
   is `4.82^2/14 = 1.6591 < 1.68`. I also verified the two inputs:
   `d_{xi_b}R_b = Amp[phi_b'(1-sigma) - phi_b sigma'] <= Amp <= 1.2` (using
   `phi_b' = sigma(xi_b/.02) <= 1` and `sigma' >= 0`), which is **exactly why the
   amplitude bracket `[.9,1.2]` matters**; and
   `C_d = lambda(D+c_eta)(1-lambda)/(beta^2(lambda-h+c_eta))` is decreasing in `c_eta`
   (since `lambda - h - D < 0`) so is maximized at `c_eta = 0`, where it is
   `~ D/beta^2 = (1/2)/(1/4) = 2` given `h/lambda -> 0`. **All verified.**
6. **[verified] Two further construction details worth recording, because they show the
   design is deliberate rather than lucky.** (a) The `log 2` in (A.10) is there so that
   `log 2 - J_0 = log 2 - log(1+eta^2) >= 0` on `|eta| <= 1`, which makes
   `l = -lambda + vartheta_f'(log 2 - J_0) <= -lambda` **one-sided** - i.e. the
   azimuthal amplitude cannot increase during profile interpolation, exactly as the paper
   says ("the amplitude cannot increase"). (b) The exponent `30` in (A.14)'s
   `e_b <= C_pre lambda^{30}` is `(1/2 + lambda) * T_w / log(1/lambda) = (1/2)(60)`,
   because `D_X log E = l - 1/2 = -lambda - 1/2` on the intermediate interval of
   logarithmic length `T_w = 60 log(1/lambda)`. **So `T_w = 60 log(1/lambda)` in (A.9) is
   reverse-engineered to give thirty powers of `lambda` of room.**

### 3.6 Summary table of the steps I judge under-justified

| # | Step | Page | Severity | Comment |
|---|---|---|---|---|
| 1 | Coverage: `I = [X_-, X_+]` contains every possible cone failure | 158 | **high** | Never stated as a single claim; assembled from A.4, A.7, A.10, B.6; the seam between the axis collar and A.4's admissible region is not displayed |
| 2 | Finiteness/closure of the list of `eta`-derivative orders needing smallness | 127, 133, 140, 157, 162 | **high** | Each use of (4.17) costs one order; the chain A.7 -> A.8 -> C.2 -> C.3 is never counted |
| 3 | `N`-dependence of radial derivatives of the modulated profile (`C_{r,m}N^{r-1}`) | 161, 163 | **medium-high** | Honestly flagged, but every downstream use of a radial profile derivative inherits `N^{r-1}`; verification lives outside Appendices A-C |
| 4 | `||B^{-1}|| = O(lambda^{-1})` for the coalescing `1, x^{-lambda}` block | 127-128, 162 | medium | One-sentence proof for a hinge used by three correction stages |
| 5 | The `o(1)` / `h << lambda` closure across Section A.5 | 134-137 | medium | Cone margins rest on `lambda`, `h/lambda` small; I could not verify that no later choice feeds back |
| 6 | Vanishing of the cubic term making `sqrt(V)` smooth at `mu = 0` | 160 | low-medium | True (`<sin^3> = 0`) but never mentioned, and without it (C.10) is not smooth through `rho = 0` |
| 7 | The four backward moment representations in Lemma A.8 | 140 | low | Asserted in one sentence; I verified all four |
| 8 | `-Z Hcal'/Hcal < h/4` for large `X_R` in (A.56) | 143 | low | Immediate from (A.36) but not shown |
| 9 | Ordering `X_a < X_an < X_-`, and the two senses of "collar" | 153, 158 | low | Ambiguous phrasing, consistent on my reading |
| 10 | Smoothness of `V_0/X` at the axis for the Proposition B.2 profile | 148 | low | Follows from (4.7) since `U`, `A_X(U)` are analytic; not stated |
| 11 | The `lambda`-exponent bookkeeping `29, 30, 26, 27, 28` in (A.14)-(A.15), (A.30) | 131, 136 | low | I verified `30`; the rest are plausible but unchecked |
| 12 | `Pi_0`'s sign properties survive reinstating the (A.11) bumps | 133-134 | low | Only the *total* increment is argued unchanged; that is all that is used, but the reader must notice |
