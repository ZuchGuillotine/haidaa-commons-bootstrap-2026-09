# Reader's report: Sections 6-8 (pp. 62-100)

**Paper.** "Finite time blowup for Navier-Stokes", author line OpenAI, compiled 2026-09-08,
166 pp. Source: `/Users/benjamincox/dsm:haidaa/navier-stokes.pdf`.
**Range assigned.** Section 6 "Auxiliary torus and separation of oscillatory supports"
(pp. 62-73), Section 7 "Oscillatory realization and correction of the residual stress"
(pp. 73-88), Section 8 "Compactly supported mean corrections" (pp. 88-100).
Read with the Read tool in two chunks, pp. 62-81 and pp. 82-100.

All math below is transcribed in plain text. Equation numbers are the paper's.
Everything I could actually check by hand is marked **[verified]**; everything I could
not is marked **[unverified]** or **[external]**.

---

## 0. Executive summary

Sections 6-8 are the "hard analysis engine room" of the construction. They contain no
statement of the form "Navier-Stokes blows up"; they build three toolkits:

1. **Section 6** — a bookkeeping/geometry layer: dyadic charts, an auxiliary 2-torus
   `Y` with a Diophantine integer matrix `J_g`, a device (Lemma 6.1) that gives distinct
   oscillations *disjoint* auxiliary supports even when their physical supports overlap,
   a "common torus" for overlapping dyadic bands (Lemma 6.2), a counting lemma
   (Lemma 6.3), and three coefficient classes `M_a` (means/moments), `W_a` (wave
   amplitudes), `S_a` (radial moments) with a product/derivative algebra
   (Definitions 6.4-6.5, Proposition 6.6).
2. **Section 7** — the wave construction: a phase `Phi`, a 2x2 amplitude ODE along a
   "pulse coordinate" `v`, a growth-then-decay envelope `P(v)` (7.12) with Gaussian
   bounds (7.16), a zero-data inverse for prescribed harmonic sources
   (Proposition 7.2), a homogeneous pulse (Lemma 7.4), a positive-coefficient
   covariance match to the target stress (Proposition 7.5), its linearization
   (Proposition 7.6), exact incompressibility by taking curls (Lemma 7.7), and a
   quadratic-remainder identity (Corollary 7.8).
3. **Section 8** — the angular-mean corrections: the mean momentum system (8.3), a
   compactly supported radial primitive (Lemma 8.2), pressure and divergence-free
   axial increments (Proposition 8.3), three integral compatibility defects
   `(P, J_theta, J_z)` (8.12), (8.15), (8.16), covariance targets `H_theta, H_z`
   (Corollary 8.5), inversion of the fast time derivative `N^{-1}` (Lemma 8.6), and a
   5x5 linear moment correction (Lemma 8.7, system (8.25)) plus its exact recomputation
   (Lemma 8.8).

**My overall read.** The material is internally coherent and, where I could check it,
*algebraically correct to a degree I did not expect* — I verified roughly a dozen
non-trivial computations exactly (list in §6). The two devices the prompt asked me to
scrutinize hardest — the auxiliary-torus separation and the amplification-then-damping
envelope — both survive scrutiny as *designs*: the separation is arranged on the
absolute torus (all deck copies), so physical evaluation `Y(r,t)` cannot reintroduce
cross terms; and the growth/damping crossover is engineered exactly (the reference
damping `d_ref` is *defined* so that the net rate vanishes at `|s| = u_*`).

However: **nothing in Sections 6-8 demonstrates a net gain**, and therefore nothing in
my range can establish convergence. Every inverse operator in my range is
gain-neutral or nearly so (Proposition 7.2: "no decrease of `alpha`"; Lemma 8.6:
"retain their `eps` exponent"; Lemma 8.2 gains only on the *cutoff remainder*). The
entire claim that "each cycle gains a fixed positive power of `eps`" is deferred to
Proposition 9.6 (stated on p. 100). Likewise all uniformity across *stages* and across
*derivative orders* — the two places where a scheme of this kind normally fails — is
asserted here with `C_{j,I}` allowed to depend on stage `j` and multi-index `I`, and
the diagonal/summation argument is deferred to Lemma 9.7 and Proposition 9.9.
So **Sections 6-8 cannot be accepted or rejected on their own**; they can only be
checked for internal consistency, which is what I did.

I found **one apparent error** (a missing weight in (7.25)) and **eight places** where
uniformity is asserted rather than proved. Details in §3.

---

## 1. Precise statements (plain-text math, with page numbers)

### 1.1 Section 6 set-up

**Dyadic charts (6.1), p. 62.** For a dyadic index `l`,
```
Q = 2^{-l},   eps = Q^h,   S_* = l^2,   R = r/sqrt(Q),   Z = z/Q^D,   T = tau/Q,   tau = 1-t
```
with `A = 1/2 + h`, `D = 1/2 - h` from (4.1). Physical velocity, pressure, residual have
chart representatives obtained by multiplication by `Q^A`, `Q^{2A}`, `Q^{2A+1/2}` (p. 63).

**The torus matrix (6.2), p. 63.**
```
J_g = [[3,1],[1,5]],  Lam_g = 4 - sqrt(2),  T_g = 4 + sqrt(2),  b_g = sqrt(2) - 1,
v_r = (1, -b_g),  v_t = (b_g, 1),  rho_g = log(Lam_g)/log(T_g),
kap_s = 10^{-5},  d_r = 2((1+h) rho_g - h kap_s) > 0.
```
`J_g v_r = Lam_g v_r`, `J_g v_t = T_g v_t`, `1 < Lam_g < T_g`. **[verified]**: trace 8,
det 14, eigenvalues `4 +- sqrt 2`; eigenvectors as displayed; `rho_g ~ 0.5626`.

**Physical evaluation map (6.3), p. 63.**
```
Y = v_r r^{d_r} + v_t t   (mod Z^2).
```
**Evaluated derivative operators (6.4), (6.6), pp. 63.**
```
N_abs = v_t . d_Y,   L_abs = v_r . d_Y,
tt = d_t + N_abs,    rr = d_r + d_r r^{d_r - 1} L_abs
t_* = Q^{1+h} tt = -eps d_T + c_i N_i,       c_i = T_g^i Q^{1+h} ~ S_*^{-1}
D_r = sqrt(Q) rr = d_R + M_i d_r R^{d_r-1} L_i,   M_i = Lam_g^i Q^{d_r/2} ~ eps^{-kap_s} S_*^{-rho_g}
D_z = sqrt(Q) d_z = eps d_Z,                  D_theta = R^{-1} d_theta
```
with band representative `Y_i = J_g^i Y (mod Z^2)`, `i = i(l) = floor(log_{T_g}(Q^{-1-h}/S_*))`
(6.5). **[verified]** the two asymptotic identities: `c_i = T_g^i Q^{1+h}` and
`T_g^i ~ Q^{-1-h}/S_*` give `c_i ~ S_*^{-1}`; and
`M_i = Lam_g^i Q^{d_r/2} = (T_g^i)^{rho_g} Q^{d_r/2} = Q^{-h kap_s} S_*^{-rho_g} = eps^{-kap_s} S_*^{-rho_g}`
using exactly the definition of `d_r`. This is a tight, non-accidental choice: `d_r` is
*defined* to make the radial-derivative loss equal to a single factor `eps^{-kap_s}`.

**Diophantine lower bounds (6.7), p. 64.** For `n` in `Z^2 \ {0}`,
```
|v_r . n| >= c/(1+|n|),    |v_t . n| >= c/(1+|n|).
```
**[verified]**: `v_r . n = (n_1+n_2) - sqrt(2) n_2`; its algebraic conjugate
`(n_1+n_2) + sqrt(2) n_2` has modulus `<= C|n|`, and the product is the nonzero integer
`(n_1+n_2)^2 - 2 n_2^2`. Same for `v_t . n = (n_2-n_1) + sqrt(2) n_1`. These bounds are
used exactly twice: (8.10) in Lemma 8.2 and (8.19) in Lemma 8.6.

**Labels and slow cutoffs (6.8), (6.9), pp. 64-65.** Squared partitions
`sum_l chi_l(q)^2 = 1` (supported where `q/Q` is in `[1/2,2]`) and, in each band, a
product squared partition `sum_a chi_{l,a}(R,Z,T)^2 = 1` with mesh `S_*^{-3}` in each
coordinate. Labels `Gam = {(l,a,sigma) : l >= l_0, a in I_l, sigma in {+,-}}`;
`eta_gam = chi_l(q) chi_{l,a}(R,Z,T)`, `K_gam = supp(eta_gam o C_l)`;
`eta_{(l,a,+)} = eta_{(l,a,-)}` and `sum_l sum_a (eta_{(l,a,+)} o C_l)^2 = 1`.

**Auxiliary rectangles (6.10), (6.11), p. 65.**
```
R_gam   = {c_gam + xi v_r + eta v_t : |xi|,|eta| < r_0}       (mod Z^2)
R_gam^+ = same with 2 r_0
R_gam^abs = pi_{i(l)}^{-1}(R_gam),   R_gam^{abs,+} = pi_{i(l)}^{-1}(R_gam^+)
Y_i - c_gam - k_copy = xi_g v_r + eta_g v_t,   v = (eta_g + r_0)/c_i,   L_s = 2 r_0/c_i ~ S_*
```
and (6.12) `D_r v = D_z v = 0`, `t_* v = 1`, `N_i xi_g = 0`. So `v` is normalized time
along the rectangle, advancing at unit speed under the *fast* time derivative, while
being invariant under the normalized spatial derivatives.

> **LEMMA 6.1 (separation of supports), p. 65.** *For the labels `Gam` and slow supports
> `K_gam` in (6.8) and (6.9), there exist centers `c_gam` and a radius `r_0 > 0`, common
> to all labels and independent of the band, such that the enlarged rectangles in (6.10)
> are injectively parametrized and*
> ```
> gam != gam',  K_gam ∩ K_gam' != empty   ==>   closure(R_gam^{abs,+}) ∩ closure(R_gam'^{abs,+}) = empty.   (6.13)
> ```

Proof in three steps (pp. 65-66): (1) the "interaction graph" (join two labels when
enlarged slow boxes meet in physical coordinates and band indices differ by at most 4;
also join the two signs of one box) has bounded degree, uniformly in the band, and
```
|i(l) - i(l')| <= 1 + (4(1+h) log 2 + 8/l_0)/log T_g =: Delta_max   for |l - l'| <= 4.   (6.14)
```
**[verified]** exactly: `i(l) = floor(((1+h) l log 2 - 2 log l)/log T_g)`, the mean value
theorem gives `2|log l - log l'| <= 8/l_0`, and flooring adds at most one.
(2) greedily color the graph with finitely many colors, one rational center `c_nu` per
color, subject to the finitely many ordered constraints
```
c_mu != J_g^Delta c_nu (mod Z^2)   for 0 <= Delta <= Delta_max,  unless (Delta,nu,mu) = (0,nu,nu).   (6.15)
```
Self-pairs with `Delta > 0` use invertibility of `J_g^Delta - I`. (3) choose one `r_0`,
small enough that a hypothetical common point would force
`c_mu - J_g^Delta c_nu = J_g^Delta e_nu - e_mu` mod `Z^2` with `|e_nu| + |e_mu| <= C r_0`,
contradicting the positive separation from (6.15).

**Consequences, p. 66.** For a family `(F_gam)` with `supp F_gam ⊂ K_gam x R_gam^{abs,+}`,
`F_gam F_gam' = 0` for `gam != gam'`, "also for derivatives of smoothly extended fields
and after evaluation on (6.3)". Transverse cutoff `chi_g in C_c^inf((-r_0,r_0))`; time
cutoff
```
psi(v) = 1 if |v - L_s/2| <= L_s/5,   supp psi ⊂ {|v - L_s/2| < L_s/3}.   (6.16)
```

> **LEMMA 6.2 (common torus for overlapping bands), p. 67.** *For a neighborhood `U`,
> `i_0 = min_{l in L(U)} i(l)`, `H = pi_{i_0}(Y)` as in (6.17), each change `H -> Y_i`
> has uniformly bounded derivative factors, and `R_gam^H` consists of
> `14^{Delta_l} <= 14^{Delta_max}` disjoint lifts of `R_gam`. Haar averages on the
> absolute, common and band tori agree whenever the function is a pullback from the
> respective torus. Radial integration at fixed `(z,t)`, torus translation, averaging,
> and directional Fourier inversion on zero-mean functions preserve common-torus descent
> and introduce no new dyadic band.*

**[verified]** the fiber count: `#(Z^2/J_g^Delta Z^2) = |det J_g|^Delta = 14^Delta`.
Compatibility (6.19) `int f(J_g^Delta H) dH = int f(Y_i) dY_i` is checked on characters
(a character of frequency `n` pulls back to `(J_g^Delta)^T n`, zero iff `n = 0`).
Chart compatibility on overlaps is (6.20).

> **LEMMA 6.3 (number of relevant labels), p. 68.** *There is `C` independent of band and
> point with `sup_{(r,z,t)} #{gam in Gam : (r,z,t) in K_gam} <= C`. For a fixed band, a
> compact chart set `B ⊂ R^3` and a bounded normalized radial interval `I_R`,*
> ```
> #{(l,a,sigma) : B_{l,a} ∩ B != empty} <= C_B S_l^9,
> sup_{Z,T} #{(l,a,sigma) : B_{l,a} ∩ (I_R x {Z} x {T}) != empty} <= C_{I_R} S_l^3.
> ```

So: **pointwise** only `O(1)` waves overlap, while **integrated** counts grow like
`S_*^9 = l^18` — polynomially, which is the whole point (polynomial growth in `S_*` is
free against any gained power of `eps = 2^{-h l}`).

**Path along a rectangle (6.21), (6.22), p. 68.** With `lam_t(w) = v_t . w/(1+b_g^2)`
(so `lam_t(v_r) = 0`, `lam_t(v_t) = 1`) and `eta_g(H) = lam_t(J_g^Delta H - c_gam - k_copy)`,
```
H_w = H + T_g^{-Delta} c_i (w - v(H)) v_t,     D_H H_w = I - v_t lam_t,   D_H v = T_g^Delta lam_t / c_i = O(S_*).
```
Higher derivatives vanish; a deck translation translates the path by the same amount.

> **DEFINITION 6.4 (mean and moment coefficients), p. 69.** *A smooth `f = f(R,Z,T,H)`
> belongs to `M_alpha` if it descends to each applicable common torus with representatives
> satisfying (6.20), extends smoothly by zero outside the active radial shell, and on
> `X_a < X < X_b`*
> ```
> d_theta f = 0,   for all I in N_0^5 there exist C_{j,I} > 0, b_{j,I}, d_{j,I} >= 0:
> |d^I f| <= C_{j,I} eps^alpha S_*^{b_{j,I}} zeta delta^{-d_{j,I}}.   (6.24)
> ```
> *For `F = F(Z,T)`:*
> ```
> F in S_alpha  <==>  |d_{Z,T}^I F| <= C_{j,I} eps^alpha S_*^{b_{j,I}}.   (6.25)
> ```

Here (6.23) is the flat edge weight
```
zeta(X) = exp(-a_a/log^2(X/X_a) - a_b/log^2(X_b/X)),   delta(X) = min{1, log(X/X_a), log(X_b/X)}
```
extended by zero, with the key property that `zeta^b delta^{-N} -> 0` faster than any
power at either boundary, for every `b > 0` and finite `N`. Radial moment normalization
(8.26 in the paper's numbering is (6.26)):
```
M_*(Z,T) = int R^e <f_*>_Y dR = Q^{a_f - (e+1)/2} int r^e <f_phys>_Y dr.   (6.26)
```
Wave data (6.27): `Phi_gam = p_gam theta + phi_gam(R,Z,T,H)`, `k = ceil(eps^{-1/2})`,
`k p_gam in Z \ {0}`, `0 < P_v <= 1` on `0 <= v <= L_s`, finite harmonic set `H_j`
independent of the band. Domains (6.28): `E_gam`, `Om_gam`, `Om_gam^cut`.

> **DEFINITION 6.5 (labelled wave coefficients), p. 70.** *For `gam in Gam`, `m in Z\{0}`, a
> smooth amplitude `a_{gam,m}` on `E_gam^+` belongs to `W_alpha` if its local formulas agree
> on lifts representing the same common-torus point, satisfy (6.20), and*
> ```
> d_theta a_{gam,m} = 0,   supp(a_{gam,m}|E_gam) ⊂ Om_gam,
> |d^I a_{gam,m}| <= C_{j,I} eps^alpha S_*^{b_{j,I}} sqrt(zeta) delta^{-d_{j,I}} P_v,  0 <= v <= L_s.   (6.29)
> ```
> *The continuation beyond `v = 0, L_s` has no envelope bound and no temporal
> zero-extension condition. A cut-off coefficient additionally has
> `supp a^cut ⊂ Om_gam^cut` and smooth zero extension in `v`. No divisibility by
> `eta_gam, chi_g, psi` is assumed.*

Note the asymmetry, which is used throughout: means carry `zeta`, waves carry
`sqrt(zeta)` — so a product of two waves lands back in the mean class with the mean's
weight. A real wave `w = sum_gam sum_{m in H_j} A_{gam,m} e^{i k_gam m Phi_gam}` is in
`W_alpha` iff every `A_{gam,m}` is (p. 71).

> **PROPOSITION 6.6 (coefficient algebra), p. 71.** *Each class is closed under finite sums
> at a fixed exponent, and*
> ```
> M_a M_b ⊂ M_{a+b},   M_a W_b ⊂ W_{a+b}.                                        (6.30)
> a_{gam,m} b_{gam,m'} in  W_{a+b} if m+m' != 0 (harmonic m+m'),
>                          M_{a+b} if m+m' = 0, after temporal cutoff.            (6.31)
> ```
> *Actual waves with distinct labels satisfy `(d^I w_gam)(d^J w_gam') = 0` for `gam != gam'`.
> For `C in {M, W}`:*
> ```
> d_{R,Z,T} : C_a -> C_a
> D_r       : C_a -> C_{a - kap_s}
> D_z = eps d_Z : C_a -> C_{a+1}
> Q^{1+h} N_abs = c_{i_0} N_{i_0} : C_a -> C_a
> -eps d_T   : C_a -> C_{a+1}                                                    (6.32)
> ```

The proof (pp. 71-72) uses `zeta^2 <= zeta`, `zeta^{3/2} P_v <= sqrt(zeta) P_v`,
`zeta P_v^2 <= zeta` (`m+m'=0`) or `<= sqrt(zeta) P_v` (`m+m' != 0`), Leibniz, and the
angular-average identity
```
<a_{gam,m} a_{gam,m'} e^{i k (m+m') Phi_gam}>_theta = a_{gam,m} a_{gam,m'} if m+m'=0, else 0
```
which holds because `k p_gam` is a nonzero integer and the amplitudes are
`theta`-independent. **[verified]** given (7.3) (`Phi` is `p theta` plus a
`theta`-independent remainder).

### 1.2 Section 7 — the waves

**Frame and growth parameters, p. 74.** Chart base velocity `b e_theta + V e_theta + G e_z`,
`F = V/R`, `g = (R F_R, G_R)` the tangential shear in order `(theta,z)`; subscript 0 means
order-zero, frozen at the representative. Then
```
N = g_0/|g_0|,  K = N^perp = (-N_z, N_theta),
lam_0^2 = -2 F_0 N_theta (2 F_0 N_theta + |g_0|) > 0,   c_0 = lam_0/(2 F_0 N_theta) < 0
```
with `g_0 = F_0(-a, b_s)`, `lam_0^2 = 2 a F_0^2 (1 - 2/v_s)` from (4.11), (4.20).

**Strict cone condition (7.1), p. 74.**
```
T_{0,*} . N < 0,     |c_0 (T_{0,*} . K)/(T_{0,*} . N)| < 1,     T_{0,*} = (Q/q)^{A+1/2} T_0.
```
On p. 74 one then chooses `u_* > 0` with `u_*/sqrt(1+u_*^2)` exceeding the supremum of
`|c_0 (T_{0,*}.K)/(T_{0,*}.N)|` over the closed annulus, and freezes `N, K, c_0` at each
representative.

**Carrier and wavenumbers (7.2), (7.3), (7.4), p. 74.**
```
k = ceil(eps^{-1/2}),   B_s^2 = lam_0/(eps k^2 (1+u_*^2)^{3/2}),   s(v) = sigma(u_*/2 + u_* v/L_s)
(p~/R_0, p_z) = B_s (K - sigma u_* g_0/(L_s |g_0|^2)),  then k p := nearest nonzero integer to k p~
Phi = p theta + p_z Z/eps + x_0 R - v(p F + p_z G),   x_0 = sigma B_s u_*/2
n_Phi := grad_* Phi = (x_0 - v(p F_R + p_z G_R),  p/R,  p_z - eps v(p F_Z + p_z G_Z)).
```

**Principal amplitude problem (7.5), (7.6), p. 75.** For a source `f_m in C^3` and
`m != 0`, find `t_m in C^3` with `n_Phi . t_m = 0` and pressure `pi_m in C`:
```
t_m' + K t_m + m^2 d t_m + i k m n_Phi pi_m = -f_m,     d = eps k^2 |n_Phi|^2,
K = [[0, -2F, 0], [2F + R F_R, 0, 0], [G_R, 0, 0]],
A_Phi = -K + n_Phi (n_Phi^T K - (n_Phi')^T)/|n_Phi|^2.
```
Geometric quantities (7.7), (7.8): `K_a = n_{Phi,tan}/|n_{Phi,tan}|`,
`N_a = ((K_a)_z, -(K_a)_theta)`, `s_a = n_{Phi,r}/|n_{Phi,tan}|`,
`U = [e_r - s_a K_a, N_a]`, `B = U [[1, 1],[c_0 sqrt(1+s^2), -c_0 sqrt(1+s^2)]]`.

> **LEMMA 7.1 (phase and frame), p. 75.** *For the labels and phases in (6.8), (7.2), (7.3)
> there is `q_* > 0`, common to all labels and all fixed derivative orders, such that every
> rectangle on `0 < q < q_*` has: `e^{i k m Phi}` single-valued in `theta` with nonzero
> angular frequency, and on `0 <= v <= L_s`*
> ```
> |n_Phi - B_s(s(v), K)| + |n_Phi'| <= C/S_*,   E_ik := (t_* + b D_r + F d_theta + G D_z) Phi = O(eps S_*^C).   (7.9)
> ```
> *`B` in (7.8) has a uniformly bounded left inverse `B^l` and*
> ```
> B^l (A_Phi B - B') = diag(lam, -lam) + E,   lam(v) = lam_0/sqrt(1+s(v)^2),   |E| <= C/S_*.   (7.10)
> ```
> *With `d_ref = eps k^2 B_s^2 (1+s^2)`,  `|d - d_ref| <= C/S_*`.   (7.11)*

**Envelope (7.12), p. 77.**
```
P(v) = exp( int_{L_s/2}^{v} (lam(w) - d_ref(w)) dw ),    P_v = P(v).
```

> **PROPOSITION 7.2 (zero-data inverse), pp. 77-78.** *Fix `alpha`, a finite set of harmonics
> `0 < |m| <= M`, and a label at level `i` on a common torus `H = Y_{i_0}`. For each `m`, let
> `f_m(x_*,H) in C^3` be a coefficient descending to `H`, smooth on the fixed enlargement of
> the label's rectangles, with fixed closed slow support in the label's slow support, compact
> transverse support in `supp chi_g`, support in the closed shell `X_a <= X <= X_b`, smooth
> zero extension across all these, and `|D^I f_m| <= C_I eps^alpha S_*^{b_I} sqrt(zeta) delta^{-a_I} P(v)`.
> Then on `0 < q < q_*` there is a unique `t_m` with `n_Phi . t_m = 0`, `t_m(0) = 0`, and a
> pressure `pi_m` such that (7.5) holds, given by*
> ```
> t_m' = A_Phi t_m - m^2 d t_m - proj_{n_Phi^perp} f_m,   t_m(0) = 0,
> pi_m = (i/(km)) (n_Phi . K t_m - n_Phi' . t_m + n_Phi . f_m)/|n_Phi|^2.   (7.13)
> ```
> *The source is evaluated along the path `H_w` of (6.21). On `X_a < X < X_b`, for every fixed `I`,*
> ```
> |D^I t_m|  <= C_I' eps^alpha       S_*^{b_I'} sqrt(zeta) delta^{-a_I'} P(v),     (7.14)
> |D^I pi_m| <= C_I' eps^{alpha+1/2} S_*^{b_I'} sqrt(zeta) delta^{-a_I'} P(v).     (7.15)
> ```
> *Constants and degrees may depend on `M, I` and the source bounds. **The threshold `q_*`
> does not depend on `M`, `alpha`, or the correction stage.** One-sided derivatives at
> `tau = 0` are inherited from the source.*

Key steps: (7.16) the two-sided Gaussian bound
```
e^{-C(v - L_s/2)^2/L_s} <= P(v) <= e^{-c(v - L_s/2)^2/L_s} <= 1,
```
derived from `(log P)' = a_net(|s(v)|)` with
```
a_net(y) = lam_0/sqrt(1+y^2) - lam_0 (1+y^2)/(1+u_*^2)^{3/2},
d/dv a_net(|s(v)|) = -(u_* lam_0 |s(v)|/L_s)((1+s^2)^{-3/2} + 2(1+u_*^2)^{-3/2}) in [-C/L_s, -c/L_s];
```
then `t_m = B z_m`, `z_m' = A_m z_m + g_m`, `A_m = diag(lam,-lam) + E - m^2 d I_2` (7.17), the
exact harmonic identity
```
V_m(v,w) = exp(-(m^2-1) int_w^v d(a) da) V_1(v,w)     (7.18)
```
and the propagator bound
```
||V_m(v,w)|| <= C P(v)/P(w),   0 <= w <= v <= L_s,  m != 0.     (7.19)
```
Derivatives: (7.20) differentiated Duhamel, with each derivative raising the polynomial
degree by at most `1 + max{b_I, max_{0<J<=I}(c_J + b_{I-J}')}`.

> **COROLLARY 7.3, p. 80.** *The solution operators `f_m -> (t_m, pi_m)` are defined at every
> finite correction stage on the same domain `0 < q < q_*`. Increasing the finite harmonic
> set or the derivative order changes constants and polynomial degrees but does not shrink
> the domain.* (Proof: (7.18) makes higher harmonics at least as damped.)

> **LEMMA 7.4 (the homogeneous pulse), p. 80.** *For each sign `sigma`, let
> `t_sigma^h = B(z_+, z_-)^T` solve `(t^h)' = A_Phi t^h - d t^h` with `z_+(0) = P(0)`,
> `z_-(0) = 0`. Writing `t_sigma^h = x(e_r - s_a K_a) + y N_a`,*
> ```
> c P(v) <= x(v) <= C P(v),    y(v)/x(v) = c_0 sqrt(1 + s(v)^2) + O(S_*^{-1}).   (7.21)
> ```
> *Every fixed derivative of this homogeneous solution in the coefficient variables or the
> pulse coordinate is bounded by `C_I S_*^{b_I} P(v)`. Its homogeneous pressure has the extra
> factor `eps^{1/2}`.*

**Energy transfer from the shear (7.22), p. 81.** For the real amplitude `t = t_sigma^h`,
```
(1/2) d|t|^2/dv = -g . (t_r(t_theta, t_z)) - eps k^2 |n_Phi|^2 |t|^2,
```
"the first term gives energy exchange with the base shear; the second gives viscous
dissipation." The flux vector is decomposed as `T = T_N N + T_K K`; the energy transfer is
`-g_0 . T = -|g_0| T_N`, and the covariance construction also prescribes `T_K` — hence
**both** components are needed.

**Scale bookkeeping, p. 81.** "At viscosity one, the base velocity and radial length give
Reynolds number of order `q^{-h}`. The carrier wavelength is `sqrt(Q)/k`, whereas the
physical amplitude of the leading oscillations has order `Q^{-A} sqrt(eps) = Q^{-1/2-h/2}`,
up to powers of `S_*`. The powers of `Q` in their product cancel. Thus `eps k^2 ~ 1` in the
pulse equation throughout the correction process."

**Covariance functionals (7.23), p. 82.**
```
C(w) = <<w_r w_tan>_theta>_Y,     B(w,v) = <<w_r v_tan + v_r w_tan>_theta>_Y,   D C(w)[v] = B(w,v).
```
For `beta = (l,a)` and `gam = (beta,sigma)`, `b_sigma = chi_g(xi_g) psi(v) t_sigma^h cos(k Phi_sigma)`;
`H = [C(b_+) | C(b_-)]` a real 2x2 matrix with rows in order `(r theta, r z)`, and
disjointness gives `C(a_+ b_+ + a_- b_-) = H (a_+^2, a_-^2)^T`.

> **PROPOSITION 7.5 (covariance matching), p. 82.** *For each slow box `beta`, use the
> homogeneous pulses of Lemma 7.4 to form `b_+-` and `H`, and the order-zero target
> `T_{0,*} = (Q/q)^{A+1/2} T_0`. After a further single decrease of `q_*`, `H` is invertible
> throughout each enlarged slow neighborhood, and the components of `H^{-1} T_{0,*}` are
> positive on its intersection with `X_a < X < X_b`. Define*
> ```
> y = H^{-1} T_{0,*},    a_sigma = sqrt(y_sigma),    W_0 = sqrt(eps) sum_{sigma=+-} a_sigma b_sigma.   (7.24)
> ```
> *The squared amplitudes satisfy, for every fixed slow multi-index `I` on `X_a < X < X_b`,*
> ```
> y_sigma >= c sqrt(S_*) zeta,      |D^I y_sigma| <= C_I S_*^{b_I} delta^{-a_I}.   (7.25)
> ```
> *The `a_sigma` extend smoothly by zero off each sign's auxiliary rectangle. The harmonic
> coefficients of `W_0` obey the estimates of `W_{1/2}`; `eta_beta W_0 in W_{1/2}`. Before
> that multiplication,*
> ```
> C(W_0) = eps T_{0,*}.   (7.26)
> ```

Proof ingredients: the column formula
```
H_sigma = (|det(v_r,v_t)|/2) (int chi_g^2 dxi_g) c_i int_0^{L_s} psi^2 x t_{sigma,tan}^h dv   (7.27)
h_sigma = (|det(v_r,v_t)|/2) (int chi_g^2 dxi_g) c_i int_0^{L_s} psi^2 x^2 dv,
c sqrt(L_s) <= int psi^2 x^2 dv <= C sqrt(L_s),   int (|v-L_s/2|/L_s) psi^2 x^2 dv <= C
h_sigma ~ c_i sqrt(L_s) ~ S_*^{-1/2}
H_sigma = h_sigma(-A_c N - sigma u_* K + e_sigma),   |e_sigma| <= C S_*^{-1/2},  A_c = -c_0 sqrt(1+u_*^2) > 0   (7.28)
h_+ y_+ = (1/2)(-T_N/A_c - T_K/u_*),   h_- y_- = (1/2)(-T_N/A_c + T_K/u_*)
|det H| >= c/S_*,   ||H^{-1}|| <= C sqrt(S_*),   c sqrt(S_*)|T| <= y_sigma <= C sqrt(S_*)|T|.   (7.29)
```

**Global assembly (7.30), p. 84.**
```
sum_{beta=(l,a)} Q^{-2A} eta_beta^2 eps T_{0,*} = q^{-A-1/2} T_0.
```

**Linearization (7.31), p. 84.** For a real two-vector `Sigma(R,Z,T) in M_alpha` *independent
of the angular and auxiliary variables*,
```
d_Sigma = H^{-1}(Sigma/eps),    delta a_sigma = (d_Sigma)_sigma/(2 a_sigma),    L Sigma = sqrt(eps) sum_sigma delta a_sigma b_sigma.
```

> **PROPOSITION 7.6, p. 84.** *For every `alpha` and every real two-vector
> `Sigma(R,Z,T) in M_alpha` independent of the angular and auxiliary variables, (7.31) defines
> a linear operation `L Sigma` whose harmonic coefficients obey the local derivative and
> envelope bounds of `W_{alpha-1/2}` with smooth zero extension at the shell edges, and*
> ```
> eta_beta L Sigma in W_{alpha-1/2},    B(W_0, L Sigma) = Sigma,    C(eta_beta L Sigma) in M_{2 alpha - 1}.   (7.32)
> ```
> *The denominator `a_sigma` remains this fixed coefficient at every correction stage. No sign
> or relative-size restriction is imposed on `Sigma`.*

Proof: `|D^I (d_Sigma)_sigma| <= C_I eps^{alpha-1} S_*^{b_I} zeta delta^{-a_I}`, hence
```
|D^I delta a_sigma| <= C_I eps^{alpha-1} S_*^{b_I} sqrt(zeta) delta^{-a_I},   (7.33)
B(W_0, L Sigma) = eps H (2 a . delta a) = eps H d_Sigma = Sigma.   (7.34)
```
Consequently (p. 85) `C(W_0 + L Sigma) = eps T_{0,*} + Sigma + C(L Sigma)`, i.e. `L` is a right
inverse of `D C(W_0)` up to a quadratic remainder. Global versions:
```
W_{0,phys}^as = sum_beta Q^{-A} eta_beta W_{0,beta},   L_phys^as sigma = sum_beta Q^{-A} eta_beta L_beta(Q^{2A} sigma),  (7.35)
W_0^as = Q^A W_{0,phys}^as in W_{1/2},   L^as Sigma = Q^A L_phys^as(Q^{-2A} Sigma) in W_{alpha-1/2},
B(W_0^as, L^as Sigma) = Sigma.   (7.36)
```

**Exact curls (7.37), (7.38), pp. 85-86.**
```
curl_* A = (R^{-1} d_theta A_z - D_z A_theta,  D_z A_r - D_r A_z,  (D_r + R^{-1})A_theta - R^{-1} d_theta A_r)
C_m = i n_Phi x t_m/(k m |n_Phi|^2),   A_m = C_m e^{i k m Phi},   curl_* A_m = (t_m + r_m) e^{i k m Phi}
r_m = (-D_z(C_m)_theta, D_z(C_m)_r - D_r(C_m)_z, (D_r + R^{-1})(C_m)_theta).
```

> **LEMMA 7.7 (exact incompressibility), p. 86.** *For every `alpha` and nonzero integer `m`,
> let `t_m in W_alpha` satisfy `n_Phi . t_m = 0`, supported in the prescribed slow and
> transverse supports, compactly inside the pulse interval, with smooth zero extension.
> Then `C_m in W_{alpha+1/2}` and `r_m in W_{alpha+1/2-kap_s}`. The velocity
> `(t_m + r_m) e^{i k m Phi} = curl_* A_m` is exactly divergence-free for the normalized
> operators and after evaluation at the phase map. For `a_m = t_m + r_m`,*
> ```
> i k m n_Phi . a_m = -((D_r + R^{-1})(a_m)_r + D_z(a_m)_z),    n_Phi . a_m in W_{alpha+1/2-kap_s}.   (7.39)
> ```
> *Physical potentials and pressures come from chart ones by `Q^{1/2-A}` and `Q^{-2A}`.*

**[verified]** `i k m n_Phi x C_m = t_m - n_Phi(n_Phi . t_m)/|n_Phi|^2 = t_m`, and the
six-term `div_* curl_* = 0` cancellation using `D_r(R^{-1}) = -R^{-2}` and
`(D_r + R^{-1})(R^{-1} f) = R^{-1} D_r f`.

**Temporal cutoff residual (7.40), p. 87.** With `t^_m = psi t_m`, `pi^_m = psi pi_m`:
```
(t^_m)' + K t^_m + m^2 d t^_m + i k m n_Phi pi^_m + f_m = (1 - psi) f_m + psi' t_m.
```
Wherever either right-hand term survives, `|v - L_s/2| >= L_s/5`, so (7.16) supplies an extra
factor `S_*^C e^{-c S_*}` in every fixed derivative estimate; and
```
q^{-N} Q^{-M} S_*^C e^{-c S_*} <= C_N exp(-c l^2 + (M+N) l log 2 + 2 C log l) -> 0.
```

> **COROLLARY 7.8 (covariance of the actual velocities), p. 87.** *Suppose `U` is a real
> divergence-free wave field, including curl remainders, with `U = W_0^as + E`, `U in W_{1/2}`,
> `E in W_rho`. Let `Sigma(R,Z,T) in M_alpha` be a real two-vector independent of the angular
> and auxiliary variables. Apply Lemma 7.7 to `L^as Sigma`, obtaining the divergence-free
> increment `V = L^as Sigma + R`, `R in W_{alpha-kap_s}`. Then*
> ```
> C(U+V) - C(U) = Sigma + B(E, L^as Sigma) + B(U, R) + C(V),          (7.41)
> B(E, L^as Sigma) in M_{rho + alpha - 1/2},   B(U,R) in M_{alpha + 1/2 - kap_s},   C(V) in M_{2 alpha - 1}.   (7.42)
> ```
> *Taking the radial divergence of each remainder decreases the `eps` exponent by at most a
> further `kap_s`.*

### 1.3 Section 8 — mean corrections

**Decomposition (8.1), (8.2), p. 88.**
```
u = (b + beta, V + v, G + gamma) + w,   <w>_theta = 0,   W_ab = <w_a w_b>_theta
(D_r + R^{-1}) beta + D_z gamma = 0
M_theta := int_0^inf R^2 <v>_Y dR = 0,   M_z := int_0^inf R <gamma>_Y dR = 0.   (8.2)
```

> **PROPOSITION 8.1 (mean momentum equations), p. 89.** *Let `u` have the decomposition (8.1),
> `p_m` its angularly invariant pressure correction, `Sigma_a = Q^{2A} T_{phys,a}` for
> `a = theta, z` the normalized base stress components from (5.41), and
> `Delta_0 = D_r^2 + R^{-1} D_r + D_z^2`. Retain the base's flat residual as a separate additive
> error. Then the angular mean of the remaining normalized residual is `(D_r p_m - g_r, E_theta, E_z)`
> where*
> ```
> E_theta = t_* v + (D_r + 2/R)(b v + beta V + beta v + W_{r theta})
>           + D_z(G v + V gamma + gamma v + W_{z theta}) - eps(Delta_0 - R^{-2}) v - (D_r + 2/R) Sigma_theta,
> E_z     = t_* gamma + (D_r + 1/R)(b gamma + beta G + beta gamma + W_{rz})
>           + D_z(2 G gamma + gamma^2 + W_{zz} + p_m) - eps Delta_0 gamma - (D_r + 1/R) Sigma_z,
> g_r     = -t_* beta - (D_r + 1/R)(2 b beta + beta^2 + W_{rr})
>           - D_z(b gamma + G beta + beta gamma + W_{zr}) + (2 V v + v^2 + W_{theta theta})/R
>           + eps(Delta_0 - R^{-2}) beta.                                              (8.3)
> ```

**[verified]** the three cylindrical transport identities displayed in the proof, and the
statement that the angular part of the vector Laplacian leaves `Delta_0 - R^{-2}` in the
radial and angular components and `Delta_0` in the axial component.

**Radial operators (8.4)-(8.6), pp. 89-90.**
```
I g(r,Y) = int_0^r g(r', Y + ((r')^{d_r} - r^{d_r}) v_r) dr'
J g(r,Y) = int_0^inf (same integrand) dr'
I_c g = I g - chi_m J g
D_e = D_r + e/R,   T_e f = R^{-e} I_c(R^e f),   A_e f = R^{-e}(d_R chi_m) J(R^e f).
```
Here `chi_m` is a fixed slow cutoff, zero on an inner collar, one on an outer collar.

> **LEMMA 8.2 (compactly supported radial primitive), p. 90.** *Fix a correction stage `j`, an
> exponent `alpha`, a band, its common torus `y = Y_{i_0}`, and `e in {0,1,2}`. Let
> `f in M_alpha` be a smooth scalar coefficient. Then `T_e f` is supported in the same active
> radial shell and in the projection of the input support to `(Z,T)`, belongs to `M_alpha`, and*
> ```
> D_e T_e f = f - A_e f.   (8.7)
> ```
> *If `int R^e <f>_Y dR = 0` at every `(Z,T)`, then `<A_e f>_Y = 0` and for every fixed `I` and
> integer `p >= 1`,*
> ```
> |d^I A_e f| <= C_{j,I,p} eps^{alpha + p kap_s} S_*^{B_{j,I,p}} zeta delta^{-B_{j,I,p}}.   (8.8)
> ```
> *If `f` is independent of `Y` and `int R^e f dR = 0` at every slow point, the cutoff remainder
> is identically zero.*

The proof uses the half-line representation (8.9), the Diophantine gain
```
|| L^{-p} H ||_{C_y^m} <= C_{m,p} || H ||_{C_y^{m+p+3}},   p >= 1   (8.10)
```
(`L = v_r . d_y`, `H` of zero Haar mean), and the exact identity
`J g = (-M^{-1})^p int_R d_U^p L^{-p} F^o(U+u, y + M u v_r) du`.
**[verified]** the derivative accounting in (8.10) and in (8.19): with
`|H^(k)| <= ||H||_{C^{m+p+3}}(1+|k|)^{-(m+p+3)}`, `m` derivatives, and a multiplier of size
`(1+|k|)^p`, the residual series is `sum_{k in Z^2}(1+|k|)^{-3} < inf`. The exponent 3 is
exactly what a 2-D lattice sum needs; the constant 4 in (8.19) is likewise exactly right.

**Mean patch and normalizations (8.11), p. 92.** `I_m = {sqrt(2X) : X in I_mean}` from
Theorem 4.6(vi); `rho_phys(r,z,t) = q^{-1/2} rho^(r/sqrt q)`, `int rho dR = 1`.
`T_0` reconstructs pressure, `T_1` an azimuthal vector potential, `T_1, T_2` recover radial
fluxes of axial and azimuthal momentum.

> **PROPOSITION 8.3, pp. 92-93.**
> *(i) Pressure. For a smooth shell-supported scalar source `g_r(R,Z,T,y)`, define*
> ```
> P = int <g_r>_Y dR,    p_m = T_0(g_r - rho P),    int rho dR = 1.   (8.12)
> Then  D_r p_m - g_r = -rho P - A_0(g_r - rho P),    d_R <p_m>_Y = <g_r>_Y - rho P.   (8.13)
> ```
> *If `g_r in M_alpha` then `P in S_alpha` and `p_m in M_alpha`; the cutoff remainder in (8.13)
> has zero auxiliary mean and is flat as `q -> 0` at every fixed derivative order.*
> *(ii) Axial velocity. Let `gamma_d(R,Z,T,y)` be a smooth shell-supported desired axial
> increment with `int R <gamma_d>_Y dR = 0` at every `(Z,T)`. With
> `Psi = r^{-1} I_c(r gamma_{d,phys})` from (8.4),*
> ```
> Psi_* = T_1 gamma_d,    Delta beta = -eps d_Z Psi_*,    Delta gamma = (D_r + 1/R) Psi_* = gamma_d - A_1 gamma_d.   (8.14)
> ```
> *The actual increment `(Delta beta, 0, Delta gamma)` is divergence-free and
> `<Delta gamma>_Y = <gamma_d>_Y`. If `gamma_d in M_alpha` then `Psi_*, Delta gamma in M_alpha`
> and `Delta beta in M_{alpha+1}`. If `gamma_d` is independent of `Y`, its zero weighted
> integral implies `Delta gamma = gamma_d` pointwise.*

**[verified]** the divergence identity `(rr + r^{-1})(-d_z Psi) + d_z((rr + r^{-1})Psi) = 0`
(this is the standard axisymmetric azimuthal-vector-potential representation) and the
factor `Q^{1/2-D} = Q^h = eps` in `Delta beta`.

**The three defects (8.12), (8.15), p. 94.**
```
J_theta = int R^2 <G v + V gamma + gamma v + W_{z theta}>_Y dR
J_z     = int R <2 G gamma + gamma^2 + W_{zz}>_Y dR - (1/2) int R^2 <g_r>_Y dR
c_rho   = (1/2) int R^2 rho dR.
```

> **PROPOSITION 8.4 (integrated equations), p. 94.** *Let `(beta,v,gamma)` and `w` be smooth
> shell-supported corrections of (8.1). Suppose the two moments (8.2) vanish at every `(Z,T)`
> and `p_m` is reconstructed from `g_r` by (8.12). Then*
> ```
> int R^2 <E_theta>_Y dR = eps d_Z J_theta,     int R <E_z>_Y dR = eps d_Z (J_z + c_rho P).   (8.16)
> ```

**[verified]** the two integration-by-parts identities
`int R^2(v_RR + R^{-1} v_R - R^{-2} v) dR = (2-1-1) int v dR = 0` and
`int R(gamma_RR + R^{-1} gamma_R) dR = 0`, and the pressure moment
`int R <p_m>_Y dR = -(1/2) int R^2 <g_r>_Y dR + c_rho P`.
Normalizations (8.17).

> **COROLLARY 8.5 (covariance targets), p. 95.** *Under the hypotheses of Proposition 8.4,
> with interior bump profiles `sigma^_theta, sigma^_z` (`int xi^2 sigma^_theta dxi = int xi sigma^_z dxi = 1`),*
> ```
> H_theta = -T_2(<E_theta>_Y - sigma_theta eps d_Z J_theta),
> H_z     = -T_1(<E_z>_Y - sigma_z eps d_Z (J_z + c_rho P)).   (8.18)
> ```
> *They are independent of the auxiliary variable, compactly supported in the active shell, and*
> ```
> (D_r + 2/R) H_theta = -<E_theta>_Y + sigma_theta eps d_Z J_theta,
> (D_r + 1/R) H_z     = -<E_z>_Y     + sigma_z     eps d_Z (J_z + c_rho P).
> ```
> *`H_theta, H_z` vanish whenever their integrands vanish for every `R`; the construction does
> not enlarge support in `Z`.*

> **LEMMA 8.6 (inverting the fast time derivative), p. 95.** *Let `N = v_t . d_y` on `T^2` and
> `F` smooth of zero Haar mean. Then `N phi = F` has a unique smooth zero-mean solution*
> ```
> N^{-1} F = sum_{k in Z^2 \ {0}} F^(k)/(2 pi i v_t . k) e^{2 pi i k . y},   ||N^{-1}F||_{C_y^m} <= C_m ||F||_{C_y^{m+4}}.   (8.19)
> ```
> *For `E_theta, E_z in M_alpha`, use this inverse in the original auxiliary coordinate `Y` to
> define the physical tangential update `-N_abs^{-1}(E_{theta,phys}^o, E_{z,phys}^o)`, with chart
> representatives*
> ```
> Delta v = -c_{i_0}^{-1} N_{i_0}^{-1} E_theta^o,   gamma_d = -c_{i_0}^{-1} N_{i_0}^{-1} E_z^o,   c_{i_0}^{-1} <= C S_*.   (8.20)
> ```
> *Both lie in `M_alpha` with zero auxiliary mean; realizing `gamma_d` by (8.14) gives*
> ```
> (Delta beta, Delta v, Delta gamma) = (-eps d_Z T_1 gamma_d, Delta v, gamma_d - A_1 gamma_d) in M_{alpha+1} x M_alpha x M_alpha,   (8.21)
> c_{i_0} N_{i_0} Delta v = -E_theta^o,    c_{i_0} N_{i_0} Delta gamma = -E_z^o + c_{i_0} N_{i_0} a,   a = -A_1 gamma_d.   (8.22)
> ```
> *The angular cancellation is exact and `c_{i_0} N_{i_0} a` is flat at every fixed order.*
> *Also `N_abs P_i = T_g^i P_i N_i`, `N_abs^{-1} P_i = T_g^{-i} P_i N_i^{-1}` (8.23).*

**Note (important, p. 96):** `N^{-1}` "acts at fixed `(R,Z,T)`, so it preserves slow and
radial support and their smooth zero extensions, **although it need not preserve
auxiliary-torus support**."

**Base power law on the mean patch (8.24), p. 97.** In `x = r/sqrt q in I_m`,
```
G_q = 0,    V_q = a(eta) x^{-1-2 lam},   |a(eta)| >= a_0 > 0,   a(eta) = 2^{1/2+lam} c_patch (1+eta^2)^{-1},  lam > 0.
```

> **LEMMA 8.7 (five-dimensional moment correction), p. 97.** *Suppose the fixed slow base
> satisfies (8.24) on the reserved interval `I_m`, with `lam > 0` and `|a(eta)| >= a_0 > 0`.
> There are three fixed azimuthal bump profiles and two fixed axial bump profiles supported in
> this interval with the following property. At each `(Z,T)`, for arbitrary scalar targets
> `(P, J_theta, J_z)`, their physical rescalings have a unique linear combination
> `(Delta v, gamma_d)` satisfying*
> ```
> int R^2 Delta v dR = 0,
> int R gamma_d dR = 0,
> int (2V/R) Delta v dR = -P,
> int R^2 (G Delta v + V gamma_d) dR = -J_theta,
> int (2 R G gamma_d - R V Delta v) dR = -J_z.        (8.25)
> ```
> *The increments are independent of the angular and auxiliary variables and depend linearly on
> the targets. For every `alpha`, targets in `S_alpha` give increments in `M_alpha` at every fixed
> derivative order, with constants permitted to depend on the fixed `lam`. The azimuthal-potential
> realization (8.14) has `Delta gamma = gamma_d` exactly and `Delta beta in M_{alpha+1}`. The two
> integral constraints are preserved exactly.*

Proof: on the patch, (8.24) splits the rows into an **angular block** with powers of `x`
equal to `2, -2-2lam, -2lam` and an **axial block** with powers `1, 1-2lam` — each list
consisting of distinct numbers for `lam > 0`. With `eta_j(x) = a_j^{-1} eta_0(x/a_j)`,
`a_j = e^{jd}`, one has `int x^p eta_j dx = mu_p e^{jdp}`, so after dividing rows by `mu_p`
each block is an ordinary Vandermonde matrix in the distinct numbers `e^{dp}`. Explicit
matrices on p. 98:
```
A_theta = [[mu_2, mu_2 e^{2d}, mu_2 e^{4d}],
           [2 a mu_{-2-2lam}, 2 a mu_{-2-2lam} e^{d(-2-2lam)}, 2 a mu_{-2-2lam} e^{2d(-2-2lam)}],
           [-a mu_{-2lam}, -a mu_{-2lam} e^{-2 d lam}, -a mu_{-2lam} e^{-4 d lam}]]
A_z     = [[mu_1, mu_1 e^{d}], [a mu_{1-2lam}, a mu_{1-2lam} e^{d(1-2lam)}]]
(u_0,u_1,u_2)^T = A_theta^{-1}(0, -P_q, -J_{z,q})^T,   (s_0,s_1)^T = A_z^{-1}(0, -J_{theta,q})^T.
```
"Since `G_q = 0`, these two blocks exhaust the five equations. ... The inverse may deteriorate
as `lam -> 0`; no uniformity in that limit is required."

**[verified]** all of it: the five power exponents; `int x^p eta_j dx = mu_p e^{jdp}`; the
Vandermonde structure and its nonvanishing determinant; the decoupling under `G_q = 0`
(with the *crossing* `A_theta <- (P, J_z)`, `A_z <- J_theta`); the row/target correspondence
in the displayed matrices.

> **LEMMA 8.8 (recomputation), pp. 98-99.** *With the base satisfying (8.24), let
> `(beta,v,gamma)` be a smooth divergence-free mean correction supported in the active shell,
> the base and wave field `w` fixed. Form `g_r` and `(P, J_theta, J_z)` by (8.3), (8.12),
> (8.15); apply Lemma 8.7 with those three defects and realize `gamma_d` by (8.14). Then with
> `(beta,v,gamma)_new = (beta+Delta beta, v+Delta v, gamma+Delta gamma)` and
> `R_g = g_{r,new} - g_r - (2V/R) Delta v`:*
> ```
> R_g = -t_* Delta beta - (D_r + 1/R)(2 b Delta beta + 2 beta Delta beta + (Delta beta)^2)
>       - D_z(b Delta gamma + G Delta beta + beta Delta gamma + gamma Delta beta + Delta beta Delta gamma)
>       + (2 v Delta v + (Delta v)^2)/R + eps(Delta_0 - R^{-2}) Delta beta,             (8.26)
> P_new       = int <R_g>_Y dR,
> (J_theta)_new = int R^2 <gamma Delta v + v Delta gamma + Delta gamma Delta v>_Y dR,
> (J_z)_new     = int R(2 gamma Delta gamma + (Delta gamma)^2)_Y dR - (1/2) int R^2 <R_g>_Y dR.   (8.27)
> ```
> *If `b` is bounded by `C_I eps S_*^{B_i}` and `V, G` by `C_I S_*^{B_i}` on the support of the
> azimuthal potential, and if `v, gamma in M_{0.9}`, `beta in M_{1.9}`, targets in `S_alpha`
> with `alpha >= 0.9`, then*
> ```
> R_g in M_{alpha + 0.9 - 2 kap_s},    (P_new, (J_theta)_new, (J_z)_new) in S_{alpha + 0.9 - 2 kap_s}.
> ```

**[verified]** the term-by-term class table on p. 99 (`-t_* Delta beta in M_{alpha+2}`,
`(D_r+1/R)(2b Delta beta) in M_{alpha+2-kap_s}`, `D_z(b Delta gamma + G Delta beta) in M_{alpha+2}`,
`2 v Delta v/R in M_{alpha+0.9}`, `(Delta v)^2/R in M_{2 alpha}`,
`eps(Delta_0 - R^{-2}) Delta beta in M_{alpha+2-2 kap_s}`) and the conclusion that the binding
exponent is `alpha + 0.9` from the `2 v Delta v/R` term, since `alpha >= 0.9`.

**Recomputation order (§8.7, p. 100).** (1) `W_ab = <w_a w_b>_theta`, then `g_r` from the third
row of (8.3); (2) `P = int <g_r>_Y dR`, `p_m = T_0(g_r - rho P)`; (3) `E_theta, E_z` from the
first two rows of (8.3); (4) `J_theta, J_z` from (8.15). Temporal update uses `E_theta^o, E_z^o`
in (8.20); moment update uses `(P, J_theta, J_z)` in (8.25).

---

## 2. Logical dependency chain

```
[EXTERNAL, Sections 3-5]
  Prop 5.5 (background u_B, p_B, T_phys; smooth one-sided extension at tau=0)
  Thm 4.6(iii) (strict stress direction / cone margin), 4.6(vi) (reserved mean patch I_mean)
  Lemma 4.5 (admissible stress cone), (4.27) leading-stress bounds, (4.30) base power law
  (5.18), (5.41) residual stress, (5.42) b = O(eps), (5.44) summed base
        |
        v
SECTION 6  (geometry + bookkeeping)
  (6.1) charts -> (6.2) J_g, v_r, v_t, rho_g, d_r -> (6.3) evaluation map -> (6.4),(6.6) operators
  (6.2) --> (6.7) Diophantine bounds ......................... used ONLY in (8.10) and (8.19)
  (6.5) covering index i(l) --> (6.14) bounded Delta --> Lemma 6.2 (common torus)
  (6.8),(6.9) labels/slow cutoffs + (6.10),(6.11) rectangles --> LEMMA 6.1 (separation)
  (6.8),(6.9) mesh S_*^{-3} --> LEMMA 6.3 (label counting)
  (6.16) time cutoff psi
  (6.21),(6.22) path along a rectangle  [needed by Prop 7.2 Step 3]
  Defs 6.4, 6.5 (M_a, W_a, S_a; weights zeta, sqrt(zeta), P_v) --> PROP 6.6 (algebra (6.30)-(6.32))
        |
        v
SECTION 7  (waves)
  (7.1) cone  <-- Lemma 4.5 / Thm 4.6(iii)
  (7.2) k, B_s, s(v);  (7.3),(7.4) phase and n_Phi
  LEMMA 7.1  <-- (7.1),(7.2),(7.3), (6.6), Lemma 6.1 (rectangles), slow-box diameter O(S_*^{-3})
      |         gives (7.9) phase defect, (7.10) diagonalization, (7.11) damping
      +--> (7.12) envelope P(v) --> (7.16) Gaussian bounds
      +--> PROP 7.2 (zero-data inverse) <-- (7.17)-(7.20), Lemma 6.2, (6.21)
      |        +--> COR 7.3 (uniform q_* over stages and harmonic sets)  <-- (7.18)
      +--> LEMMA 7.4 (homogeneous pulse)  <-- (7.10),(7.16),(7.19)
               +--> (7.22) energy transfer (why both T_N and T_K are needed)
               +--> PROP 7.5 (covariance match)  <-- LEMMA 6.1 (disjointness), (7.1), (4.27)
                        +--> (7.30) global leading stress
                        +--> PROP 7.6 (linearized L)  <-- Prop 6.6
                                 +--> (7.35),(7.36) global W_0^as, L^as
  LEMMA 7.7 (exact curls)  <-- Prop 6.6, Lemma 7.1
  (7.40) cutoff residual  <-- Prop 7.2, (6.16), (7.16)
  COR 7.8 (covariance of actual velocities)  <-- Lemma 7.7, Prop 7.6, Prop 6.6
        |
        v
SECTION 8  (means)
  (8.1),(8.2) decomposition + two moments;  (5.41) --> PROP 8.1 (mean system (8.3))
  (6.7) + (8.10) --> LEMMA 8.2 (radial primitive T_e, cutoff remainder A_e, flatness (8.8))
      +--> PROP 8.3(i) pressure T_0     ] <-- Thm 4.6(vi) for the mean patch
      +--> PROP 8.3(ii) axial velocity T_1 ]
      +--> PROP 8.4 (integrated defects (8.16))  <-- Prop 8.1, Prop 8.3, (8.2)
             +--> COR 8.5 (covariance targets H_theta, H_z via T_2, T_1)
                    ... these are the Sigma fed to PROP 7.6 / COR 7.8
  (6.7) --> LEMMA 8.6 (N^{-1}) --> (8.20),(8.21),(8.22)   <-- Prop 8.3(ii) for realization
  (8.24) [<-- (4.30), (5.44), Prop 5.5] --> LEMMA 8.7 (five-equation system (8.25))
      +--> LEMMA 8.8 (exact recomputation (8.26),(8.27))  <-- Prop 8.1, Lemma 8.7
  §8.7 recomputation order
        |
        v
[EXTERNAL, Section 9]  Prop 9.1 (transport/viscous residual), PROP 9.6 (the cycle and
  its eps-gain), Lemma 9.7 (one domain for all stages), Prop 9.9 (u_loc) --> Thm 3.1 --> Thm 1.1
```

**The load-bearing coupling** is: Corollary 8.5 produces two-vector targets `Sigma`
independent of the angular and auxiliary variables; Proposition 7.6 requires exactly that
independence ("This auxiliary independence is an additional requirement here; the class
`M_alpha` alone allows auxiliary dependence", p. 84); Lemma 8.6 supplies the piece with zero
auxiliary average that Corollary 8.5 cannot handle. The three-way division of the mean
residual (compactly supported defects -> Lemma 8.7; auxiliary-averaged remainder ->
Corollary 8.5 -> waves; zero-auxiliary-average part -> Lemma 8.6) is the architectural
crux of Section 8 and looks complete on its own terms.

---

## 3. Load-bearing estimates; what is under-justified

### 3.1 The tightest quantitative balances (these are where the paper lives or dies)

**(L1) `|E| <= C/S_*` in (7.10) against `L_s ~ S_*`.** The frame error is integrated over the
whole pulse, of length `L_s = 2 r_0/c_i ~ S_*`. So `int_0^{L_s} |E| dv = O(1)` — a *bounded
multiplicative* distortion of the envelope, and no better. If the bound in (7.10) were
`C S_*^{-1+kap}` for any `kap > 0` the argument would collapse: the envelope `P` would be off
by `e^{c S_*^{kap}}`, destroying both the two-sided bound (7.16) and (through it) the
"exponentially small tails" claim. **This is the single sharpest inequality in Section 7,
and it is exactly saturated.** The paper does derive it (p. 76: slow box diameter
`O(S_*^{-3})`, base comparison `O(eps^2)`, rounding `O(k^{-1})`, then multiplication by
`v = O(S_*)`), and I **[verified]** each of the three contributions is `<= C/S_*` for
`l >= l_0(h)`: `S_*^{-3} S_* = S_*^{-2}`, `eps^2 S_* = 2^{-2hl} l^2`, `k^{-1} S_* = 2^{-hl/2} l^2`.
But note the consequence: **`l_0` must be chosen depending on `h`, and the paper never
displays that dependence.** It also means the frame comparison must be done with the
*frozen* `N, K, c_0`, which is why the slow mesh is as fine as `S_*^{-3}` (and hence why
Lemma 6.3's counts are `S_*^9`).

**(L2) `eps k^2 ~ 1` (p. 74, p. 81).** `k = ceil(eps^{-1/2})`, so viscosity is *exactly
marginal* in the pulse equation: the oscillations sit at Reynolds number `O(1)`. This is not
an inviscid construction with viscosity treated as a perturbation. Everything therefore
depends on the algebraic fact that shear amplification can still beat `O(1)` viscous
damping over a window — encoded in `a_net(y) > 0` for `|y| < u_*`. **[verified]** exactly:
with `lam(v) = lam_0/sqrt(1+s^2)` from (7.10) and
`d_ref = eps k^2 B_s^2 (1+s^2) = lam_0(1+s^2)/(1+u_*^2)^{3/2}` from (7.2),
```
a_net(y) = lam(y) - d_ref(y) = lam_0/sqrt(1+y^2) - lam_0(1+y^2)/(1+u_*^2)^{3/2},   a_net(u_*) = 0,
```
and the displayed `d/dv` formula on p. 78 is exactly right. So `B_s` is *defined* to place
the crossover at `|s| = u_*` — this is engineered, not lucky. The pulse grows for
`v < L_s/2` and decays for `v > L_s/2`, with amplification factor
`P(L_s/2)/P(0) = e^{Theta(S_*)} = e^{Theta(l^2)}`.

**(L3) The `S_*^{-1/2}` concentration loss (7.28), (7.29).** The pulse is "on" only in a
Gaussian window of width `sqrt(L_s)` out of a rectangle of length `L_s`, so
`h_sigma ~ c_i sqrt(L_s) ~ S_*^{-1/2}`, hence `|det H| >= c/S_*`,
`||H^{-1}|| <= C sqrt(S_*)`, `y_sigma ~ sqrt(S_*)|T|`, and the wave amplitudes carry a
factor `S_*^{1/4}`. **[verified]** all of these, including `h_sigma ~ c_i sqrt(L_s)` from
`c sqrt(L_s) <= int psi^2 x^2 dv <= C sqrt(L_s)` and `c_i L_s = 2 r_0`. Two consequences the
reader must notice:
- the "frozen `s = sigma u_*`" approximation in (7.28) is good only to `O(S_*^{-1/2})`,
  because `s(v) - sigma u_* = O(u_* sqrt(L_s)/L_s)`. Hence **the strict cone margin
  `eta_c` from Theorem 4.6(iii) must exceed `C S_*^{-1/2}`** — another `l_0` condition,
  this time depending on an external quantity from Section 4. Stated only as "the errors in
  (7.28) preserve this strict inequality" (p. 83).
- `psi = 1` must cover the Gaussian core, i.e. `sqrt(L_s) << L_s/5`. True for `l >= l_0`.

**(L4) The cone inequality and the choice of `u_*`.** **[verified]** The realizable cone is
strictly smaller than the cone in (7.1): the two-column decomposition
`h_+ y_+ = (1/2)(-T_N/A_c - T_K/u_*)`, `h_- y_- = (1/2)(-T_N/A_c + T_K/u_*)` is positive iff
`|T_K| < (u_*/A_c)(-T_N)`, i.e. iff
```
|c_0 T_K / T_N| < u_*/sqrt(1 + u_*^2)  < 1,
```
whereas (7.1) only gives `|c_0 T_K/T_N| < 1`. The paper closes the gap by taking `u_*` large
enough (p. 74) — legitimate, since the supremum over the closed annulus is `< 1` strictly.
**[verified]** `A_c = -c_0 sqrt(1+u_*^2)`, and I reconstructed the reference-cone inequality
independently. I also **[verified]** the internal consistency of the frame constants:
`c_0^2 = (v_s - 2)/2` follows exactly from `g_0 = F_0(-a, b_s)` and
`lam_0^2 = 2 a F_0^2 (1 - 2/v_s)` **provided** `v_s = (a^2 + b_s^2)/a`; the factors `F_0`
cancel. I could not check that last identification (it lives in (4.11), (4.20)) —
**[unverified, external]**. Equivalently `lam_0^2 > 0  <==>  v_s > 2  <==>  c_0^2 > 0`.

**(L5) "Cutoff tails are exponentially small" (7.16) + (6.16) + p. 87.** **[verified] and I
consider this claim correct and quantitatively robust.** Where `psi' != 0` or `psi != 1`,
`|v - L_s/2| >= L_s/5`, so `P <= e^{-c L_s/25} = e^{-c' S_*} = e^{-c' l^2}`. The conversion to
physical variables costs `q^{-N} Q^{-M} S_*^C = 2^{(M+N)l} l^{2C}`, and
`exp(-c' l^2 + (M+N) l log 2 + 2C log l) -> 0`, so *for each fixed `N`* the physical
derivatives are `O(q^N)`. This works because `S_* = l^2` grows *superlinearly* in `l` while
the price is only exponential-linear in `l`. This is the one place where the peculiar choice
`S_* = l^2` (rather than, say, `S_* = l`) is essential and visibly so. Caveats:
- the constant `c` in the upper bound of (7.16) is `inf|d a_net/dv| L_s`-derived and depends
  on `u_*`, `lam_0`, hence on the profile data of Section 4. **[unverified]**
- the required `l_0` and the constant `C_N` depend on `N` (and on `M`, `C`). So the
  residual is flat to *each* finite order with constants depending on the order; infinite
  flatness of the assembled field is a *diagonal* statement deferred to Section 9.
- derivatives of `psi` cost only `L_s^{-1} ~ S_*^{-1}` — harmless. **[verified]**

**(L6) `A_wave ~ q^{-1/2-h/2}` versus the background `~ q^{-1/2-h}`.** **[verified]** The
chart-to-physical factor is `Q^{-A}` with `A = 1/2+h`, and `W_0 in W_{1/2}` means chart
amplitude `~ eps^{1/2} = Q^{h/2}`, so physical wave amplitude
`Q^{-1/2-h} Q^{h/2} = Q^{-1/2-h/2}`, i.e. **the waves are smaller than the background by the
factor `eps^{1/2} = Q^{h/2}`** (recall `Q -> 0`, so `Q^{-1/2-h} > Q^{-1/2-h/2}`). Their
quadratic covariance is then `~ eps` times the background's quadratic scale, which is exactly
the size of the target annular stress: **[verified]** the scale identity on p. 84,
```
Q^{-2A} eps T_{0,*} = Q^{-2A+h}(Q/q)^{A+1/2} T_0 = q^{-A-1/2} T_0,
```
where the `Q`-exponent cancels identically (`-2A + h + A + 1/2 = 0`). So the amplitude scale
is forced, in the standard convex-integration way, as the square root of the stress. I also
**[verified]** the consistency of the wavelength: `sqrt(Q)/k = Q^{(1+h)/2}`, and
amplitude x wavelength `= Q^{-1/2-h/2} Q^{(1+h)/2} = 1` — exactly the statement "the powers
of `Q` in their product cancel", i.e. the oscillations sit at unit Reynolds number.
**The relation of `A_wave` to the background is therefore internally consistent and is a
genuine smallness (`eps^{1/2}`), not a fudge.** But note the "up to powers of `S_*`" rider on
p. 81: from (L3) the actual factor is `eps^{1/2} S_*^{1/4}`, i.e.
`Q^{-1/2-h/2} l^{1/2}` — the paper is honest about this but the reader should carry it.

**(L7) The `kap_s = 10^{-5}` budget.** Each `D_r` costs one factor `eps^{-kap_s}` (6.32). All
uses in my range are finite in number per stage, so the total cost is `eps^{-C kap_s}` with
`C` small: (7.42) loses `kap_s` on the curl remainder and one more on the radial divergence;
Lemma 8.8 loses `2 kap_s`. Since the claimed per-cycle gain is `0.9` (Lemma 8.8) or
`min(rho - 1/2, 1/2 - kap_s, alpha - 1)` (Corollary 7.8), the `kap_s` budget is comfortable.
**[verified]** arithmetic of the class table on p. 99. But the number of radial derivatives
per *derivative order* `|I|` grows with `|I|` — so `eps^{-|I| kap_s}` losses appear at high
derivative order. This is fine for a fixed order but is another place where "uniform in
derivative order" cannot be literally true; the paper's `C_{j,I}, b_{j,I}, d_{j,I}` formalism
accommodates it.

### 3.2 The separation-of-supports device: does `Y(r,t)` reintroduce interactions?

**My conclusion: no, and the design is correct — but for a reason the paper states only in
passing, and which a reader should make explicit.**

The separation in (6.13) is asserted for the **absolute** rectangles
`R_gam^{abs,+} = pi_{i(l)}^{-1}(R_gam^+)`, i.e. for *all* deck copies on `T^2`, not merely for
the band representatives. Physical evaluation (6.3) sends each `(r,t)` to a *single* point
`Y in T^2`. Two fields whose absolute-torus supports are disjoint therefore have literally
disjoint physical supports, and their pointwise product vanishes identically. **Evaluation
cannot reintroduce cross terms.** This is why the lemma is stated on the absolute torus and
why Lemma 6.2 (Haar averages agree on absolute / common / band tori for pullbacks) is
needed.

What this *means* physically, and what the paper never says in plain words: the two wave
families are not "magically non-interacting"; they occupy **disjoint, finely interleaved
physical space-time slabs**. Since `t_* v = 1` with `t_* = -eps d_T + c_i N_i` and
`c_i ~ S_*^{-1}`, the auxiliary phase sweeps a full rectangle in physical time of order
`c_i^{-1} r_0 ~ S_*`, in normalized units — so pulses that "overlap" at the coarse slow
scale in fact fire in alternating fast time windows. There is no free lunch: the price of
disjointness is that each family is on only a fraction of the time, and that price is
precisely the `S_*^{-1/2}` in `h_sigma` (L3), paid back by boosting amplitudes by `S_*^{1/4}`.
I think this is a legitimate and rather elegant device; a reader should nonetheless confirm
independently that the two families' Reynolds stresses genuinely *add* (rather than the
time-sharing halving the effective stress in a way not accounted for) — the accounting is
in (7.27), where the `c_i int_0^{L_s} ... dv` factor is exactly the duty cycle, and I
**[verified]** the normalization `c_i L_s = 2 r_0` closes.

**Where I would push back on the write-up:**

1. **(U1) Lemma 6.1, Step 1, "bounded degree" (p. 66) is a three-sentence argument for a
   uniform-in-band combinatorial claim.** The pieces are checkable — `|l - l'| <= 2` from
   `q/Q, q/Q' in [1/2, 2]`, mesh-ratio `(l'/l)^6 <= (1 + 4/l_0)^6` for `|l-l'| <= 4`, and
   `Q/Q'` bounded — but the paper compresses "the changes of physical scale in each
   coordinate are bounded, as are the ratios `l^2/(l')^2`" into one clause covering three
   different coordinate rescalings (`Q^{1/2}, Q^D, Q`). **Terse but I believe correct.**
2. **(U2) The greedy coloring of a countable graph (p. 66) is invoked without a well-order
   or a bound on the number of colors.** A bounded-degree graph is `(D+1)`-colorable
   greedily, but the paper needs the *color set to be finite and fixed before `r_0` is
   chosen* (because `r_0` is then chosen small relative to the finitely many separations in
   (6.15)). This ordering of choices is correct as written but is easy to get wrong and
   deserves a sentence.
3. **(U3) "Each excluded equality is a proper closed constraint on the finite tuple of
   centers. ... Their complement is open and dense and contains a rational tuple." (p. 66)**
   For `Delta = 0`, `mu != nu`, the constraint is `c_mu != c_nu` — fine. For `Delta > 0`,
   `mu = nu` it uses `det(J_g^Delta - I) != 0`; **[verified]** since the eigenvalues
   `(4 +- sqrt 2)^Delta != 1`. For `Delta > 0`, `mu != nu` the constraint set is a coset of a
   finite subgroup, again proper and closed. So the claim is right, but the paper does not
   distinguish these three cases and the reader has to.
4. **(U4) The `14^{i(l)}` filaments.** A band-`l` wave's absolute support is a union of
   `14^{i(l)}` rectangles, each of side `~ r_0 Lam_g^{-i}` by `~ r_0 T_g^{-i}` — total measure
   `~ r_0^2`, constant, but the sets become long thin near-dense filaments as `i -> inf`.
   Disjointness of two such filament families is a genuinely delicate statement, and the
   proof handles it correctly *only because* `Delta = |i(l) - i(l')| <= Delta_max` (6.14) for
   the pairs that matter; the error term `J_g^Delta e_nu` is then bounded by
   `||J_g||^{Delta_max} C r_0`. **[verified]** the logic. This is the crux of Lemma 6.1 and,
   to the paper's credit, it is where the `T_g` / `Lam_g` split of `J_g` earns its keep.
   I would still call the presentation *compressed relative to its importance*.
5. **(U5) "The same identity holds for derivatives of smoothly extended fields and after
   evaluation on (6.3)." (p. 66)** For derivatives this is trivial (support does not grow),
   and the `2 r_0` enlargement provides room, but the phrase "and after evaluation" is doing
   the essential work described above and is stated as an aside.
6. **(U6, a real strength, not a gap)** The paper is careful that the disjointness applies
   only to *cut-off* fields: "Before that cutoff, the coefficient algebra is used only on its
   labelled band rectangle" (p. 72), because Proposition 7.2's solution `t_m` does *not*
   vanish outside the rectangle (Definition 6.5: "The continuation beyond `v = 0, L_s` has no
   envelope bound"). Multiplication by `psi` (whose support is strictly inside) restores
   containment. This is the kind of detail whose omission would have been fatal, and it is
   handled. **[verified]** that `supp psi ⊂ {|v - L_s/2| < L_s/3}` sits strictly inside
   `v in (0, L_s)`, so `psi t_m` and `psi' t_m` are both supported inside the rectangle.
7. **(U7) `N^{-1}` breaks the localization.** Lemma 8.6, p. 96: the multiplier "need not
   preserve auxiliary-torus support." So the mean corrections `Delta v, gamma_d` are *not*
   confined to rectangles. This is fine (they are means, and `M_a W_b ⊂ W_{a+b}` controls
   their products with waves), but it means **the separation device protects only wave x wave
   products; all other interactions are controlled purely by `eps` powers.** A reader should
   check that no step silently assumes disjointness for a field that has passed through
   `N^{-1}`, `T_e`, or `A_e`. I did not find such a step in my range, but I could not check
   Section 9.

### 3.3 Does amplification-then-damping control **all** derivatives?

**Partially — the structure is right, the argument for arbitrary order is an outline.**

What is proved carefully: (7.19) `||V_m(v,w)|| <= C P(v)/P(w)` for values, from
`d_v |z| <= (lam - d_ref + C/S_*)|z|` plus `L_s/S_*` bounded; and (7.18) makes every higher
harmonic at least as damped as `m = 1`. **[verified]** both, including that `exp(C L_s/S_*)`
is `O(1)` — which is (L1) again.

What is outlined: the derivative bounds (7.14), (7.15). The mechanism is (7.20),
```
(D^I z_m)' = A_m D^I z_m + D^I g_m + sum_{0 < J <= I} binom(I,J) (D^J A_m) D^{I-J} z_m,
```
solved with the same propagator, each order raising the polynomial degree in `S_*` by
`1 + max{b_I, max_{0<J<=I}(c_J + b'_{I-J})}` and each integration contributing `L_s = O(S_*)`.
The key hypothesis is `|D^I A_m| <= C_I S_*^{c_I}` — polynomial, not small. Concerns:
- **(U8) the polynomial degrees grow with `|I|`, and the paper never displays the growth.**
  Concretely: `t_sigma^h(0)` has `z_+(0) = P(0) = exp(-Theta(S_*))`, so each slow derivative
  of the initial datum costs a factor `~ S_* x S_*^{c}` (differentiating the exponent
  `int(lam - d_ref)`, itself of size `S_*`). So `b_I` grows at least linearly in `|I|`. That
  is *allowed* by Definition 6.5 (`b_{j,I}` may depend on `I`), and harmless because
  `eps^{alpha}` beats every fixed power of `S_*`; but it means the `l_0` needed for a given
  smallness *grows with the derivative order*, so no estimate in Section 7 is uniform in
  `|I|`. Lemma 7.1's phrasing — `q_*` "common to all labels **and fixed derivative orders**"
  — is stronger than what is argued and should be read as "the smallness conditions that fix
  `q_*` are finitely many and involve no derivative", which is what p. 77 actually proves
  (`S_*^2(eps + eps^2 + k^{-1}) <= 1` — **[verified]** this holds for `l >= l_0(h)`, since
  `l^4 2^{-hl/2} -> 0`). A reader should confirm that no later step needs uniformity in `I`
  of `q_*` itself.
- **(U9) the induction in (7.20) is asserted, not carried out.** "An induction using (7.19)
  increases the polynomial degree by at most ..." — the statement is plausible and standard,
  but it is the *only* justification offered for an infinite family of estimates that the
  whole `C^infinity` claim rests on. I cannot certify it from the text; I can say the
  structure (linear ODE, bounded propagator, polynomially bounded coefficient derivatives,
  interval length `O(S_*)`) is the standard one and I see no obstruction.
- **(U10) Lemma 7.4's derivative claim** ("every fixed derivative ... bounded by
  `C_I S_*^{b_I} P(v)`") is proved by "Differentiating the coordinate equation and using
  (7.19) proves the derivative estimates just as in (7.20)". Same status as (U9). Note that
  `v`-derivatives are *easy* (from the ODE, and `ds/dv = u_*/L_s` is small); the slow
  derivatives are the substantive ones.
- **(U11) The one thing I would call a likely error.** In (7.25), p. 82, the second bound is
  printed as
  ```
  |D^I y_sigma| <= C_I S_*^{b_I} delta^{-a_I}
  ```
  with **no factor `zeta`**, while the first is `y_sigma >= c sqrt(S_*) zeta`. The very next
  page uses the chain rule
  ```
  D^I a_sigma = sum_{n=1}^{|I|} sum c_{I_1..I_n} y_sigma^{1/2 - n} prod_{v=1}^n D^{I_v} y_sigma
  ```
  and concludes `|D^I a_sigma| <= C_I S_*^{b_I} sqrt(zeta) delta^{-a_I}`. With (7.25) as
  printed, the `n`-th term is `O((sqrt(S_*) zeta)^{1/2 - n} (S_*^{b})^n)`, i.e. it carries
  `zeta^{1/2 - n}`, which **diverges** at the shell edges for `n >= 1`. The conclusion
  requires `|D^I y_sigma| <= C_I S_*^{b_I} zeta delta^{-a_I}`. And the proof of
  Proposition 7.6 confirms this reading verbatim: "Every derivative of `y_sigma^{-1/2}` is a
  sum of terms of the form `C y_sigma^{-k-1/2} prod D^{I_v} y_sigma`. **The `k` numerator
  weights cancel all but `zeta^{-1/2}` in the denominator**" (p. 84) — that cancellation is
  possible only if each `D y_sigma` carries a `zeta`. So I read (7.25) as a **typo: a missing
  factor `zeta`** (which is also consistent with `y_sigma <= C sqrt(S_*)|T|` in (7.29) and
  with `|T| <~ zeta` from (4.27)). Harmless if intended, but it is the kind of omission that
  makes the statement as printed false, and it should be corrected. **[flagged; the
  correction is available from the paper's own proof text.]**

### 3.4 Other uniformity assertions I could not verify

- **(U12) "constants ... are uniform in band, label, rectangle copy, and point" (p. 69, and
  repeated on pp. 73, 78, 90, 97).** This is the single most repeated and least
  demonstrated claim in my range. It is *plausible* because every construction is expressed
  in normalized chart variables with `Q, eps, S_*` treated as constants and all
  band-dependence funnelled through the bounded matrices of Lemma 6.2 and the bounded
  `Delta <= Delta_max`. But it is never proved by a uniform argument; it is asserted
  per-construction.
- **(U13) The number of pulses.** Uniformity in the number of labels is handled by
  Lemma 6.3 (pointwise `O(1)`, integrated `O(S_*^9)`) plus "the sum over these labels
  therefore permits polynomial growth in `S_*`" (p. 68) and "the common-torus lifts adds at
  most the fixed factor `14^{Delta_max}`". I **[verified]** the counting arithmetic
  (`O(S_*^3)` positions per coordinate, three coordinates, two signs). What is *not* shown
  in my range is that the accumulated polynomial factors are beaten by the per-cycle gain —
  again Proposition 9.6.
- **(U14) Stage-uniformity.** Corollary 7.3 is the one place where stage-independence is
  actually *proved* (via the exact damping identity (7.18)), and it is a good argument.
  But `C_{j,I}` are allowed to depend on the stage `j` throughout Definitions 6.4-6.5, so
  the per-stage constants may grow arbitrarily; only Section 9 can control that.
- **(U15) Lemma 8.7's dependence on `lam` and on `G_q = 0`.** The invertibility of the 5x5
  system relies on **two** external facts: `lam > 0` (so the five power-law exponents are
  distinct) and `G_q = 0` on the mean patch (so the system *decouples* into 3x3 and 2x2
  blocks). Both come from (8.24) / (4.30) / (5.44) / Proposition 5.5 — **[external,
  unverified]**. If `G_q` were merely small rather than zero, the decoupling argument as
  written would not apply and invertibility would need a perturbation argument. The paper
  concedes "The inverse may deteriorate as `lam -> 0`; no uniformity in that limit is
  required" — fine, `lam` is fixed, but it means the constants of Lemma 8.7 (and hence of
  the whole mean correction) depend on Section 4's profile data in an undisplayed way.
- **(U16) Lemma 8.2's flatness order.** (8.8) holds "for every ... integer `p >= 1`" with
  constants `C_{j,I,p}` depending on `p`. Any *fixed* flatness order `N` is achieved by
  choosing `p` with `h(alpha + p kap_s) > N + L`; the paper is explicit that "A fixed finite
  regularity bound would give only a fixed finite flatness order" (p. 92). So infinite
  flatness needs a diagonal argument — deferred. **This is honestly flagged by the authors.**
- **(U17) No net gain in my range.** Proposition 7.2 explicitly "proves every fixed-order
  derivative estimate with **no decrease of `alpha`**" (p. 79); Lemma 8.6's increments
  "retain their `eps` exponent" (p. 96). The only gains visible in Sections 6-8 are
  structural: Corollary 7.8's remainders in `M_{rho + alpha - 1/2}`, `M_{alpha+1/2-kap_s}`,
  `M_{2 alpha - 1}` (gains `rho - 1/2`, `1/2 - kap_s`, `alpha - 1`) and Lemma 8.8's
  `M_{alpha + 0.9 - 2 kap_s}` (gain `0.9 - 2 kap_s`). I **[verified]** that these are
  genuine gains **provided** `rho > 1/2` and `alpha > 1`, i.e. provided the accumulated
  correction is already smaller than the primary wave and the residual already smaller than
  `eps`. Whether the iteration can be started and closed with those inequalities is
  Proposition 9.6's business. **The reader should not mistake Sections 6-8 for a convergence
  proof; they are a toolbox plus an exact bookkeeping identity.**

---

## 4. Plain-language account of the mechanism

**The target.** By the time Section 6 begins, Sections 4-5 have built a smooth axisymmetric
background that concentrates self-similarly at `t = 1`, and have computed that its momentum
equation fails to close: there is a leftover stress `T` supported on a thin annulus around
the collapsing vortex core. `T` is a *radial flux of azimuthal and axial momentum* —
physically, the background needs momentum transported across cylinders at a rate it cannot
supply by itself. Sections 6-8 supply that transport with a bath of tiny, very fast
oscillations whose *time-averaged* Reynolds stress equals `T`.

**How shear amplifies the pulses.** This is the classical Orr / Kelvin transient-growth
mechanism, engineered. Put a plane wave `e^{i k Phi}` into a background with radial shear
`g = (r F_r, G_r)`. The shear tilts the wave's phase gradient: in (7.3) the phase carries a
term `-v(p F + p_z G)`, so as the "pulse clock" `v` advances the radial component of
`n_Phi = grad Phi` sweeps linearly, `n_{Phi,r} ~ B_s s(v)` with
`s(v) = sigma(u_*/2 + u_* v/L_s)` (7.2). Two things happen simultaneously:
- **Growth.** While the wavevector leans *against* the shear, the shear does work on the
  wave: (7.22) says `(1/2) d|t|^2/dv = -g . (t_r(t_theta, t_z)) - eps k^2 |n_Phi|^2 |t|^2`.
  The first term is exactly the shear feeding the wave, at rate
  `lam(v) = lam_0/sqrt(1 + s(v)^2)` — largest when the wavevector is most nearly radial and
  falling off as it tilts.
- **Damping.** The same tilting makes `|n_Phi|` grow, and viscosity costs
  `d = eps k^2 |n_Phi|^2 ~ d_ref = lam_0(1 + s^2)/(1 + u_*^2)^{3/2}` — *growing* as the wave
  tilts.

Growth wins early and loses late. The crossover is placed *by construction*: `B_s` in (7.2)
is chosen so that `lam = d_ref` exactly when `|s| = u_*`, i.e. exactly at the midpoint
`v = L_s/2`. The result is the envelope `P(v)` of (7.12) with the two-sided Gaussian bound
(7.16): the amplitude rises by a factor `e^{c l^2}` from `v = 0` to the midpoint, then falls
back by the same factor. That is what "a pulse" means in this paper.

This shape is not an aesthetic choice, it is what makes the construction *compactly
supported in time*. Because the pulse is exponentially tiny at both ends of its window, you
may multiply it by a smooth time cutoff `psi` (6.16) and the equation you actually solved is
violated only by `(1 - psi)f + psi' t` (7.40) — a quantity of size `e^{-c l^2}`, which beats
every power of `q` (p. 87). So the pulses can be switched on and off cleanly. Without the
"grow then decay" structure you would have to either accept a non-compactly-supported wave
or a cutoff error comparable to the wave itself.

**Why the Reynolds stress can be prescribed.** A real oscillation `a b_sigma` contributes to
the averaged momentum flux at order `a^2`. Averaging over the angle and over the auxiliary
torus, the flux contributed by family `sigma` is `a^2 H_sigma`, with `H_sigma` given
explicitly by (7.27)-(7.28) as
```
H_sigma = h_sigma(-A_c N - sigma u_* K + small),    h_sigma > 0.
```
Because the *direction* `H_sigma/|H_sigma|` is determined by the frozen frame `(N, K)` and by
`u_*`, and the *magnitude* is a free nonnegative number `a^2`, the set of achievable stresses
is a cone. Setting `y = H^{-1} T_{0,*}` and `a_sigma = sqrt(y_sigma)` (7.24) inverts the map
exactly, giving `C(W_0) = eps T_{0,*}` (7.26) — the averaged quadratic flux equals the
prescribed stress, on the nose. Summing over the squared partition of unity gives the global
identity (7.30). That is the sense in which "the Reynolds stress can be prescribed": it is
literally a 2x2 linear solve, and the only issue is whether the solution has *nonnegative*
components.

**Why two families are needed.** Three reasons, all the same reason:
1. **Squares are nonnegative.** A real wave contributes `a^2 H_sigma` with `a^2 >= 0`. One
   family gives a *ray*, not a plane. The target `T = T_N N + T_K K` is a general two-vector
   in the tangential plane, so you need at least two independent rays whose nonnegative span
   contains it.
2. **Both components of the momentum flux matter.** (7.22) shows the shear controls only the
   `N`-component: the energy transfer is `-|g_0| T_N`. The `K`-component `T_K` is transverse
   to the shear direction and is not fixed by the energetics — "the covariance construction
   also prescribes the transverse component `T_K`. We therefore need both components of the
   covariance to match the prescribed momentum flux" (p. 81). The two families supply it: the
   sign `sigma` flips the `-sigma u_* K` term in (7.28), so the `+` family transports `+K`
   momentum and the `-` family `-K`. Physically the two families lean opposite ways against
   the shear.
3. **The cone must be strict.** Solving gives
   `h_+ y_+ = (1/2)(-T_N/A_c - T_K/u_*)`, `h_- y_- = (1/2)(-T_N/A_c + T_K/u_*)`, both positive
   iff `T_N < 0` and `|T_K| < (u_*/A_c)(-T_N)`. That is why the background construction had to
   deliver a stress in a *strictly interior* cone (7.1) — a boundary stress would need a
   vanishing amplitude, and the amplitude appears under a square root, so all derivative
   bounds would degenerate. `u_*` is then chosen large enough that
   `u_*/sqrt(1 + u_*^2)` exceeds the (strictly sub-unit) cone ratio.

**Why the auxiliary torus.** Two pulses whose slow envelopes overlap would, in the naive
construction, produce a cross Reynolds stress `~ a_1 a_2` of the *same size* as the intended
diagonal terms, and pointing in an uncontrolled direction. The fix is to give each label a
small rectangle on an auxiliary 2-torus and to make the rectangles disjoint (Lemma 6.1);
since the physical fields are obtained by evaluating at the single point
`Y(r,t) = v_r r^{d_r} + v_t t`, disjoint auxiliary supports means disjoint physical supports,
and all cross products vanish identically. The arithmetic of `J_g = [[3,1],[1,5]]` (an
integer matrix with irrational eigen-directions and eigenvalues `4 +- sqrt 2`) is what makes
the disjointness arrangeable at every scale with one fixed rectangle radius `r_0`: the
covering index `i(l)` grows with the band, but only *differences* of covering indices matter,
and those are bounded by `Delta_max` (6.14). The cost of the device is the duty cycle: each
pulse is on only a `S_*^{-1/2}` fraction of its window, so amplitudes must be `S_*^{1/4}`
larger — a polynomial price paid in a currency (powers of `S_* = l^2`) that is always beaten
by a gained power of `eps = 2^{-hl}`.

**Why Section 8 exists at all.** The pulses fix the *angular-mean radial flux*, but the
angular mean of the full residual has other pieces: a pressure that must be reconstructed, a
part that depends on the auxiliary torus even after angular averaging, and three global
*integral* obstructions. The last are the interesting ones. Radial integration of the mean
momentum equations, with the two moment constraints (8.2) imposed, kills every term except
three numbers `(P, J_theta, J_z)` (8.12), (8.15), (8.16). A compactly supported correction
cannot change these numbers unless it is designed to: they are the flux through the annulus.
Lemma 8.7 is the finite-dimensional fix — five bumps on a reserved "mean patch" where the
background is a clean power law `V_q = a x^{-1-2lam}`, `G_q = 0`, chosen so the resulting 5x5
matrix splits into two Vandermonde blocks and is therefore invertible. Everything else in
Section 8 is inverse operators (a compactly supported radial primitive `T_e` in Lemma 8.2, a
Fourier inverse `N^{-1}` of the fast time derivative in Lemma 8.6) built so that they do not
enlarge supports and do not spoil the class bounds.

---

## 5. Five concrete verification tasks for a newcomer

**(V1) Numerical test of the linearized pulse amplification in a frozen background shear.**
*Two levels.*
(a) *ODE level.* Integrate `z' = (diag(lam(v), -lam(v)) + E - m^2 d(v) I_2) z` on `[0, L_s]`
with `lam(v) = lam_0/sqrt(1 + s(v)^2)`, `d(v) = eps k^2 B_s^2 (1 + s(v)^2)`,
`s(v) = sigma(u_*/2 + u_* v/L_s)`, `B_s^2 = lam_0/(eps k^2 (1 + u_*^2)^{3/2})`, `eps k^2 = 1`,
`z(0) = (P(0), 0)`. Check: (i) `x(v) := z_+ + z_-` satisfies `c P(v) <= x <= C P(v)` with `P`
from (7.12) and the Gaussian bounds (7.16) — extract the constants `c, C` as functions of
`u_*`; (ii) `y/x -> c_0 sqrt(1 + s(v)^2) + O(S_*^{-1})`, i.e. (7.21); (iii) the amplification
ratio `P(L_s/2)/P(0)` really is `e^{Theta(L_s)}`; (iv) **the sensitivity to `E`** — set
`|E| = c_E/S_*` and confirm the envelope changes by an `O(1)` multiplicative factor
`exp(O(c_E L_s/S_*))`, then set `|E| = c_E S_*^{-1/2}` and watch the two-sided bound (7.16)
fail. This last item is the sharpest check available on (L1).
(b) *PDE level.* Linearize incompressible Navier-Stokes about a frozen axisymmetric shear
`u_B = V(r) e_theta + G(r) e_z` and evolve a single Kelvin mode with the paper's `(p, p_z)`
and `k = eps^{-1/2}`. Verify the transient growth, the crossover at `|s| = u_*`, and the
Reynolds stress `<u_r u_tan>` direction against (7.28). A frozen-shear 2-D/3-D spectral run
at `Re = O(1)` for the mode is enough. **This is the test I would run first**: it checks the
physical mechanism independently of all the bookkeeping.

**(V2) Symbolic verification of the whole of §7.1.** With a CAS, from
`g_0 = F_0(-a, b_s)`, `N = g_0/|g_0|`, `K = N^perp`, `v_s := (a^2 + b_s^2)/a`, derive
`lam_0^2 = -2 F_0 N_theta(2 F_0 N_theta + |g_0|) = 2 a F_0^2(1 - 2/v_s)` and
`c_0^2 = (v_s - 2)/2` (I verified these by hand; confirm the `v_s` identification against
(4.11), (4.20) in the paper, which I could not read). Then verify: `a_net(u_*) = 0`; the
`d a_net/dv` formula on p. 78; the two-column solve
`h_+ y_+ = (1/2)(-T_N/A_c - T_K/u_*)` etc.; the equivalence
`|c_0 T_K/T_N| < u_*/sqrt(1+u_*^2)  <==>  y_+ , y_- > 0`; and that (7.1)'s bound `< 1` is
strictly weaker, so `u_*` must be chosen large. **Deliverable: a table of the achievable cone
half-angle as a function of `u_*`, compared with the margin `eta_c` of Theorem 4.6(iii).**

**(V3) Direct test of Lemma 6.1.** (a) Numerically compute, for `J_g = [[3,1],[1,5]]` and
`Delta = 1, ..., Delta_max`, the minimum over a candidate finite set of rational centers
`{c_nu}` of `dist(c_mu - J_g^Delta c_nu, Z^2)`, and hence the largest admissible `r_0`; check
that it does *not* degenerate as `Delta_max` grows only logarithmically (use (6.14) with a
realistic `h`, `l_0`). (b) Independently verify the "bounded degree" claim by explicitly
building, for `l_0 = 20` say, the mesh boxes of bands `l_0, ..., l_0 + 4` with mesh
`S_l^{-3} = l^{-6}` in each of `(R, Z, T)` and counting how many boxes of band `l'` meet a
given box of band `l`; confirm the count is bounded independent of `l`. (c) Then verify the
*physical* consequence directly: sample `(r,t)`, compute `Y = v_r r^{d_r} + v_t t mod Z^2`,
and check that at most one label's rectangle contains it among a set of labels with
overlapping slow supports. **This turns Lemma 6.1 from a Diophantine assertion into a finite
computation.**

**(V4) Quantify every hidden `l_0` condition.** The paper fixes "a sufficiently large lower
band index `l_0`" once (p. 64) and thereafter uses at least six distinct smallness
requirements. Extract and tabulate them as explicit inequalities in `(h, l)`:
(i) `S_*^2(eps + eps^2 + k^{-1}) <= 1` (p. 77), i.e. `l^4 2^{-h l/2} <= 1`;
(ii) `S_*^{-3} S_* + eps^2 S_* + k^{-1} S_* <= C/S_*` for (7.9)/(L1);
(iii) `sqrt(L_s) <= L_s/5` for the `psi = 1` core to cover the Gaussian;
(iv) `C S_*^{-1/2} < eta_c` for the cone margin in (7.28)-(7.29);
(v) `e^{-c L_s/25} q^{-N} Q^{-M} S_*^C <= 1` for each `(M, N, C)` (p. 87);
(vi) `h(alpha + p kap_s) > N + L` for the flatness order in Lemma 8.2 (p. 92).
**Deliverable: `l_0(h, N)` in closed form, and a check that (v)-(vi) do not force `l_0` to
depend on quantities that are only fixed later** (this is the ordering-of-quantifiers risk in
a construction with this many parameters, and it is the failure mode I would look for
hardest).

**(V5) Verify the moment system (8.25) numerically and probe its conditioning.** Take
`G_q = 0`, `V_q = a x^{-1-2lam}`, a bump `eta_0` supported in `[x_0 - h_0, x_0 + h_0]` inside
the patch, `eta_j(x) = a_j^{-1} eta_0(x/a_j)`, `a_j = e^{jd}`, `j = 0, 1, 2`. Build `A_theta`
(3x3) and `A_z` (2x2) as displayed on p. 98, and:
(a) confirm the five power exponents `2, -2-2lam, -2lam` and `1, 1-2lam` by direct
integration; (b) confirm `int x^p eta_j dx = mu_p e^{jdp}`; (c) compute `det A_theta`,
`det A_z` and `||A^{-1}||` as functions of `(lam, d)`; (d) verify the decoupling and the
*crossing* of targets (`A_theta` takes `(0, -P, -J_z)`, `A_z` takes `(0, -J_theta)`);
(e) map the conditioning as `lam -> 0` and as `d -> 0`, and check the three bump supports stay
inside the reserved patch `I_m`; (f) solve (8.25) for random targets and verify the two
constraints `int R^2 Delta v dR = int R gamma_d dR = 0` hold to machine precision.
**Then, as a bonus, check the recomputation identity (8.26)-(8.27) symbolically** by
substituting `(beta + Delta beta, v + Delta v, gamma + Delta gamma)` into (8.3) and
confirming that the linear-in-base terms cancel exactly against the last three rows of
(8.25) — this is Lemma 8.8's whole content and is a pure algebra check.

*(Sixth task, if the reader has time: settle (U11). Attempt to prove
`|D^I a_sigma| <= C_I S_*^{b_I} sqrt(zeta) delta^{-a_I}` from (7.25) exactly as printed. I
believe it cannot be done, and that the intended hypothesis carries a factor `zeta`.)*

---

## 6. Overall assessment, with explicit uncertainty

**What I checked and believe.** I verified by hand, exactly, the following: the eigenstructure
and determinant of `J_g` and its eigenvectors `v_r, v_t` (6.2); the Diophantine bounds (6.7)
via algebraic conjugates; the identity `M_i = eps^{-kap_s} S_*^{-rho_g}` (6.6), which pins
down `d_r`; the covering-index bound (6.14) including the `8/l_0` term; the fiber count
`14^Delta` and the character computation (6.19); the label counts in Lemma 6.3; the three
cylindrical transport identities in Proposition 8.1; the vanishing of `a_net` at `y = u_*` and
its `v`-derivative formula (p. 78); the duty-cycle normalization `c_i L_s = 2 r_0` and
`h_sigma ~ S_*^{-1/2}` (7.27); the two-column positive solve and the reference cone
(7.28)-(7.29); the covariance identities (7.26), (7.34), and `C(W_0 + L Sigma)` on p. 85; the
curl identity `i k m n_Phi x C_m = t_m` and the `div curl = 0` cancellation in Lemma 7.7; the
scale identity `Q^{-2A} eps T_{0,*} = q^{-A-1/2} T_0` (7.30) and the amplitude/wavelength
cancellation on p. 81; the exponential-tail estimate on p. 87; the derivative-counting in
(8.10) and (8.19) against a 2-D lattice sum; both integration-by-parts identities and the
pressure moment in Proposition 8.4; the divergence-freeness of the azimuthal-potential
increment (8.14); all five power exponents, the moment formula `mu_p e^{jdp}`, the Vandermonde
structure, and the decoupling under `G_q = 0` in Lemma 8.7; and the term-by-term class table
in Lemma 8.8. **Every one of these came out right.** For a 166-page paper this is a
non-trivial signal: the internal algebra of Sections 6-8 is not the algebra of a document
that was assembled carelessly.

**What I found wrong.** One thing: **(7.25) as printed is missing a factor `zeta` in the
derivative bound**, and the statement it is used to prove (p. 83) does not follow without it.
The paper's own proof of Proposition 7.6 confirms the intended form. I rate this a typo with
high confidence (~90%).

**What I could not verify, and why it matters.**
- **Everything from Sections 3-5** that Sections 6-8 consume: Proposition 5.5 (the
  background, and in particular the flat residual, `b = O(eps)`, `G_q = 0` on the mean
  patch, and the one-sided endpoint derivatives), Theorem 4.6(iii) (the *strict* stress cone
  with a *uniform margin*), Theorem 4.6(vi) (the reserved mean patch), (4.27), (4.30). The
  cone margin is used in an essential quantitative way (it must exceed `C S_*^{-1/2}`), and
  `G_q = 0` is essential to the *invertibility* of (8.25). If either of these is weaker than
  advertised, Sections 6-8 do not close. **This is my largest single uncertainty about my own
  range.**
- **Everything in Section 9**, and in particular Proposition 9.6. **No net gain of any kind
  is demonstrated in pages 62-100.** The residual-improvement claim on p. 100 ("each cycle
  gains a fixed positive power of `eps`") is exactly the assertion that the whole scheme
  converges, and it is stated, not proved, in my range. Corollary 7.8's (7.42) and
  Lemma 8.8's exponent table are consistent with such a gain *provided* the iteration enters
  with `alpha > 1` and `rho > 1/2`; I could not check that it can.
- **The infinite-derivative-order and infinite-flatness diagonal arguments** (Lemma 9.7,
  Proposition 9.9). Sections 6-8 prove, for each fixed multi-index `I` and stage `j`, a bound
  with constants `C_{j,I}` and polynomial degrees `b_{j,I}` that are allowed to grow without
  any displayed rate. That is the correct architecture, but it defers the actual `C^infinity`
  statement entirely.
- The frame constants `a, b_s, v_s` from (4.11), (4.20), hence the numerical value of
  `lam_0`, `c_0`, and the constants in (7.16).

**My judgement.** Restricted to what Sections 6-8 claim to do, I did not find a structural
obstruction. The two devices I was asked to attack hardest both hold up: the auxiliary-torus
separation is arranged on the absolute torus and therefore *cannot* be undone by the physical
evaluation `Y(r,t)` — evaluation is a point map, and disjoint auxiliary supports give
literally disjoint physical supports (the pulses time-share, and the `S_*^{-1/2}` duty-cycle
price is correctly accounted for in (7.27)-(7.29)); and the amplification-then-damping
envelope is engineered exactly, with `B_s` defined so that the growth/damping crossover sits
at the pulse midpoint, which is what makes the temporal cutoff error `e^{-c l^2}` and hence
flat. The wave amplitude scale `q^{-1/2-h/2}` is genuinely `eps^{1/2}` smaller than the
background `q^{-1/2-h}`, with the quadratic covariance landing exactly on the target stress —
verified by an exponent computation that cancels identically. The "exponentially small tails"
claim is correct and is one of the cleanest arguments in my range; it works because
`S_* = l^2` grows superlinearly.

Against that: the paper's uniformity claims are **assertions supported by outlines**, not
proofs, in at least eight places (U1-U3, U8-U10, U12-U14, U16). The two that would worry me
most as a referee are **(U9)** — the induction on derivative order in (7.20), which is the
sole support for an infinite family of estimates on which the `C^infinity` conclusion rests —
and **(U17)** — that no gain is exhibited anywhere in pages 62-100, so the reader has no way
to sanity-check the iteration's arithmetic from within this range.

**Confidence, stated plainly.**
- That Sections 6-8 are *internally* consistent and that their algebra is correct: **high**
  (~85%), modulo the (7.25) typo, based on ~15 independent hand verifications.
- That the auxiliary-torus separation device works as claimed and is not circumvented by
  physical evaluation: **high** (~90%).
- That the amplification-then-damping argument controls values and low derivatives:
  **high** (~85%). That it controls *all* derivative orders as asserted: **moderate**
  (~60%) — the structure is standard but the induction is not written out.
- That Sections 6-8 as written are sufficient input for the theorem, i.e. that no hypothesis
  they need from Sections 4-5 is stronger than what those sections deliver: **cannot
  assess** — this is precisely what I could not read.
- That the paper as a whole proves Theorem 1.1: **cannot assess from this range, and I would
  caution against inferring anything from the quality of Sections 6-8.** A construction of
  this type can be locally impeccable for a hundred pages and fail at the single point where
  the gain is summed. That point is Proposition 9.6, and it is outside pages 62-100. Given
  the extraordinary nature of the claim, and the author line, the burden of the verification
  is on Sections 4-5 and 9, and a reader should go there next.

**Also worth noting as a matter of record**, since I was asked to mark what I could not
verify: the PDF's author line reads simply "OpenAI" with no individual authors, no
institutional affiliation, no arXiv identifier visible in my range, and a compile date of
2026-09-08 (today). I have no way from within this file to confirm provenance, and nothing in
my reading should be taken as bearing on it.
