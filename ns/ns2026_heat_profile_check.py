"""
Numerical audit of Equation (4.29) and Lemma 4.8(iii) of
"Finite time blowup for Navier-Stokes" (OpenAI, compiled 2026-09-08), pages 33 and 35.

The paper defines, for A = 1/2 + h,

    Hcal(Z) = (1/Gamma(1+h)) * int_0^infty e^{-v} v^h (1 + Z v)^{-h} dv,        (4.29)

and asserts:
  (a)  E = c_inf X^{-A} Hcal(2 d / X) for X >= X_b, with d = 1 - eta^2;
  (b)  the physical swirl there is the exact radial heat solution
         K(r,t) = c_inf s^{-A} Hcal(2 tau / s),   s = r^2/2,  tau = 1 - t,
       i.e. d_t K = (d_rr + r^{-1} d_r - r^{-2}) K;
  (c)  Hcal is positive and smooth for Z >= 0 including the one-sided endpoint Z = 0;
  (d)  K_r < 0 for r > 0.

By hand (see the accompanying record) the substitution K = s^{-A} Hcal(Z), Z = 2 tau / s
reduces (b) to the ODE

    Z^2 Hcal'' + ((2 + 2h) Z + 1) Hcal' + h(1+h) Hcal = 0.                     (*)

This script checks (*), Hcal(0) = 1, the full PDE (b) by centered finite differences in
(r,t), the sign claims (c),(d), and two consequences used elsewhere in Section 4:

  (e)  the radial shear coefficient a = 2 - 2 D_X log H, H = sqrt(2X) E,
       equals 4/5 on the reference inner branch E = P_* f(eta) x^{1/10} of Lemma 4.8(i)
       (the value a(X_i,eta) = 4/5 quoted on page 37, and the window .7 <= a <= .9 on p.38);
  (f)  on the heat exterior a = 2 + 2h + 2 Z Hcal'/Hcal, so a <= 2 + 2h always, and
       a > 2 + h iff -Z Hcal'/Hcal < h/2 -- the inequality (A.56) asserts in the
       stronger form -Z Hcal'/Hcal < h/4.  We report where that holds.

No claim in this script is a proof; it is a numerical cross-check of hand algebra.

Run:  ./venv/bin/python commons/ns/artifacts/ns2026_heat_profile_check.py
Output: ns2026_heat_profile_check_output.txt
"""

import sys
import numpy as np
import scipy
from scipy.integrate import quad
from scipy.special import gamma as Gamma

OUT = []


def emit(s=""):
    OUT.append(s)
    print(s)


def Hcal(Z, h, order=0):
    """order-th derivative of Hcal at Z, by differentiating under the integral sign.

    Hcal^{(k)}(Z) = (-1)^k h(h+1)...(h+k-1)/Gamma(1+h)
                    * int_0^inf e^{-v} v^{h+k} (1+Zv)^{-h-k} dv
    """
    pref = 1.0
    for j in range(order):
        pref *= -(h + j)
    pref /= Gamma(1.0 + h)
    k = order

    def integrand(v):
        return np.exp(-v) * v ** (h + k) * (1.0 + Z * v) ** (-h - k)

    val, err = quad(integrand, 0.0, np.inf, limit=400)
    return pref * val, abs(pref) * err


def main():
    emit("environment: python %s" % sys.version.split()[0])
    emit("numpy %s   scipy %s" % (np.__version__, scipy.__version__))
    emit("")

    # ------------------------------------------------------------------
    # 1. Hcal(0) = 1 and positivity / monotonicity
    # ------------------------------------------------------------------
    emit("=== 1. Hcal(0) = 1, Hcal > 0, Hcal' < 0 ===")
    for h in (1e-3, 1e-2, 0.05, 0.2, 0.45):
        H0, _ = Hcal(0.0, h, 0)
        emit("  h=%.4g   Hcal(0) = %.15f   (target 1)" % (h, H0))
    emit("")
    h_show = 0.01
    for Z in (0.0, 0.01, 0.1, 1.0, 10.0, 100.0):
        H, _ = Hcal(Z, h_show, 0)
        Hp, _ = Hcal(Z, h_show, 1)
        emit("  h=%.3g Z=%-8.4g Hcal=%.12f  Hcal'=%.6e" % (h_show, Z, H, Hp))
    emit("")

    # ------------------------------------------------------------------
    # 2. the reduced ODE (*)
    # ------------------------------------------------------------------
    emit("=== 2. ODE  Z^2 H'' + ((2+2h)Z + 1) H' + h(1+h) H = 0 ===")
    emit("    (residual normalised by max(1,|each term|) scale)")
    worst = 0.0
    for h in (1e-3, 1e-2, 0.05, 0.2, 0.45):
        for Z in (0.0, 0.05, 0.5, 2.0, 20.0, 200.0):
            H, eH = Hcal(Z, h, 0)
            H1, e1 = Hcal(Z, h, 1)
            H2, e2 = Hcal(Z, h, 2)
            res = Z * Z * H2 + ((2.0 + 2.0 * h) * Z + 1.0) * H1 + h * (1.0 + h) * H
            scale = max(abs(Z * Z * H2), abs(((2 + 2 * h) * Z + 1) * H1),
                        abs(h * (1 + h) * H), 1e-300)
            rel = abs(res) / scale
            worst = max(worst, rel)
            emit("  h=%-6.4g Z=%-7.4g  residual=%+.3e  relative=%.3e" % (h, Z, res, rel))
    emit("  worst relative ODE residual over the grid: %.3e" % worst)
    emit("")

    # ------------------------------------------------------------------
    # 3. the physical swirl heat equation, by finite differences
    # ------------------------------------------------------------------
    emit("=== 3. d_t K = (d_rr + r^{-1} d_r - r^{-2}) K  for K = s^{-A} Hcal(2 tau / s) ===")
    emit("    centered differences, step chosen by scale; c_inf = 1")

    def K(r, t, h):
        A = 0.5 + h
        s = 0.5 * r * r
        tau = 1.0 - t
        Z = 2.0 * tau / s
        val, _ = Hcal(Z, h, 0)
        return s ** (-A) * val

    worst3 = 0.0
    for h in (1e-3, 1e-2, 0.05, 0.2):
        for (r, t) in ((0.7, 0.3), (1.0, 0.5), (2.0, 0.9), (0.35, 0.99), (5.0, 0.0)):
            dr = 1e-4 * r
            dt = 1e-5
            Kc = K(r, t, h)
            Krr = (K(r + dr, t, h) - 2 * Kc + K(r - dr, t, h)) / dr ** 2
            Kr = (K(r + dr, t, h) - K(r - dr, t, h)) / (2 * dr)
            Kt = (K(r, t + dt, h) - K(r, t - dt, h)) / (2 * dt)
            lhs = Kt
            rhs = Krr + Kr / r - Kc / r ** 2
            rel = abs(lhs - rhs) / max(abs(lhs), abs(rhs), 1e-300)
            worst3 = max(worst3, rel)
            emit("  h=%-6.4g r=%-5.3g t=%-5.3g  d_tK=%+.8e  (Lap-1/r^2)K=%+.8e  rel=%.2e"
                 % (h, r, t, lhs, rhs, rel))
            if abs(Kr) > 0:
                pass
    emit("  worst relative PDE residual: %.3e   (finite-difference truncation dominated)" % worst3)
    emit("")

    # ------------------------------------------------------------------
    # 4. K_r < 0
    # ------------------------------------------------------------------
    emit("=== 4. K_r < 0 for r > 0 ===")
    emit("    hand identity: K_s = -s^{-A-1} (A Hcal + Z Hcal'),")
    emit("    A Hcal + Z Hcal' = (1/Gamma(1+h)) int e^{-v} v^h (1+Zv)^{-h-1} [A + (A-h) Z v] dv > 0")
    allneg = True
    for h in (1e-3, 1e-2, 0.05, 0.2):
        A = 0.5 + h
        for (r, t) in ((0.5, 0.5), (1.0, 0.5), (3.0, 0.2), (0.2, 0.95)):
            dr = 1e-4 * r
            Kr = (K(r + dr, t, h) - K(r - dr, t, h)) / (2 * dr)
            s = 0.5 * r * r
            Z = 2.0 * (1.0 - t) / s
            H, _ = Hcal(Z, h, 0)
            H1, _ = Hcal(Z, h, 1)
            comb = A * H + Z * H1
            allneg = allneg and (Kr < 0)
            emit("  h=%-6.4g r=%-5.3g t=%-5.3g  K_r=%+.6e   A*H+Z*H'=%+.8f" % (h, r, t, Kr, comb))
    emit("  all sampled K_r negative: %s" % allneg)
    emit("")

    # ------------------------------------------------------------------
    # 5. shear coefficient a on the two reference branches
    # ------------------------------------------------------------------
    emit("=== 5. a = 2 - 2 D_X log H, H = sqrt(2X) E ===")
    emit("  (a) reference inner branch of Lemma 4.8(i): E = P_* f(eta) x^{1/10}")
    emit("      H propto X^{1/2 + 1/10} = X^{3/5} so l = D_X log H = 3/5 exactly and a = 2 - 6/5 = 4/5")
    # numeric confirmation with an arbitrary P_*, eta, X_R
    Pstar, eta, XR = 7.0, 0.31, 123.0
    f = 1.0 / (1.0 + eta ** 2)

    def a_inner(X):
        E = Pstar * f * (X / XR) ** 0.1
        return np.sqrt(2.0 * X) * E

    for X in (1.0, 10.0, 123.0, 1000.0):
        dX = 1e-5 * X
        Hh = a_inner(X)
        l = X * (a_inner(X + dX) - a_inner(X - dX)) / (2 * dX) / Hh
        emit("      X=%-8.4g  l=%.12f  a=2-2l=%.12f" % (X, l, 2 - 2 * l))
    emit("")
    emit("  (b) heat exterior: E = c_inf X^{-A} Hcal(2d/X); D_X Z = -Z so")
    emit("      l = -h - Z Hcal'/Hcal  and  a = 2 + 2h + 2 Z Hcal'/Hcal <= 2 + 2h.")
    emit("      Lemma 4.9 needs 2 + h < a, i.e. -Z Hcal'/Hcal < h/2;")
    emit("      (A.56) asserts the stronger -Z Hcal'/Hcal < h/4.")
    for h in (1e-3, 1e-2, 0.05):
        emit("      h = %g   (h/4 = %.6g)" % (h, h / 4))
        for Z in (1e-4, 1e-3, 1e-2, 0.1, 0.5, 1.0):
            H, _ = Hcal(Z, h, 0)
            H1, _ = Hcal(Z, h, 1)
            m = -Z * H1 / H
            a = 2 + 2 * h - 2 * m
            emit("        Z=%-8.4g  -Z H'/H = %.6e   a = %.10f   < h/4 ? %s"
                 % (Z, m, a, m < h / 4))
    emit("")
    emit("  Note: Z = 2d/X with d = 1 - eta^2 <= 1, so Z <= 2/X; the condition")
    emit("  -Z Hcal'/Hcal < h/4 is a largeness requirement on the outer radius,")
    emit("  consistent with Lemma 4.8's 'for every X_R >= R_*'.  For small h the")
    emit("  admissible Z window itself shrinks like h, which the table shows.")
    emit("")

    with open(__file__.replace(".py", "_output.txt"), "w") as fh:
        fh.write("\n".join(OUT) + "\n")


if __name__ == "__main__":
    main()
