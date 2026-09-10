"""Arithmetic checks on Sections 9-10 of "Finite time blowup for Navier-Stokes"
(OpenAI, compiled 2026-09-08), pages 100-126.

Everything here is elementary bookkeeping of exponents printed in the paper.
No claim of the paper is proved or disproved by this script; it only checks that
the displayed numbers are mutually consistent and computes two quantities the
paper does not display (the stage threshold of Lemma 9.8 and the Type II
exponents implied by A = 1/2 + h).

Run:  ./venv/bin/python commons/ns/artifacts/ns2026_cl_exponent_checks.py
Output file: ns2026_cl_exponent_checks_output.txt
"""

from fractions import Fraction as F
import sys
import numpy as np

out = []
def p(s=""):
    out.append(str(s))

p("environment: python %s, numpy %s" % (sys.version.split()[0], np.__version__))
p("source: OpenAI 2026, Finite time blowup for Navier-Stokes, pp. 100-126")
p("")

# ---------------------------------------------------------------- 1. schedule
p("1. Correction-cycle schedule (9.8): sigma_j = 1/5 + j/10, B_j = 1/2 + sigma_j,")
p("   C*_j = 1 + sigma_j.  Proposition 9.5 asserts B_0 = 0.7 and C*_0 = 1.2;")
p("   Proposition 9.6 asserts B_{j+1} = B_j + 1/10 and C*_{j+1} = C*_j + 1/10.")
sig = [F(1, 5) + F(j, 10) for j in range(11)]
B = [F(1, 2) + s for s in sig]
C = [F(1) + s for s in sig]
p("   j   sigma_j   B_j     C*_j")
for j in range(11):
    p("   %-3d %-9s %-7s %s" % (j, float(sig[j]), float(B[j]), float(C[j])))
p("   B_0 == 0.7 : %s      C*_0 == 1.2 : %s" % (B[0] == F(7, 10), C[0] == F(6, 5)))
p("   increments B_{j+1}-B_j all == 1/10 : %s" %
  all(B[j + 1] - B[j] == F(1, 10) for j in range(10)))
p("   increments C*_{j+1}-C*_j all == 1/10 : %s" %
  all(C[j + 1] - C[j] == F(1, 10) for j in range(10)))
p("   C*_j - B_j == 1/2 for all j : %s" % all(C[j] - B[j] == F(1, 2) for j in range(11)))
p("")

# ------------------------------------------------- 2. closure inequalities p.111
k = F(1, 100000)          # kappa_s = 1e-5
Bv = F(7, 10)             # the proof uses B >= 0.7
Cv = Bv + F(1, 2)         # C_* = B + 1/2
H = Cv - 2 * k            # (9.14): H = C_* - 2 kappa_s
p("2. The five closure inequalities displayed on p. 111, at kappa_s = 1e-5, B = 0.7.")
i1_lhs = min(F(1, 2) - 3 * k, F(1, 2) - k, Bv - k, F(2, 5))
i2_lhs = min(F(1, 2) - 4 * k, F(2, 5) - k, F(1, 2) - 2 * k, Bv - 3 * k)
rows = [
    ("min{1/2-3k, 1/2-k, B-k, 0.4} >= 0.4", i1_lhs, F(2, 5)),
    ("min{1/2-4k, 0.4-k, 1/2-2k, B-3k} >= 0.4-k", i2_lhs, F(2, 5) - k),
    ("H - B = 1/2-2k > 0.1", H - Bv, F(1, 10)),
    ("min{0.17, 1-4k} > 0.1", min(F(17, 100), 1 - 4 * k), F(1, 10)),
    ("0.9 - 4k > 0.1", F(9, 10) - 4 * k, F(1, 10)),
]
for name, lhs, rhs in rows:
    p("   %-45s lhs=%-12s rhs=%-10s margin=%-10s holds=%s" %
      (name, float(lhs), float(rhs), float(lhs - rhs), lhs >= rhs))
p("   binding (margin zero): %s" %
  ", ".join(n for n, l, r in rows if l == r))
p("   H - B = 1/2 - 2k exactly : %s" % (H - Bv == F(1, 2) - 2 * k))
p("   operative per-cycle demand is 0.1; smallest achieved margin over the last")
p("   three lines is %s, so the induction is not tight." %
  float(min(l - r for n, l, r in rows[2:])))
p("")

# ---------------------------------------------- 3. the p. 108 Step 1 table row
p("3. p. 108, Step 1 table of new nonzero-harmonic exponents, evaluated at B = 0.7:")
p("     linear error                        B + 1/2 - 3 kappa_s = %s" % float(Bv + F(1, 2) - 3 * k))
p("     cross with old exact waves          B + 1/2 - kappa_s   = %s" % float(Bv + F(1, 2) - k))
p("     particular-correction self-interact 2B - kappa_s        = %s" % float(2 * Bv - k))
p("     cross with total wave correction    B + 0.4             = %s" % float(Bv + F(2, 5)))
mu = F(9, 10)             # (9.9): v, gamma, p_m in M_{0.9}
p("   Lemma 9.2(i)  wave x mean gain  = mu - 1/2 with mu = 0.9 : %s" % float(mu - F(1, 2)))
p("   Lemma 9.2(ii) wave x wave gain  = 1/2 - kappa_s          : %s" % float(F(1, 2) - k))
p("   printed row value 0.4 equals the 9.2(i) value : %s ; equals the 9.2(ii) value : %s" %
  (F(2, 5) == mu - F(1, 2), F(2, 5) == F(1, 2) - k))
p("   p. 109, Step 2, prints the same gain with the label 'mean interaction',")
p("   exponent B + 0.4 - kappa_s = %s" % float(Bv + F(2, 5) - k))
p("   consistency of (9.9) with the schedule: B_0 - kappa_s = %s > 0.68 : %s" %
  (float(Bv - k), (Bv - k) > F(68, 100)))
p("")

# --------------------------------------- 4. Lemma 9.8 stage threshold (not in paper)
p("4. Lemma 9.8: gain g_j = h j / 10, loss l_m = 2A + (m+1)(1 + 3h/2), A = 1/2 + h.")
p("   Stage j carries a positive net power of q at Cartesian derivative order m only")
p("   when g_j > l_m, i.e. j > (10/h) [2A + (m+1)(1 + 3h/2)].  The paper never")
p("   displays this threshold; every earlier stage is controlled by the cutoff alone.")
p("   h        m    l_m           minimal j with h j/10 > l_m")
for h in (0.01, 0.005, 0.001):
    A = 0.5 + h
    for m in (0, 1, 2, 4, 8):
        lm = 2 * A + (m + 1) * (1 + 1.5 * h)
        jmin = int(np.floor(10 * lm / h)) + 1
        p("   %-8s %-4d %-13.6f %d" % (h, m, lm, jmin))
p("   The threshold is linear in m and of order 1/h in size.")
p("")

# ------------------------------------------------ 5. Type II exponents / energy
p("5. Growth exponents.  A = 1/2 + h and D = 1/2 - h are definitions (notation table,")
p("   p. 20); 0 < h < 1/100 (Theorem 3.1).  Along (10.20)-(10.21) u_theta ~ tau^{-A}.")
p("   h        tau^{1/2}|u| exponent   core energy   dissipation   L^3 norm   2h/D")
for h in (0.001, 0.01, 1.0 / 6.0 - 1e-9, 0.2):
    A = 0.5 + h
    t1 = 0.5 - A                      # tau^{1/2} ||u||_inf  ~ tau^{-h}
    ecore = 1.5 - h - 2 * A           # volume tau^{3/2-h} times tau^{-2A}
    diss = 1.5 - h - 2 * A - 1        # volume times squared radial derivative
    l3 = (1.5 - h - 3 * A) / 3.0      # (volume * |u|^3)^{1/3}
    p("   %-8.5f tau^%-19.5f tau^%-11.5f tau^%-11.5f tau^%-9.5f %.5f" %
      (h, t1, ecore, diss, l3, 2 * h / (0.5 - h)))
p("   core energy exponent 1/2-3h > 0 and the dissipation integral converges")
p("   exactly when h < 1/6; the exterior-tail condition 2h/D < 1 is the same")
p("   condition (2h < 1/2 - h).  Both hold with margin at h < 1/100.")
p("   L^3 exponent -4h/3 checks: (3/2 - h - 3A)/3 = -4h/3 for every h : %s" %
  bool(np.allclose([(1.5 - h - 3 * (0.5 + h)) / 3.0 + 4 * h / 3.0
                    for h in (0.001, 0.01, 0.05)], 0.0)))
p("   tau^{1/2} ||u||_inf ~ tau^{-h} is unbounded, so the constructed blowup is")
p("   Type II in the pointwise sense, for every h > 0.")
p("")

# ------------------------------------------------------ 6. viscosity rescaling
p("6. Viscosity rescaling (10.22)-(10.23): u_nu(x,t) = nu^{1/2} u(x/nu^{1/2}, t).")
p("   Energy scaling ||u_nu||_2^2 = nu^{5/2} ||u||_2^2 comes from nu^{1/2 * 2} times")
p("   the Jacobian nu^{3/2}: exponent 1 + 3/2 = %s (paper: 5/2)." % (1 + F(3, 2)))
p("   Dissipation: nu * |grad u_nu|^2 * dx scales as nu^{1} * nu^{1-1} * nu^{3/2}")
p("   = nu^{5/2}, matching (10.23).")
p("   Force derivatives: d_x^alpha d_t^m f_nu = nu^{(1-|alpha|)/2} (d_y^alpha d_t^m f).")
for a in range(4):
    p("     |alpha| = %d  ->  factor nu^{%s}" % (a, F(1 - a, 2)))
p("   So spatial derivative norms of the force grow like nu^{-(|alpha|-1)/2} as")
p("   nu -> 0, while the spatial support shrinks to K_nu = nu^{1/2} K and the time")
p("   support K x [0,2] is unchanged.")

text = "\n".join(out) + "\n"
with open(__file__.replace(".py", "_output.txt"), "w") as fh:
    fh.write(text)
print(text)
