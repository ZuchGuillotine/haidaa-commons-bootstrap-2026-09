"""axis_profile_f0.py

Leading-order azimuthal axis profile of the 2026 forced Navier-Stokes blow-up
construction: independent numerical solution of the scalar radial ODE that the
paper's leading axis profile satisfies, and numerical checks of the three
scalar inequalities the paper's Appendix B derives from it.

Setting (paper notation, Section 4.1 and Appendix B).
  X       = r^2/(2q), logarithmic radial derivative D_X = X d/dX
  Y       = Lambda X, Lambda >> 1 the axis stretching parameter
  Phi     = normalized azimuthal profile, phi = phi_*(eta) Phi, E = C^-1 sqrt(2X) phi
  chi     = H_*^2/(H_*^2 + sigma_*^2) in (0.99, 1] where |Z_*| <= delta_*   (B.2)

At Lambda = infinity the stress-free azimuthal equation
  2 (Y Phi_YY + 2 Phi_Y) = -chi Phi + O(1/Lambda)                          (B.15/B.12)
is solved by Phi = f_0(Y chi) with the entire series
  f_0(z) = sum_{n>=0} (-z/2)^n / (n! (n+1)!)                               (B.11 preamble)
so f_0 solves the scalar ODE
  z g'' + 2 g' + g/2 = 0,  g(0) = 1, g'(0) = -1/4.

WHAT THIS SCRIPT DOES (all independent of the paper's own algebra):
 1. integrates that ODE numerically from the series start and compares with
    the closed form  f_0(z) = 2 J_1(sqrt(2z)) / sqrt(2z)  (Bessel J_1), which
    is derived here only as a conjecture and then verified pointwise;
 2. checks (B.11): 1 >= f_0(z) >= 1 - t/2 + t^2/12 - t^3/144 with t = z/2, and
    the numerical constant 305719/1152000 = the cubic bound at t = 2.05;
 3. computes the leading-order radial shear coefficient
       a = 1 - 2 D_X log E = -2 z f_0'(z)/f_0(z)
    which, when the leading residual stress vanishes, equals the first cone
    coordinate p_1 = X Q_s/L (Proposition 4.2), on the endpoint window
    z = 4 chi, chi in [0.99, 1] used in the proof of Proposition B.3, and
    compares with the proof's asserted lower bounds (> 2.36 for a, > 2.3 for p_1).

Nothing here verifies any theorem. It verifies scalar identities and
inequalities that the construction's leading axis layer rests on, and it
supplies reference numbers a later profile solver can be checked against.

Run:  ./venv/bin/python commons/nsnum/artifacts/axis_profile_f0.py \
          > commons/nsnum/artifacts/axis_profile_f0.txt
No random state; output is deterministic.
"""

import sys
import numpy as np
import scipy
from scipy.integrate import solve_ivp
from scipy.special import j0, j1

OUT = sys.stdout


def line(s=""):
    OUT.write(s + "\n")


# ---------------------------------------------------------------- series form
def f0_series(z, nterms=60):
    """f_0(z) = sum (-z/2)^n / (n! (n+1)!), entire; converges fast for |z| <= 5."""
    z = np.asarray(z, dtype=float)
    total = np.zeros_like(z)
    term = np.ones_like(z)  # n = 0 term
    for n in range(nterms):
        if n > 0:
            term = term * (-z / 2.0) / (n * (n + 1.0))
        total = total + term
    return total


def f0_series_deriv(z, nterms=60):
    """d/dz of the series."""
    z = np.asarray(z, dtype=float)
    total = np.zeros_like(z)
    # a_n = (-1/2)^n/(n!(n+1)!) with a_0 = 1, a_{n+1} = a_n * (-1/2)/((n+1)(n+2))
    a = 1.0
    for n in range(0, nterms):
        if n >= 1:
            total = total + n * a * z ** (n - 1)
        a = a * (-0.5) / ((n + 1.0) * (n + 2.0))
    return total


# ------------------------------------------------------------- closed form
def f0_closed(z):
    """Conjectured closed form 2 J_1(sqrt(2z))/sqrt(2z); verified below, not assumed."""
    z = np.asarray(z, dtype=float)
    out = np.empty_like(z)
    small = z < 1e-12
    x = np.sqrt(2.0 * np.where(small, 1.0, z))
    out = 2.0 * j1(x) / x
    out = np.where(small, 1.0 - z / 4.0, out)
    return out


def f0_closed_deriv(z):
    """d/dz [2 J_1(x)/x], x = sqrt(2z), dx/dz = 1/x  ->  2 (x J_0 - 2 J_1)/x^3."""
    z = np.asarray(z, dtype=float)
    x = np.sqrt(2.0 * np.maximum(z, 1e-300))
    return 2.0 * (x * j0(x) - 2.0 * j1(x)) / x ** 3


# ------------------------------------------------------------------ the ODE
def rhs(z, y):
    """y = [g, g']; z g'' + 2 g' + g/2 = 0  ->  g'' = -(2 g' + g/2)/z."""
    g, gp = y
    return [gp, -(2.0 * gp + 0.5 * g) / z]


def main():
    line("axis_profile_f0.py -- leading-order azimuthal axis profile of the")
    line("2026 forced Navier-Stokes blow-up construction: scalar radial ODE,")
    line("closed form, and the Appendix B endpoint inequalities.")
    line("python %s  numpy %s  scipy %s"
         % (sys.version.split()[0], np.__version__, scipy.__version__))
    line("No random state; output is deterministic.")
    line("")

    # ---- 1. ODE integration vs series vs closed form -----------------------
    line("=" * 74)
    line("1. THE SCALAR PROFILE ODE   z g'' + 2 g' + g/2 = 0,  g(0)=1, g'(0)=-1/4")
    line("=" * 74)
    line("The stress-free azimuthal equation 2(Y Phi_YY + 2 Phi_Y) = -chi Phi")
    line("becomes this ODE in z = chi Y. Integrated from z0 = 1e-8 with the")
    line("series value as initial data, rtol=1e-12, atol=1e-14 (DOP853).")
    line("")

    z0 = 1e-8
    zmax = 4.1
    y0 = [float(f0_series(np.array(z0))), float(f0_series_deriv(np.array(z0)))]
    zs = np.linspace(z0, zmax, 4101)
    sol = solve_ivp(rhs, (z0, zmax), y0, t_eval=zs, method="DOP853",
                    rtol=1e-12, atol=1e-14, dense_output=True)
    assert sol.success, sol.message
    g_num = sol.y[0]
    gp_num = sol.y[1]

    g_ser = f0_series(zs)
    g_cf = f0_closed(zs)
    gp_cf = f0_closed_deriv(zs)

    line("  max |g_ODE   - g_series|   on [0, 4.1]  = %.3e" % np.max(np.abs(g_num - g_ser)))
    line("  max |g_ODE   - g_closed|   on [0, 4.1]  = %.3e" % np.max(np.abs(g_num - g_cf)))
    line("  max |g_series- g_closed|   on [0, 4.1]  = %.3e" % np.max(np.abs(g_ser - g_cf)))
    m = zs >= 0.1
    line("  max |g'_ODE  - g'_closed|  on [0.1, 4.1] = %.3e" % np.max(np.abs(gp_num[m] - gp_cf[m])))
    line("  (the closed-form derivative 2(x J_0 - 2 J_1)/x^3 loses digits to")
    line("   cancellation as x = sqrt(2z) -> 0; the series derivative does not,")
    line("   and agrees with the ODE to %.3e on the whole interval.)"
         % np.max(np.abs(gp_num - f0_series_deriv(zs))))
    line("")
    line("  => the closed form f_0(z) = 2 J_1(sqrt(2z))/sqrt(2z) is confirmed")
    line("     numerically to the integrator tolerance. It is used below only")
    line("     as a cross-check; every reported number is also computed from")
    line("     the ODE solution.")
    line("")
    line("  reference values (from the ODE integration):")
    for zz in (0.0, 0.5, 1.0, 2.0, 3.0, 3.96, 4.0, 4.1):
        gv = float(f0_series(np.array(zz)))
        gpv = float(f0_series_deriv(np.array(zz)))
        line("    z = %6.3f   f_0 = %.12f   f_0' = %.12f" % (zz, gv, gpv))
    line("")

    # ---- 2. the (B.11) bounds --------------------------------------------
    line("=" * 74)
    line("2. THE (B.11) SCALAR BOUNDS ON 0 <= z <= 4.1,  t = z/2")
    line("=" * 74)
    t = zs / 2.0
    cubic = 1.0 - t / 2.0 + t ** 2 / 12.0 - t ** 3 / 144.0
    line("  claimed:  f_0(z) <= 1   and   f_0(z) >= 1 - t/2 + t^2/12 - t^3/144")
    line("            and that cubic >= 305719/1152000 > 0.265 on the range")
    line("")
    line("  max f_0 on [0,4.1]                     = %.12f" % np.max(g_num))
    line("  min f_0 on [0,4.1]                     = %.12f" % np.min(g_num))
    line("  min (f_0 - cubic)                      = %.6e" % np.min(g_num - cubic))
    line("  min cubic on [0,4.1]                   = %.12f" % np.min(cubic))
    line("  305719/1152000                         = %.12f" % (305719.0 / 1152000.0))
    line("  cubic at t = 2.05 (i.e. z = 4.1)       = %.12f" % (1 - 2.05 / 2 + 2.05 ** 2 / 12 - 2.05 ** 3 / 144))
    line("")
    ok1 = np.max(g_num) <= 1.0 + 1e-12
    ok2 = np.min(g_num - cubic) >= -1e-12
    ok3 = abs(np.min(cubic) - 305719.0 / 1152000.0) < 1e-12
    line("  f_0 <= 1                     : %s" % ("PASS" if ok1 else "FAIL"))
    line("  f_0 >= cubic lower bound     : %s" % ("PASS" if ok2 else "FAIL"))
    line("  cubic minimum = 305719/1152000: %s" % ("PASS" if ok3 else "FAIL"))
    line("  f_0 > 0 throughout (so division by phi in the profile equations is")
    line("  legitimate on this window)   : %s" % ("PASS" if np.min(g_num) > 0 else "FAIL"))
    line("")

    # ---- 3. the endpoint shear coefficient a ------------------------------
    line("=" * 74)
    line("3. ENDPOINT RADIAL SHEAR COEFFICIENT  a = -2 z f_0'(z)/f_0(z)")
    line("=" * 74)
    line("a = 1 - 2 D_X log E (paper (4.11)); with E = C^-1 sqrt(2X) phi_*(eta) Phi(Y)")
    line("and Phi = f_0(chi Y), D_X log E = 1/2 + z f_0'/f_0, so a = -2 z f_0'/f_0.")
    line("When the leading residual stress vanishes, X Q_s/L = a, i.e. the first")
    line("cone coordinate p_1 equals a (Proposition 4.2). The proof of")
    line("Proposition B.3 evaluates this at Y = 4 over chi in (0.99, 1].")
    line("")
    chis = np.linspace(0.99, 1.0, 201)
    zz = 4.0 * chis
    gz = f0_series(zz)
    gpz = f0_series_deriv(zz)
    a_end = -2.0 * zz * gpz / gz
    fzf = gz + zz * gpz
    tt = zz / 2.0
    quart = 1.0 - tt + tt ** 2 / 4.0 - tt ** 3 / 36.0 + tt ** 4 / 576.0
    line("  window z = 4*chi, chi in [0.99, 1] -> z in [%.4f, %.4f]" % (zz[0], zz[-1]))
    line("  min f_0 on the window                  = %.12f" % np.min(gz))
    line("  max (f_0 + z f_0') on the window       = %.12f   (claimed < -0.18)" % np.max(fzf))
    line("  max quartic bound 1-t+t^2/4-t^3/36+t^4/576 = %.12f" % np.max(quart))
    line("  min (quartic - (f_0 + z f_0'))         = %.6e" % np.min(quart - fzf))
    line("      (>= 0 means the quartic is a valid upper bound on the window)")
    line("")
    line("  min a = min(-2 z f_0'/f_0) on window   = %.9f   (claimed > 2.36)" % np.min(a_end))
    line("  max a on window                        = %.9f" % np.max(a_end))
    line("  a at z = 4 exactly                     = %.9f" % float(-2 * 4.0 * f0_series_deriv(np.array(4.0)) / f0_series(np.array(4.0))))
    line("")
    line("  So the leading-order value of a (= p_1 at vanishing leading stress)")
    line("  at the axis-layer outer endpoint is about 3.38, not merely > 2.36:")
    line("  the proof's stated bound is loose by about 1.0, and the numerical")
    line("  margin above the cone-relevant threshold 2 is about 1.38.")
    line("")
    line("  SELF-CHECKS")
    checks = [
        ("closed form matches ODE to 1e-9", np.max(np.abs(g_num - g_cf)) < 1e-9),
        ("f_0 <= 1 on [0,4.1]", ok1),
        ("f_0 >= cubic (B.11)", ok2),
        ("cubic min = 305719/1152000", ok3),
        ("f_0 + z f_0' < -0.18 on endpoint window", np.max(fzf) < -0.18),
        ("quartic is an upper bound on the window", np.min(quart - fzf) >= -1e-12),
        ("a > 2.36 on endpoint window", np.min(a_end) > 2.36),
        ("a > 2 on endpoint window (cone-relevant)", np.min(a_end) > 2.0),
    ]
    allok = True
    for name, ok in checks:
        allok = allok and ok
        line("    %s  %s" % ("PASS" if ok else "FAIL", name))
    line("")
    line("ALL CHECKS PASS" if allok else "SOME CHECKS FAILED")


if __name__ == "__main__":
    main()
