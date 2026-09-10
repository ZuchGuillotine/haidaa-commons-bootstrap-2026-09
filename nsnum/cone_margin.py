"""cone_margin.py

Independent numerical evaluator for the "admissible stress cone condition" of
the 2026 forced Navier-Stokes blow-up construction, and numerical checks of the
scalar inequalities the paper's Lemma 4.5 and Appendix A.5 rest on.

Definitions used (paper Section 4.3, equations (4.20), (4.21), (4.22), (4.35)).
Inputs are four real numbers at a point of the pulse annulus:
  a       > 0   radial azimuthal shear coefficient, a = 1 - 2 D_X log E
  b_s           radial axial shear coefficient,     b_s = 2 D_X U / E
  p_1, p_2      the two inviscid stress coordinates p_s = (X Q_s/L, X N_s/(L E))
Derived:
  t_s = -b_s/a,  v_s = a (1 + t_s^2) = a + b_s^2/a
  P_c = p_1 + t_s p_2 = p_1 - (b_s/a) p_2
  J_c = p_2 - t_s p_1 = p_2 + (b_s/a) p_1
Relaxed cone condition:      P_c > 2  and  v_s < U(P_c, J_c), where
  U(P_c, J_c) = P_c + J_c^2/4 - |J_c| sqrt((P_c - 2)/2 + J_c^2/16)
Admissible stress cone condition: the relaxed condition together with v_s > 2.
Lemma 4.5 (first assertion): for v_s > 2 the admissible condition is equivalent to
  P_c > v_s   and   (v_s - 2) J_c^2 < 2 (P_c - v_s)^2.                     (4.22)
Cone-gap vector (4.35): Psi = (a, v_s - 2, P_c - v_s, 2(P_c - v_s)^2 - (v_s-2) J_c^2);
admissibility is exactly positivity of all four components.

WHAT THIS SCRIPT DOES
 1. implements both forms and tests their equivalence on 400000 pseudorandom
    parameter points (seed 20260909), reporting disagreements;
 2. checks the algebraic identity used in the proof of Lemma 4.5's second
    assertion, which converts the cone test into a condition on (a, b_s, w)
    with w = p_2/p_1 alone in the large-p_1 limit:
      (v_s - 2)(w + b_s/a)^2 - 2(1 - b_s w/a)^2
        = (1 + b_s^2/a^2) ( (a-2) w^2 + 2 b_s w + b_s^2/a - 2 );
 3. from that identity defines the large-p_1 margin
      m(a, b_s, w) = 2 - [ 2 b_s w + b_s^2/a + (a-2) w^2 ]
    (positive iff the relaxed cone is reachable by taking p_1 large enough)
    and maps it over (b_s, w) at the leading-order axis endpoint value
    a = 3.389130712 computed in artifact axis_profile_f0.txt;
 4. checks the two scalar suprema the paper's Appendix A.5 uses on the axial
    pulse interval:  sup_{R>=0} (-2 R^2 + 2.41 R) < 0.74  and
    sup_{R>=0} (-3.5 R^2 + 4.82 R) < 1.68 < 2.

Nothing here verifies any theorem, and nothing here uses the paper's actual
profiles: the audit of a computed profile against this cone test is the open
task. This script is the checked evaluator that such an audit would call.

Run:  ./venv/bin/python commons/nsnum/artifacts/cone_margin.py \
          > commons/nsnum/artifacts/cone_margin.txt
Seed: 20260909 (numpy default_rng). Output is deterministic.
"""

import sys
import numpy as np
import scipy

OUT = sys.stdout


def line(s=""):
    OUT.write(s + "\n")


def cone_coords(a, b_s, p1, p2):
    t_s = -b_s / a
    v_s = a * (1.0 + t_s ** 2)
    P_c = p1 + t_s * p2
    J_c = p2 - t_s * p1
    return t_s, v_s, P_c, J_c


def U_fun(P_c, J_c):
    inner = (P_c - 2.0) / 2.0 + J_c ** 2 / 16.0
    inner = np.where(inner < 0.0, np.nan, inner)
    return P_c + J_c ** 2 / 4.0 - np.abs(J_c) * np.sqrt(inner)


def admissible_definition(a, b_s, p1, p2):
    """(4.21) form: P_c > 2, v_s < U(P_c,J_c), v_s > 2."""
    _, v_s, P_c, J_c = cone_coords(a, b_s, p1, p2)
    with np.errstate(invalid="ignore"):
        Uv = U_fun(P_c, J_c)
        return (P_c > 2.0) & (v_s < Uv) & (v_s > 2.0)


def admissible_lemma45(a, b_s, p1, p2):
    """(4.22) form: v_s > 2, P_c > v_s, (v_s-2) J_c^2 < 2 (P_c - v_s)^2."""
    _, v_s, P_c, J_c = cone_coords(a, b_s, p1, p2)
    return (v_s > 2.0) & (P_c > v_s) & ((v_s - 2.0) * J_c ** 2 < 2.0 * (P_c - v_s) ** 2)


def psi_vector(a, b_s, p1, p2):
    """Cone-gap vector (4.35). All four components positive iff admissible."""
    _, v_s, P_c, J_c = cone_coords(a, b_s, p1, p2)
    return np.stack([a * np.ones_like(v_s),
                     v_s - 2.0,
                     P_c - v_s,
                     2.0 * (P_c - v_s) ** 2 - (v_s - 2.0) * J_c ** 2])


def margin_large_p1(a, b_s, w):
    """m = 2 - [2 b_s w + b_s^2/a + (a-2) w^2]; > 0 iff cone reachable for large p_1."""
    return 2.0 - (2.0 * b_s * w + b_s ** 2 / a + (a - 2.0) * w ** 2)


def main():
    line("cone_margin.py -- admissible stress cone evaluator for the 2026 forced")
    line("Navier-Stokes blow-up construction, with checks of Lemma 4.5 and of the")
    line("Appendix A.5 scalar suprema.")
    line("python %s  numpy %s  scipy %s"
         % (sys.version.split()[0], np.__version__, scipy.__version__))
    line("Seed 20260909 (numpy default_rng). Output is deterministic.")
    line("")

    rng = np.random.default_rng(20260909)
    N = 400000

    # ---- 1. equivalence of the two forms ---------------------------------
    line("=" * 74)
    line("1. EQUIVALENCE OF (4.21) AND (4.22)  [Lemma 4.5, first assertion]")
    line("=" * 74)
    a = 10.0 ** rng.uniform(-1.0, 1.5, N)          # a > 0
    b_s = rng.uniform(-6.0, 6.0, N)
    p1 = rng.uniform(-8.0, 30.0, N)
    p2 = rng.uniform(-15.0, 15.0, N)

    d = admissible_definition(a, b_s, p1, p2)
    l = admissible_lemma45(a, b_s, p1, p2)
    mism = np.where(d != l)[0]
    line("  samples                                 = %d" % N)
    line("  admissible by (4.21) definition         = %d" % int(d.sum()))
    line("  admissible by (4.22) Lemma 4.5 form     = %d" % int(l.sum()))
    line("  disagreements                           = %d" % mism.size)
    if mism.size:
        # report how close to the boundary the disagreements sit
        _, v_s, P_c, J_c = cone_coords(a[mism], b_s[mism], p1[mism], p2[mism])
        gap = np.abs(2.0 * (P_c - v_s) ** 2 - (v_s - 2.0) * J_c ** 2)
        line("  max |gap| at a disagreement (boundary)  = %.3e" % gap.max())
    line("")
    line("  Psi (4.35) all-positive vs (4.22): disagreements = %d"
          % int(np.sum((psi_vector(a, b_s, p1, p2) > 0).all(axis=0) != l)))
    line("")

    # ---- 2. the algebraic identity ---------------------------------------
    line("=" * 74)
    line("2. THE IDENTITY BEHIND LEMMA 4.5's SECOND ASSERTION")
    line("=" * 74)
    line("  (v_s - 2)(w + b_s/a)^2 - 2(1 - b_s w/a)^2")
    line("      = (1 + b_s^2/a^2) ( (a-2) w^2 + 2 b_s w + b_s^2/a - 2 )")
    w = rng.uniform(-5.0, 5.0, N)
    v_s = a + b_s ** 2 / a
    lhs = (v_s - 2.0) * (w + b_s / a) ** 2 - 2.0 * (1.0 - b_s * w / a) ** 2
    rhs = (1.0 + b_s ** 2 / a ** 2) * ((a - 2.0) * w ** 2 + 2.0 * b_s * w + b_s ** 2 / a - 2.0)
    rel = np.abs(lhs - rhs) / np.maximum(np.abs(lhs), 1.0)
    line("")
    line("  max relative residual over %d samples  = %.3e" % (N, rel.max()))
    line("  identity holds to machine precision      : %s"
          % ("PASS" if rel.max() < 1e-11 else "FAIL"))
    line("")
    line("  Consequence used below: in the limit p_1 -> +infinity with")
    line("  p_2 = w p_1, admissibility of the relaxed cone is governed by the")
    line("  sign of  m(a,b_s,w) = 2 - [2 b_s w + b_s^2/a + (a-2) w^2]  alone.")
    line("")

    # ---- 3. margin map at the leading-order axis endpoint value of a -----
    line("=" * 74)
    line("3. LARGE-p_1 CONE MARGIN AT THE AXIS-LAYER ENDPOINT VALUE OF a")
    line("=" * 74)
    a_end = 3.389130712   # from artifact axis_profile_f0.txt, z = 4, chi = 1
    line("  a = %.9f  (leading-order value of a = -2 z f_0'/f_0 at z = 4;" % a_end)
    line("  see artifact axis_profile_f0.txt). At vanishing leading residual")
    line("  stress this same number is the first cone coordinate p_1 = X Q_s/L.")
    line("")
    line("  m(a,b_s,w) = 2 - [2 b_s w + b_s^2/a + (a-2) w^2] on a grid:")
    line("")
    bs_grid = np.array([-1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5])
    w_grid = np.array([-1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5])
    hdr = "     b_s \\ w  " + "".join("%9.2f" % ww for ww in w_grid)
    line(hdr)
    for bb in bs_grid:
        row = "     %8.2f " % bb
        for ww in w_grid:
            row += "%9.4f" % margin_large_p1(a_end, bb, ww)
        line(row)
    line("")
    line("  (positive entries: the relaxed cone is reachable at that (b_s, w)")
    line("   by taking the stress magnitude p_1 large enough; negative entries:")
    line("   it is not reachable at any p_1.)")
    line("")
    # admissible w-interval at b_s = 0 and the optimal w
    aa = a_end
    line("  at b_s = 0 the admissible interval is |w| < sqrt(2/(a-2)) = %.6f"
          % np.sqrt(2.0 / (aa - 2.0)))
    line("  best w for given b_s is w* = -b_s/(a-2), giving")
    line("      m_max = 2 + 2 b_s^2/(a(a-2)),  e.g. b_s = 1 -> m_max = %.6f"
          % (2.0 + 2.0 / (aa * (aa - 2.0))))
    line("  so for every a > 2 there is a stress direction with margin >= 2;")
    line("  the cone condition constrains the DIRECTION w = p_2/p_1, not the")
    line("  existence of some admissible direction. This is why the audit that")
    line("  matters is the one that uses the profiles' own w, not a sweep.")
    line("")
    # sensitivity of the admissible w-interval to a
    line("  sensitivity of the admissible w-interval half-width to a (b_s = 0):")
    for av in (2.2, 2.6, 3.0, 3.389130712, 4.0, 6.0):
        line("      a = %11.6f   |w|_max = %.6f" % (av, np.sqrt(2.0 / (av - 2.0))))
    line("")

    # ---- 4. the Appendix A.5 scalar suprema ------------------------------
    line("=" * 74)
    line("4. THE APPENDIX A.5 SCALAR SUPREMA ON THE AXIAL PULSE INTERVAL")
    line("=" * 74)
    line("  The paper bounds b_s w <= sup_{R>=0}(-2 R^2 + 2.41 R) and states")
    line("  the value is < 0.74; and bounds the full sufficient-condition")
    line("  quantity 2 b_s w + b_s^2/a + (a-2) w^2 by sup_{R>=0}(-3.5 R^2 + 4.82 R)")
    line("  and states the value is < 1.68, hence < 2 with a margin.")
    line("")
    s1 = 2.41 ** 2 / 8.0
    s2 = 4.82 ** 2 / 14.0
    R = np.linspace(0.0, 3.0, 3000001)
    s1n = np.max(-2.0 * R ** 2 + 2.41 * R)
    s2n = np.max(-3.5 * R ** 2 + 4.82 * R)
    line("  sup(-2 R^2 + 2.41 R)  exact 2.41^2/8   = %.9f" % s1)
    line("                        grid maximum      = %.9f" % s1n)
    line("                        claimed bound     < 0.74     %s"
          % ("PASS" if s1 < 0.74 else "FAIL"))
    line("  sup(-3.5 R^2 + 4.82 R) exact 4.82^2/14 = %.9f" % s2)
    line("                        grid maximum      = %.9f" % s2n)
    line("                        claimed bound     < 1.68     %s"
          % ("PASS" if s2 < 1.68 else "FAIL"))
    line("  margin below the cone threshold 2      = %.9f" % (2.0 - s2))
    line("  i.e. on that interval the sufficient cone test is satisfied with")
    line("  about 17 percent of headroom in the o(1) terms that were dropped.")
    line("")

    line("  SELF-CHECKS")
    checks = [
        ("(4.21) and (4.22) agree on all samples", mism.size == 0),
        ("Psi (4.35) positivity agrees with (4.22)",
         int(np.sum((psi_vector(a, b_s, p1, p2) > 0).all(axis=0) != l)) == 0),
        ("Lemma 4.5 identity to 1e-11", rel.max() < 1e-11),
        ("sup(-2R^2+2.41R) < 0.74", s1 < 0.74),
        ("sup(-3.5R^2+4.82R) < 1.68", s2 < 1.68),
        ("grid maxima match closed forms to 1e-9",
         abs(s1 - s1n) < 1e-9 and abs(s2 - s2n) < 1e-9),
        ("margin at a=3.3891, b_s=0, w=0 equals 2", abs(margin_large_p1(aa, 0.0, 0.0) - 2.0) < 1e-12),
    ]
    allok = True
    for name, ok in checks:
        allok = allok and ok
        line("    %s  %s" % ("PASS" if ok else "FAIL", name))
    line("")
    line("ALL CHECKS PASS" if allok else "SOME CHECKS FAILED")


if __name__ == "__main__":
    main()
