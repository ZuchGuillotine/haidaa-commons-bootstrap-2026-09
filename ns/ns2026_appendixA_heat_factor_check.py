"""
Independent numerical check of the Appendix A heat-factor identities of the preprint
"Finite time blowup for Navier-Stokes" (author line OpenAI, compiled 2026-09-08),
pages 137-139 and 142-143.

Checked here, for several values of h in (0, 1/2):

  (A.32)  Hcal(Z) = (1/Gamma(a_K)) int_0^inf e^{-v} v^{a_K-1} (1+Zv)^{-h} dv,  a_K = 1+h
  (A.34)  Hcal^{(m)}(Z) = ((-1)^m (h)_m / Gamma(1+h)) int_0^inf e^{-v} v^{h+m} (1+Zv)^{-h-m} dv
  (A.35)  Hcal^{(m)}(0) = (-1)^m (h)_m (1+h)_m
  (A.36)  Hcal(Z) = 1 - h(1+h) Z + O(Z^2)
  (A.37)  Z^2 Hcal'' + (1 + 2 a_K Z) Hcal' + a_K(a_K-1) Hcal = 0
          two-sided bound 0 <= -Z Hcal'/Hcal < h
  (A.33)  K(r,t) = c_inf s^{-A} Hcal(2 tau / s), s = r^2/2, tau = 1-t, A = 1/2+h,
          solves d_t K = (d_rr + r^{-1} d_r - r^{-2}) K, and K_r < 0.

The point of the exercise: the choice a_K = 1 + h is forced by A = 1/2 + h through the
two exponent identities 2A+1 = 2 a_K and A^2 - 1/4 = a_K(a_K-1); the radial heat
operator then reduces exactly to the confluent-type ODE (A.37).

No claim of proof is made here; this is a numerical consistency check of algebra that
was also done by hand.
"""

import numpy as np
from scipy.integrate import quad
from scipy.special import gamma as Gamma

HS = [0.001, 0.01, 0.05, 0.2, 0.45]
ZS = [0.0, 1e-3, 0.05, 0.5, 2.0, 10.0]


def hcal(Z, h):
    aK = 1.0 + h
    f = lambda v: np.exp(-v) * v ** (aK - 1.0) * (1.0 + Z * v) ** (-h)
    val, err = quad(f, 0.0, np.inf, limit=400)
    return val / Gamma(aK), err


def hcal_deriv(Z, h, m):
    """(A.34): m-th derivative of Hcal."""
    if m == 0:
        return hcal(Z, h)[0]
    poch_h = np.prod([h + i for i in range(m)])          # (h)_m
    f = lambda v: np.exp(-v) * v ** (h + m) * (1.0 + Z * v) ** (-h - m)
    val, _ = quad(f, 0.0, np.inf, limit=400)
    return ((-1.0) ** m) * poch_h * val / Gamma(1.0 + h)


def poch(x, m):
    return np.prod([x + i for i in range(m)])


print("environment: python %s, numpy %s, scipy %s"
      % (".".join(map(str, __import__("sys").version_info[:3])),
         np.__version__, __import__("scipy").__version__))
print()

print("== Hcal(0) = 1 ==")
for h in HS:
    v, _ = hcal(0.0, h)
    print("  h=%.4g  Hcal(0)=%.12f  |Hcal(0)-1|=%.2e" % (h, v, abs(v - 1.0)))
print()

print("== (A.34) vs finite differences of (A.32), m=1,2 ==")
for h in HS:
    for Z in [0.05, 0.5, 2.0]:
        e = 1e-5
        d1_fd = (hcal(Z + e, h)[0] - hcal(Z - e, h)[0]) / (2 * e)
        d2_fd = (hcal(Z + e, h)[0] - 2 * hcal(Z, h)[0] + hcal(Z - e, h)[0]) / e ** 2
        d1 = hcal_deriv(Z, h, 1)
        d2 = hcal_deriv(Z, h, 2)
        print("  h=%.4g Z=%.4g  H'=%.10g (fd %.10g, dif %.1e)  H''=%.10g (fd %.10g, dif %.1e)"
              % (h, Z, d1, d1_fd, abs(d1 - d1_fd), d2, d2_fd, abs(d2 - d2_fd)))
print()

print("== (A.35) Hcal^{(m)}(0) = (-1)^m (h)_m (1+h)_m, m=1,2,3 ==")
for h in HS:
    for m in [1, 2, 3]:
        num = hcal_deriv(0.0, h, m)
        exact = ((-1.0) ** m) * poch(h, m) * poch(1.0 + h, m)
        print("  h=%.4g m=%d  numeric=%.12g  formula=%.12g  dif=%.2e"
              % (h, m, num, exact, abs(num - exact)))
print()

print("== (A.36) Hcal(Z) = 1 - h(1+h) Z + O(Z^2) ==")
for h in HS:
    for Z in [1e-4, 1e-3, 1e-2]:
        v, _ = hcal(Z, h)
        lin = 1.0 - h * (1.0 + h) * Z
        print("  h=%.4g Z=%.1e  Hcal=%.12f  linear=%.12f  |dif|/Z^2=%.6g"
              % (h, Z, v, lin, abs(v - lin) / Z ** 2))
print()

print("== (A.37) residual Z^2 H'' + (1 + 2 a_K Z) H' + a_K(a_K-1) H ==")
worst = 0.0
for h in HS:
    aK = 1.0 + h
    for Z in ZS:
        H = hcal(Z, h)[0]
        H1 = hcal_deriv(Z, h, 1)
        H2 = hcal_deriv(Z, h, 2)
        res = Z ** 2 * H2 + (1.0 + 2.0 * aK * Z) * H1 + aK * (aK - 1.0) * H
        worst = max(worst, abs(res))
        print("  h=%.4g Z=%.4g  residual=%+.3e   (scale H=%.4g)" % (h, Z, res, H))
print("  worst |residual| over the grid: %.3e" % worst)
print()

print("== two-sided bound 0 <= -Z Hcal'/Hcal < h ==")
for h in HS:
    for Z in ZS:
        H = hcal(Z, h)[0]
        H1 = hcal_deriv(Z, h, 1)
        q = -Z * H1 / H
        print("  h=%.4g Z=%.4g  -Z H'/H = %.10g   in [0,h)? %s"
              % (h, Z, q, (0.0 <= q < h)))
print()

print("== exponent identities forcing a_K = 1 + h, given A = 1/2 + h ==")
for h in HS:
    A = 0.5 + h
    aK = 1.0 + h
    print("  h=%.4g   2A+1=%.12g  2a_K=%.12g   A^2-1/4=%.12g  a_K(a_K-1)=%.12g"
          % (h, 2 * A + 1, 2 * aK, A ** 2 - 0.25, aK * (aK - 1.0)))
print()

print("== (A.33) K = s^{-A} Hcal(2 tau/s) solves d_t K = (d_rr + r^{-1} d_r - r^{-2}) K ==")


def K(r, t, h):
    s = r * r / 2.0
    tau = 1.0 - t
    A = 0.5 + h
    return s ** (-A) * hcal(2.0 * tau / s, h)[0]


for h in [0.01, 0.2, 0.45]:
    for (r, t) in [(0.7, 0.3), (1.0, 0.5), (2.0, 0.9), (3.5, 0.1)]:
        er, et = 1e-4, 1e-5
        dt = (K(r, t + et, h) - K(r, t - et, h)) / (2 * et)
        drr = (K(r + er, t, h) - 2 * K(r, t, h) + K(r - er, t, h)) / er ** 2
        dr = (K(r + er, t, h) - K(r - er, t, h)) / (2 * er)
        rhs = drr + dr / r - K(r, t, h) / r ** 2
        print("  h=%.4g r=%.3g t=%.3g  d_tK=%+.8g  Lap K=%+.8g  rel dif=%.2e  K_r=%+.4g (<0? %s)"
              % (h, r, t, dt, rhs, abs(dt - rhs) / max(abs(dt), 1e-300), dr, dr < 0))
