"""
Quantifying the asserted inequality behind (A.56) / Lemma 4.9 of
"Finite time blowup for Navier-Stokes" (OpenAI, compiled 2026-09-08), page 36 and page 143.

Lemma 4.9 (p. 36) asserts 2 + h < a <= 2 + 2h on the outer collar e^{1/2} X_tail <= X < X_b.
With H = sqrt(2X) E and E = c_inf X^{-A} Hcal(2d/X), A = 1/2 + h, d = 1 - eta^2, and
D_X Z = -Z for Z = 2d/X, the shear coefficient a = 2 - 2 D_X log H satisfies exactly

    a = 2 + 2h + 2 Z Hcal'(Z)/Hcal(Z),          m(Z,h) := -Z Hcal'(Z)/Hcal(Z) > 0,

so a <= 2 + 2h always, and a > 2 + h holds exactly when m(Z,h) < h/2.  The appendix
asserts the stronger m(Z,h) < h/4 "for large X_R" without displaying a threshold.

This script finds, for a range of h, the largest Z with m(Z,h) <= h/4 and with
m(Z,h) <= h/2, and reports the corresponding lower bound on X = 2d/Z (worst case d = 1).
It also checks the small-Z asymptotics m(Z,h) = h Z + O(Z^2).

Run:  ./venv/bin/python commons/ns/artifacts/ns2026_heat_shear_threshold.py
Output: ns2026_heat_shear_threshold_output.txt
"""

import sys
import numpy as np
import scipy
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import gamma as Gamma

OUT = []


def emit(s=""):
    OUT.append(s)
    print(s)


def H_and_Hp(Z, h):
    def i0(v):
        return np.exp(-v) * v ** h * (1.0 + Z * v) ** (-h)

    def i1(v):
        return np.exp(-v) * v ** (h + 1.0) * (1.0 + Z * v) ** (-h - 1.0)

    a0 = quad(i0, 0.0, np.inf, limit=400)[0] / Gamma(1.0 + h)
    a1 = -h * quad(i1, 0.0, np.inf, limit=400)[0] / Gamma(1.0 + h)
    return a0, a1


def m(Z, h):
    H, Hp = H_and_Hp(Z, h)
    return -Z * Hp / H


def main():
    emit("environment: python %s   numpy %s   scipy %s"
         % (sys.version.split()[0], np.__version__, scipy.__version__))
    emit("")
    emit("m(Z,h) = -Z Hcal'(Z)/Hcal(Z);   a = 2 + 2h - 2 m(Z,h);   Z = 2d/X, d = 1-eta^2 <= 1")
    emit("")
    emit("small-Z asymptotics  m/(h Z) -> 1 :")
    for h in (1e-3, 1e-2, 0.05, 0.2):
        row = ["h=%-7.4g" % h]
        for Z in (1e-4, 1e-3, 1e-2):
            row.append("m/(hZ)@Z=%.0e: %.6f" % (Z, m(Z, h) / (h * Z)))
        emit("  " + "   ".join(row))
    emit("")
    emit("thresholds Z* where m(Z,h) = h/4 and = h/2, and the implied X >= 2/Z* (d=1):")
    emit("  %-10s %-14s %-14s %-14s %-14s" % ("h", "Z*(h/4)", "X>=2/Z*", "Z*(h/2)", "X>=2/Z*"))
    for h in (1e-4, 1e-3, 1e-2, 0.05, 0.2, 0.45):
        z4 = brentq(lambda Z: m(Z, h) - h / 4.0, 1e-8, 1e4, xtol=1e-10, rtol=1e-12)
        z2 = brentq(lambda Z: m(Z, h) - h / 2.0, 1e-8, 1e4, xtol=1e-10, rtol=1e-12)
        emit("  %-10.4g %-14.8f %-14.6f %-14.8f %-14.6f" % (h, z4, 2.0 / z4, z2, 2.0 / z2))
    emit("")
    emit("Reading: the threshold is essentially independent of h (Z* about 0.45 for m=h/4")
    emit("and about 1.0 for m=h/2), so the asserted inequality -Z Hcal'/Hcal < h/4 on the")
    emit("outer collar is the concrete requirement X > about 4.4 d, i.e. a fixed lower bound")
    emit("on the outer radius, NOT a further smallness condition on h.  Since the paper")
    emit("places the collar at X >= e^{1/2} X_tail with X_tail > X_v > X_R = 110 (C P_*)^{10}")
    emit("and P_* > e^{e^{M_d}+10}, the requirement is met with enormous room.")

    with open(__file__.replace(".py", "_output.txt"), "w") as fh:
        fh.write("\n".join(OUT) + "\n")


if __name__ == "__main__":
    main()
