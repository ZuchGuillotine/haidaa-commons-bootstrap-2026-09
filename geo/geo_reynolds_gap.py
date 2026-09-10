"""
geo_reynolds_gap.py

Radial and azimuthal Reynolds numbers of real geophysical vortex cores, computed
from molecular viscosity, and the decade separation from the forced Navier-Stokes
blow-up construction's defining regime (radial Reynolds number of order one while
the azimuthal Reynolds number diverges).

Definitions used throughout (order-of-magnitude, core values):
    Re_r     = |u_r| * l_r / nu          radial Reynolds number
    Re_theta = |u_theta| * l_r / nu      azimuthal Reynolds number
    S_w      = |u_theta| / |u_r|         swirl ratio (core, not chamber, definition)
    gap      = log10(Re_r) - log10(Re_r_construction), with Re_r_construction = 1

l_r is the core radius (radius of maximum azimuthal wind); u_r is the mean radial
inflow speed at that radius; u_theta is the peak azimuthal speed. nu is the
MOLECULAR kinematic viscosity: this script deliberately does NOT use an eddy
viscosity. The eddy-viscosity case is a separate calculation (geo_re_eff_grid.py).

No fitting, no random numbers, no seeds needed: this is closed-form arithmetic on
tabulated scale estimates. The scale estimates themselves are order-of-magnitude
values and each is labelled with its basis.

Run:
    ./venv/bin/python commons/geo/artifacts/geo_reynolds_gap.py
"""

import sys
import numpy as np

NU_AIR = 1.5e-5     # m^2 / s, dry air near sea level, ~15 C
NU_WATER = 1.0e-6   # m^2 / s, sea water, ~20 C

# name, l_r [m], u_r [m/s], u_theta [m/s], nu [m^2/s], basis for the scales
CASES = [
    ("tornado core (violent, EF4-5)", 1.0e2, 1.0e1, 8.0e1, NU_AIR,
     "mobile-radar core radii of order 100 m; peak azimuthal wind 70-90 m/s; inflow order 10 m/s"),
    ("tornado core (weak, EF0-1)", 3.0e1, 3.0e0, 3.0e1, NU_AIR,
     "small-radius weak vortex; chosen as the least favourable air case for the gap"),
    ("tropical cyclone eyewall", 3.0e4, 1.0e1, 6.0e1, NU_AIR,
     "radius of maximum wind 30 km; boundary-layer inflow order 10 m/s; peak wind 60 m/s"),
    ("dust devil", 5.0e0, 2.0e0, 1.0e1, NU_AIR,
     "smallest coherent atmospheric vortex class; deliberately included as the minimum-Re case"),
    ("stratospheric polar vortex", 3.0e6, 1.0e0, 6.0e1, NU_AIR,
     "planetary-scale; nu at 10 hPa is larger than sea-level nu, so this Re_r is a LOWER bound"),
    ("mesoscale ocean eddy", 5.0e4, 1.0e-2, 5.0e-1, NU_WATER,
     "altimetry-scale eddy; molecular nu of sea water"),
]

# Sensitivity probe: how far would the scales have to be wrong for the gap to close
# to two decades (the threshold at which the seeding assessment said the framing
# would need revision)?
GAP_REVISION_THRESHOLD = 2.0


def main():
    print("geo_reynolds_gap.py")
    print("environment: python %s, numpy %s" % (sys.version.split()[0], np.__version__))
    print("command: ./venv/bin/python commons/geo/artifacts/geo_reynolds_gap.py")
    print("no random numbers used; no seed required")
    print("")
    print("Construction regime (OpenAI 2026, Finite time blowup for Navier-Stokes, pp. 4, 8):")
    print("  Re_r = O(1) held fixed as tau -> 0;  Re_theta ~ tau^{-h} -> infinity, h < 1/100")
    print("  reference value used for the gap: Re_r_construction = 1")
    print("")
    hdr = ("%-32s %10s %8s %10s %12s %12s %8s %8s"
           % ("case", "l_r [m]", "u_r", "u_theta", "Re_r", "Re_theta", "S_w", "decades"))
    print(hdr)
    print("-" * len(hdr))

    rows = []
    for name, lr, ur, ut, nu, basis in CASES:
        re_r = ur * lr / nu
        re_t = ut * lr / nu
        sw = ut / ur
        gap = np.log10(re_r)
        rows.append((name, lr, ur, ut, nu, re_r, re_t, sw, gap, basis))
        print("%-32s %10.3g %8.3g %10.3g %12.3g %12.3g %8.2f %8.2f"
              % (name, lr, ur, ut, re_r, re_t, sw, gap))

    gaps = np.array([r[8] for r in rows])
    print("")
    print("decade separation in Re_r from the construction's O(1) regime:")
    print("  minimum  %.2f decades  (%s)" % (gaps.min(), rows[int(np.argmin(gaps))][0]))
    print("  maximum  %.2f decades  (%s)" % (gaps.max(), rows[int(np.argmax(gaps))][0]))
    print("  median   %.2f decades" % float(np.median(gaps)))
    print("")

    print("basis for each scale estimate:")
    for r in rows:
        print("  %-32s %s" % (r[0], r[9]))
    print("")

    print("robustness: factor by which the product u_r*l_r would have to be REDUCED")
    print("for Re_r to fall to 10^%.0f (the threshold at which the seeding assessment"
          % GAP_REVISION_THRESHOLD)
    print("said the Common's framing would need revision toward small-scale vortices):")
    for r in rows:
        need = r[5] / (10.0 ** GAP_REVISION_THRESHOLD)
        print("  %-32s factor %.3g" % (r[0], need))
    print("")

    print("Re_theta / Re_r for each case (this equals the swirl ratio S_w by construction,")
    print("and is the quantity the blow-up construction drives to infinity as tau^{-h}):")
    for r in rows:
        print("  %-32s %.2f" % (r[0], r[6] / r[5]))
    print("")
    print("With h < 1/100, over a tornado intensification during which tau shrinks by a")
    print("factor of 60 (one hour to one minute), the construction's swirl ratio changes by")
    print("60^0.01 = %.4f, i.e. %.1f percent." % (60.0 ** 0.01, 100.0 * (60.0 ** 0.01 - 1.0)))
    print("Real cores sit at S_w of order 5-60 and do not evolve as a power of time-to-peak.")
    print("")
    print("CONCLUSION (arithmetic only, no claim about the proof): every geophysical vortex")
    print("class tabulated here sits between %.1f and %.1f decades above the construction's"
          % (gaps.min(), gaps.max()))
    print("radial Reynolds number. There is no literal transfer of the construction's regime")
    print("to any of them. Note the dust devil case is the closest at %.2f decades, slightly"
          % gaps.min())
    print("below the '6 to 11 decades' range quoted in the seeding assessment; the assessment's")
    print("lower bound should be read as 6 only if dust devils are excluded.")


if __name__ == "__main__":
    main()
