#!/usr/bin/env python3
"""Reproduces the computed tables in hemodynamics.md (sections 1.2, 1.4b, 1.4c, 1.5).

Inputs: the blowup paper's stated exponents (pp. 4, 8) -- l_r ~ tau^(1/2),
l_z ~ tau^(1/2-h), |u_theta|,|u_z| ~ tau^(-1/2-h), h < 1/100 -- plus standard
blood properties. Nothing here is taken from the literature except the Carreau
parameters (Cho & Kensey-type, provenance unverified).
"""
import math

def anisotropy_reachability(h, target):
    """log10(tau) needed for aspect ratio l_z/l_r = tau^-h to reach `target`."""
    return -math.log10(target) / h

def damage_exponent(alpha, h):
    """Exponent of tau in the instantaneous volume-integrated damage rate.
    gamma_dot ~ tau^-(1+h), V ~ tau^(3/2-h)  =>  int tau_s^alpha dV ~ tau^(3/2-h-alpha(1+h)).
    Instantaneous integral diverges if < 0; cumulative time integral diverges if <= -1."""
    return 1.5 - h - alpha * (1.0 + h)

def alpha_crit_cumulative(h):
    """Critical hemolysis exponent for cumulative (time-integrated) divergence."""
    return (2.5 - h) / (1.0 + h)

def carreau(g, mu0=0.056, muinf=0.00345, lam=3.313, n=0.3568):
    return muinf + (mu0 - muinf) * (1.0 + (lam * g) ** 2) ** ((n - 1) / 2)

if __name__ == "__main__":
    print("1.2 anisotropy reachability (log10 tau required)")
    for h in (0.01, 0.001):
        row = "  h=%-7g " % h + " ".join(
            "AR%-5d:1e%-6.0f" % (t, anisotropy_reachability(h, t)) for t in (2, 10, 100, 4000))
        print(row)

    print("\n1.5 damage-functional divergence")
    for a in (1.991, 2.416, 2.5, 3.0):
        e = damage_exponent(a, 0.0)
        print("  alpha=%-6.3f exponent %+0.4f  instantaneous %-8s cumulative %s"
              % (a, e, "DIVERGES" if e < 0 else "bounded",
                 "DIVERGES" if e <= -1 else "converges"))
    for h in (0.0, 0.005, 0.01):
        print("  h=%-6g alpha_crit(cumulative) = %.4f" % (h, alpha_crit_cumulative(h)))

    print("\n1.4b Carreau viscosity vs shear rate")
    for g in (0.1, 1, 10, 100, 1000, 5000, 20000):
        mu = carreau(g)
        print("  gamma_dot=%-8g mu=%7.3f mPa.s  (%6.3f x mu_inf)" % (g, mu * 1e3, mu / 0.00345))

    print("\n1.4c length-scale floors (FDA nozzle throat D=4mm, nu=3.3e-6 m^2/s)")
    nu, D = 3.3e-6, 4e-3
    for Re in (500, 2000, 3500, 5000, 6500):
        print("  Re=%-6d U=%5.2f m/s  eta ~ D*Re^-3/4 = %6.1f um" % (Re, Re * nu / D, D * Re ** -0.75 * 1e6))
    print("  RBC diameter 8.0 um; continuum floor (5-10 RBC) = 40-80 um")
