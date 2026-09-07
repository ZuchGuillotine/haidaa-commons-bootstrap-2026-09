"""
Isotropic-equivalent radio energy scaling for fast radio bursts.

Purpose
-------
Put the 2020 April 28 radio burst from the Galactic magnetar SGR 1935+2154 on
the same isotropic-energy axis as cosmological FRBs, so that the size of the
energy gap can be stated with its input assumptions made explicit. This is a
transparent unit conversion, NOT a new measurement and NOT a physical model.

Formula
-------
    E_iso = 4 * pi * D_L**2 * F_nu * B_eff / (1 + z)

where F_nu is the specific fluence (erg cm^-2 Hz^-1), B_eff an assumed
effective emission bandwidth (Hz), D_L the luminosity distance (cm), z the
redshift. The (1+z) factor accounts for the bandwidth being measured in the
observer frame. Different papers adopt different bandwidth conventions
(observing bandwidth, burst bandwidth, or nu_c); results scale linearly with
B_eff, so the choice matters at the factor-of-a-few level.

Cosmology: flat LCDM, H0 = 67.7 km/s/Mpc, Omega_m = 0.31 (Planck-like).
Comoving distance integrated numerically with scipy.integrate.quad.

Determinism: no random numbers are used anywhere in this script, so no seed is
required; output is fully determined by the constants below.

Run:
    python frb_energetics_scaling.py > frb_energetics_scaling_out.txt
"""

import sys

import numpy as np
import scipy
from scipy.integrate import quad

# ---------------------------------------------------------------- constants
C_KMS = 299792.458          # speed of light, km/s
H0 = 67.7                   # km/s/Mpc
OMEGA_M = 0.31
OMEGA_L = 1.0 - OMEGA_M
MPC_CM = 3.0856775814913673e24   # cm per Mpc
KPC_CM = 3.0856775814913673e21   # cm per kpc
JY_CGS = 1.0e-23            # erg s^-1 cm^-2 Hz^-1 per Jy


def comoving_distance_mpc(z):
    """Line-of-sight comoving distance in Mpc for flat LCDM."""
    if z <= 0.0:
        return 0.0
    integrand = lambda zp: 1.0 / np.sqrt(OMEGA_M * (1.0 + zp) ** 3 + OMEGA_L)
    val, _err = quad(integrand, 0.0, z, epsabs=1e-10, epsrel=1e-10)
    return (C_KMS / H0) * val


def luminosity_distance_cm(z):
    return (1.0 + z) * comoving_distance_mpc(z) * MPC_CM


def e_iso(fluence_jy_ms, distance_cm, bandwidth_hz, z):
    """Isotropic-equivalent radio energy in erg."""
    f_nu = fluence_jy_ms * 1.0e-3 * JY_CGS      # Jy ms -> erg cm^-2 Hz^-1
    return 4.0 * np.pi * distance_cm ** 2 * f_nu * bandwidth_hz / (1.0 + z)


# ------------------------------------------------------------------- cases
# Each case: (label, fluence [Jy ms], distance spec, bandwidth [Hz], note)
# Fluence values are the approximate published/representative values named in
# the note; they are inputs to this conversion, not results of it.

CASES = [
    dict(label="SGR 1935+2154 (2020-04-28 radio burst), d = 9 kpc",
         fluence_jy_ms=1.5e6, z=0.0, dist_cm=9.0 * KPC_CM, band_hz=1.0e9,
         note="fluence ~1.5 MJy ms (STARE2, Bochenek+2020); distance poorly known"),
    dict(label="SGR 1935+2154, d = 4.4 kpc (low distance estimate)",
         fluence_jy_ms=1.5e6, z=0.0, dist_cm=4.4 * KPC_CM, band_hz=1.0e9,
         note="same fluence, lower distance bound sometimes quoted"),
    dict(label="SGR 1935+2154, d = 12.5 kpc (high distance estimate)",
         fluence_jy_ms=1.5e6, z=0.0, dist_cm=12.5 * KPC_CM, band_hz=1.0e9,
         note="same fluence, higher distance bound sometimes quoted"),
    dict(label="Nearby repeater in M81 (FRB 20200120E), d = 3.6 Mpc, F = 1 Jy ms",
         fluence_jy_ms=1.0, z=0.0, dist_cm=3.6 * MPC_CM, band_hz=1.0e9,
         note="representative bright burst; M81 distance ~3.6 Mpc"),
    dict(label="Modest FRB at z = 0.1, F = 1 Jy ms",
         fluence_jy_ms=1.0, z=0.1, dist_cm=None, band_hz=1.0e9,
         note="representative low-z cosmological FRB"),
    dict(label="Typical FRB at z = 0.5, F = 1 Jy ms",
         fluence_jy_ms=1.0, z=0.5, dist_cm=None, band_hz=1.0e9,
         note="representative cosmological FRB"),
    dict(label="Bright FRB at z = 0.5, F = 10 Jy ms",
         fluence_jy_ms=10.0, z=0.5, dist_cm=None, band_hz=1.0e9,
         note="bright end of the common population"),
    dict(label="FRB 20220610A-like at z = 1.016, F = 45 Jy ms",
         fluence_jy_ms=45.0, z=1.016, dist_cm=None, band_hz=1.0e9,
         note="z from Ryder+2023; fluence approximate, VERIFY against source"),
]


def main():
    print("frb_energetics_scaling.py -- isotropic radio energy conversions")
    print("python: %s" % sys.version.split()[0])
    print("numpy: %s   scipy: %s" % (np.__version__, scipy.__version__))
    print("cosmology: flat LCDM, H0=%.1f, Om=%.2f" % (H0, OMEGA_M))
    print("assumed effective bandwidth: 1.0e9 Hz for every case")
    print("no random numbers used; output is deterministic")
    print("")
    print("%-62s %12s %14s" % ("case", "D_L (Mpc)", "E_iso (erg)"))
    print("-" * 92)

    results = {}
    for c in CASES:
        z = c["z"]
        d_cm = c["dist_cm"] if c["dist_cm"] is not None else luminosity_distance_cm(z)
        e = e_iso(c["fluence_jy_ms"], d_cm, c["band_hz"], z)
        results[c["label"]] = e
        print("%-62s %12.4g %14.3e" % (c["label"], d_cm / MPC_CM, e))
        print("    note: %s" % c["note"])

    print("")
    print("ratios (dimensionless, same bandwidth convention throughout):")
    sgr = results["SGR 1935+2154 (2020-04-28 radio burst), d = 9 kpc"]
    typ = results["Typical FRB at z = 0.5, F = 1 Jy ms"]
    m81 = results["Nearby repeater in M81 (FRB 20200120E), d = 3.6 Mpc, F = 1 Jy ms"]
    hiz = results["FRB 20220610A-like at z = 1.016, F = 45 Jy ms"]
    print("  E(typical z=0.5 FRB) / E(SGR 1935 burst)      = %.3e  (%.2f dex)"
          % (typ / sgr, np.log10(typ / sgr)))
    print("  E(M81 repeater burst) / E(SGR 1935 burst)     = %.3e  (%.2f dex)"
          % (m81 / sgr, np.log10(m81 / sgr)))
    print("  E(z=1.016, 45 Jy ms) / E(SGR 1935 burst)      = %.3e  (%.2f dex)"
          % (hiz / sgr, np.log10(hiz / sgr)))
    print("")
    print("caveats: E_iso scales as d^2 and linearly with assumed bandwidth;")
    print("the SGR 1935+2154 distance is uncertain by roughly a factor 3 in d,")
    print("i.e. ~1 dex in energy. These are order-of-magnitude framings only.")


if __name__ == "__main__":
    main()
