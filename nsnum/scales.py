"""
scales.py -- scale-separation calculator for the 2026 forced Navier-Stokes
blow-up construction (leading-order layer), in the notation of the paper's
similarity variables.

Reproduces three tables of the CFD feasibility assessment:
  (0.2) the tau at which each small parameter of the construction first
        reaches 0.1, and the state of the flow at tau = 1e-6;
  (1.2) closed-form resolution requirements N_r, N_z, N_theta, N_t as a
        function of tau = 10^-k;
  (1.5) the exponent-discrimination table (amplitude excess over Leray
        scaling per decade of dynamic range) and the surrogate-h table.

Definitions used (all from the leading-order scalings; q ~ tau on compact
sets of similarity coordinates away from eta = +-1):

  tau            = 1 - t, time to the singular time
  s              = -ln tau, the similarity time
  A = 1/2 + h    velocity growth exponent, |u|_inf ~ tau^-(1/2+h)
  D = 1/2 - h    axial length exponent, l_z ~ tau^(1/2-h)
  l_r ~ tau^(1/2)                     core radial scale
  l_z ~ tau^(1/2-h)                   core axial scale
  aspect  = l_z/l_r = tau^-h
  Re_theta = tau^-h                   azimuthal Reynolds number
  q^(2h)  = tau^(2h)                  background-correction expansion parameter
  q^(h/2) = tau^(h/2)                 pulse-to-core scale ratio l_wave/l_r
  waves per core radius = tau^(-h/2)

NOTHING here is a simulation and nothing here tests whether the theorem is
correct. These are algebraic consequences of the stated exponents.

Run:
  <venv>/bin/python scales.py
Output is deterministic; there is no random state.
"""

import sys

import numpy as np

H_PROOF = 1.0 / 200.0  # h = 1/200, the midpoint-ish value used throughout
                       # the assessment; the proof requires 0 < h < 1/100.


def sep(title):
    print()
    print("=" * 78)
    print(title)
    print("=" * 78)


def tau_of_k(k):
    """tau = 10^-k, returned as the exponent only (tau itself underflows)."""
    return -k


def k_for_target(exponent, target):
    """Smallest k with 10^(-k*exponent) <= target, for exponent > 0."""
    return np.log10(1.0 / target) / exponent


def table_0_2():
    sep("SECTION 0.2 -- where the asymptotic regime begins (h = 1/200)")
    h = H_PROOF
    print("Small parameters of the construction, all powers of q^h with q ~ tau.")
    print()
    print(f"{'target':<44}{'requires tau =':>18}{'s = -ln tau':>16}")
    rows = [
        ("q^(2h) = 0.1  (background expansion small)", k_for_target(2 * h, 0.1)),
        ("Re_theta = tau^-h = 10", k_for_target(h, 0.1)),
        ("q^(h/2) = 0.1 (pulses genuinely sub-core)", k_for_target(h / 2, 0.1)),
    ]
    for label, k in rows:
        print(f"{label:<44}{'10^-' + format(k, '.0f'):>18}{k * np.log(10.0):>16.0f}")
    print()
    print("IEEE-754 double underflow (smallest normal) is about 1e-308, so")
    print("tau = 1e-400 cannot be stored as a double at all; tau = 1e-100 and")
    print("tau = 1e-400 are s = 230 and s = 921, both exactly ordinary numbers.")
    print()
    print("State of the flow at the deepest tau a physical-space time-marching")
    print("computation could plausibly reach, tau = 1e-6, at h = 1/200:")
    k = 6.0
    print(f"  Re_theta = tau^-h                 = {10 ** (k * h):.4f}")
    print(f"  aspect ratio l_z/l_r = tau^-h     = {10 ** (k * h):.4f}")
    print(f"  q^(2h)  (expansion parameter)     = {10 ** (-k * 2 * h):.4f}")
    print(f"  q^(h/2) = l_wave/l_r              = {10 ** (-k * h / 2):.4f}")
    print(f"  waves per core radius             = {10 ** (k * h / 2):.4f}")
    print(f"  amplitude excess over tau^-1/2    = {100 * (10 ** (k * h) - 1):.1f} %")
    print()
    print("i.e. an O(1)-Reynolds-number viscous vortex with about one carrier")
    print("wavelength across the core: no scale separation of any kind.")


def table_1_2():
    sep("SECTION 1.2 -- resolution requirements vs tau (h = 1/200)")
    h = H_PROOF
    print("N_r = n*tau^-(1/2)          core radial")
    print("N_z = n*tau^-(1/2-h)        core axial")
    print("N_pulse = n*tau^-(1/2+h/2)  pulse-carrier resolved (both directions)")
    print("N_theta = 3*m_max, m_max ~ 2*pi*tau^(-h/2)   [INFERRED, see below]")
    print("N_t = 2n*(finest resolved inverse scale)     advective CFL, IMEX diffusion")
    print()
    for n, label in ((16, "n = 16 points per scale (spectral)"),
                     (48, "n = 48 points per scale (2nd-order finite volume)")):
        print(label)
        print(f"  {'tau':>8}{'N_r':>12}{'N_z':>12}{'N_pulse':>12}"
              f"{'N_theta':>10}{'N_t (core)':>13}{'N_t (pulse)':>13}")
        for k in (2.0, 4.0, 6.0, 8.0):
            n_r = n * 10 ** (k * 0.5)
            n_z = n * 10 ** (k * (0.5 - h))
            n_p = n * 10 ** (k * (0.5 + h / 2))
            m_max = 2 * np.pi * 10 ** (k * h / 2)
            n_th = 3 * m_max
            n_t_core = 2 * n * 10 ** (k * 0.5)
            n_t_pulse = 2 * n * 10 ** (k * (0.5 + h / 2))
            print(f"  {'1e-' + format(k, '.0f'):>8}{n_r:>12.2e}{n_z:>12.2e}"
                  f"{n_p:>12.2e}{n_th:>10.1f}{n_t_core:>13.2e}{n_t_pulse:>13.2e}")
        print()
    print("Exterior: the exterior swirl obeys the 1D radial heat-type equation")
    print("-d_tau K = d_rr K + K_r/r - K/r^2, so it costs O(1e2) radial points as a")
    print("1D problem and must not be discretized as 3D Navier-Stokes.")
    print()
    print("CAVEAT on N_theta: the azimuthal estimate assumes the carrier is")
    print("isotropic, m ~ 2*pi*r/l_wave ~ tau^(-h/2), giving m_max = 6.5 at")
    print("tau = 1e-6. That is an INFERENCE, not a quotation from the paper.")
    print("If instead m ~ 1/l_wave ~ tau^-(1/2+h/2), N_theta and the 3D cost grow")
    print("by a further factor tau^-(1/2), i.e. about 1e3 at tau = 1e-6.")


def table_1_5():
    sep("SECTION 1.5 -- exponent discrimination and the surrogate-h table")
    h = H_PROOF
    print("The numerical claim to be discriminated is |u|_inf ~ tau^-(1/2+h)")
    print("against the ordinary Leray scaling tau^-(1/2). Over k decades of tau")
    print("the amplitude ratio is 10^(k*h).")
    print()
    print(f"  {'dynamic range':>16}{'amplitude excess':>20}{'log10 excess':>16}")
    for k in (2.0, 4.0, 6.0, 60.0):
        exc = 10 ** (k * h)
        print(f"  {format(k, '.0f') + ' decades':>16}"
              f"{100 * (exc - 1):>19.2f}%{k * h:>16.4f}")
    print()
    print("So at h = 1/200 the measurement is a log-log slope discrimination of")
    print(f"  0.5000 vs {0.5 + h:.4f}  -- a relative exponent difference of "
          f"{100 * h / 0.5:.1f}%.")
    print("Over 6 decades that is 0.030 in log10 amplitude, so a 5-sigma")
    print("measurement needs the TOTAL systematic error budget (spatial +")
    print("temporal truncation + outer-boundary/cutoff + initialization")
    print("transient + fitting-window bias) below about 0.3% in amplitude,")
    print("sustained over 6 decades.")
    print()
    print("Surrogate-h table, evaluated at tau = 1e-4:")
    print(f"  {'h':>8}{'Re_theta':>12}{'aspect':>10}{'q^(h/2)':>10}"
          f"{'waves/core':>13}{'slope 1/2+h':>13}")
    k = 4.0
    for hh in (1.0 / 200.0, 0.01, 0.05, 0.10, 0.15):
        print(f"  {hh:>8.3f}{10 ** (k * hh):>12.2f}{10 ** (k * hh):>10.2f}"
              f"{10 ** (-k * hh / 2):>10.2f}{10 ** (k * hh / 2):>13.2f}"
              f"{0.5 + hh:>13.3f}")
    print()
    print("h in (1/100, 1/6) is OUTSIDE the interval the proof's estimates")
    print("require (0 < h < 1/100) though inside the energy-integrability")
    print("constraint h < 1/6. A run at such h is a mechanism study, never a")
    print("verification of the theorem.")


def self_check():
    """Assert the values quoted in the assessment text, so a misread exponent
    is caught immediately rather than propagating into the cost model."""
    h = H_PROOF
    checks = [
        ("k for q^(2h)=0.1", k_for_target(2 * h, 0.1), 100.0, 1e-9),
        ("k for Re_theta=10", k_for_target(h, 0.1), 200.0, 1e-9),
        ("k for q^(h/2)=0.1", k_for_target(h / 2, 0.1), 400.0, 1e-9),
        ("Re_theta at tau=1e-6", 10 ** (6 * h), 1.0715, 1e-3),
        ("q^(h/2) at tau=1e-6", 10 ** (-6 * h / 2), 0.9661, 1e-3),
        ("amplitude excess at 6 decades", 10 ** (6 * h) - 1, 0.0715, 1e-3),
        ("amplitude excess at 60 decades", 10 ** (60 * h) - 1, 0.9953, 1e-3),
        ("s at tau=1e-100", 100 * np.log(10.0), 230.26, 1e-2),
        ("s at tau=1e-400", 400 * np.log(10.0), 921.03, 1e-2),
        ("Re_theta at h=0.15, tau=1e-4", 10 ** (4 * 0.15), 3.981, 1e-3),
    ]
    sep("SELF-CHECK against the values quoted in the assessment")
    ok = True
    for name, got, want, tol in checks:
        good = abs(got - want) <= tol * max(1.0, abs(want))
        ok = ok and good
        print(f"  {'PASS' if good else 'FAIL'}  {name:<34} "
              f"computed {got:>12.5f}  quoted {want:>10.5f}")
    print()
    print("ALL CHECKS PASS" if ok else "SOME CHECKS FAILED")
    return ok


def main():
    print("scales.py -- scale separation for the 2026 forced Navier-Stokes")
    print("blow-up construction, leading-order scalings only.")
    print(f"python {sys.version.split()[0]}  numpy {np.__version__}")
    print("No random state; output is deterministic.")
    table_0_2()
    table_1_2()
    table_1_5()
    ok = self_check()
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
