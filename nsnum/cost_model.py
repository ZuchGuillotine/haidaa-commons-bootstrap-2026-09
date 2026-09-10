"""
cost_model.py -- explicit order-of-magnitude cost model for simulating the
leading-order layer of the 2026 forced Navier-Stokes blow-up construction,
in physical space and in the construction's own similarity variables.

Reproduces the four cost tables of section 1.4 of the CFD feasibility
assessment:
  (A) axisymmetric, physical space, core-resolved (pulses unresolved)
  (B) axisymmetric, physical space, pulse-carrier resolved
  (C) full 3D, physical space, pulse-resolved, N_theta = 32
  (D) similarity variables (X, eta, theta, s = -ln tau), fixed grid

EVERY NUMBER PRINTED HERE IS AN ESTIMATE FROM THE MODEL BELOW, NOT A
MEASUREMENT. The model is stated so that it can be disagreed with
quantitatively. Nothing here is a simulation, and nothing here evaluates
whether the theorem is correct.

Model assumptions (assessment section 1.3):
  - 300 flop/point/step for an incompressible solver with IMEX diffusion and
    a multigrid or FFT pressure solve;
  - 1000 flop/point/step for implicit sparse-spectral (Dedalus-style) with
    transforms, used for the similarity-variable rows;
  - 1 GFlop/s sustained fp64 per server core (memory-bandwidth-bound stencil
    and multigrid work; peak is 30-100x higher and irrelevant);
  - 1 TFlop/s sustained fp64 per data-centre GPU, so 1 GPU-hour is about
    1e3 core-hours for well-tuned spectral-element or FFT kernels;
  - 10 doubles per grid point resident (u, p, RHS, work arrays);
  - snapshot = 5 fields * 8 B * cells, 200 snapshots per run;
  - no AMR in the physical-space rows. The similarity rows already capture
    the AMR benefit exactly, and better, because the map is known
    analytically.

Run:
  <venv>/bin/python cost_model.py
Output is deterministic; there is no random state.
"""

import sys

import numpy as np

H_PROOF = 1.0 / 200.0
N_SPECTRAL = 16          # points per finest resolved structure scale
FLOP_FV = 300.0          # flop / point / step, IMEX incompressible solver
FLOP_SPECTRAL = 1000.0   # flop / point / step, implicit sparse spectral
CORE_FLOPS = 1.0e9       # sustained fp64 per core
GPU_PER_CORE = 1.0e3     # core-hours per GPU-hour
BYTES_RESIDENT = 10 * 8  # 10 doubles per point
SNAP_BYTES = 5 * 8       # 5 fields per snapshot point
N_SNAP = 200
N_THETA_3D = 32          # 3D azimuthal modes, from m_max ~ 6.5 rounded up


def core_hours(cells, steps, flop_per_point_step):
    return cells * steps * flop_per_point_step / CORE_FLOPS / 3600.0


def gb(nbytes):
    return nbytes / 1e9


def fmt_bytes(nbytes):
    if nbytes >= 1e12:
        return f"{nbytes / 1e12:.2g} TB"
    if nbytes >= 1e9:
        return f"{nbytes / 1e9:.2g} GB"
    return f"{nbytes / 1e6:.2g} MB"


def sep(title):
    print()
    print("=" * 84)
    print(title)
    print("=" * 84)


def table_A(h=H_PROOF, n=N_SPECTRAL):
    sep("(A) Axisymmetric, physical space, core-resolved (pulses unresolved) [EST]")
    print(f"{'tau':>8}{'N_r':>10}{'N_z':>10}{'cells':>12}{'steps':>10}"
          f"{'core-hours':>13}{'memory':>12}")
    for k in (2.0, 4.0, 6.0):
        n_r = n * 10 ** (k * 0.5)
        n_z = n * 10 ** (k * (0.5 - h))
        cells = n_r * n_z
        steps = 2 * n * 10 ** (k * 0.5)
        ch = core_hours(cells, steps, FLOP_FV)
        print(f"{'1e-' + format(k, '.0f'):>8}{n_r:>10.1e}{n_z:>10.1e}"
              f"{cells:>12.1e}{steps:>10.1e}{ch:>13.1e}"
              f"{fmt_bytes(cells * BYTES_RESIDENT):>12}")


def table_B(h=H_PROOF, n=N_SPECTRAL):
    sep("(B) Axisymmetric, physical space, pulse-carrier resolved [EST]")
    print(f"{'tau':>8}{'N_r = N_z':>12}{'cells':>12}{'steps':>10}"
          f"{'core-hours':>13}{'memory':>12}")
    for k in (2.0, 4.0, 6.0):
        n_p = n * 10 ** (k * (0.5 + h / 2))
        cells = n_p * n_p
        steps = 2 * n * 10 ** (k * (0.5 + h / 2))
        ch = core_hours(cells, steps, FLOP_FV)
        print(f"{'1e-' + format(k, '.0f'):>8}{n_p:>12.1e}{cells:>12.1e}"
              f"{steps:>10.1e}{ch:>13.1e}"
              f"{fmt_bytes(cells * BYTES_RESIDENT):>12}")


def table_C(h=H_PROOF, n=N_SPECTRAL):
    sep(f"(C) Full 3D, physical space, pulse-resolved, N_theta = {N_THETA_3D} [EST]")
    print("N_theta = 32 comes from m_max ~ 2*pi*tau^(-h/2) = 6.5 at tau = 1e-6,")
    print("rounded up; that m-scaling is an INFERENCE, not a quotation from the")
    print("paper. If m ~ 1/l_wave instead, every 3D figure below is low by ~1e3.")
    print()
    print(f"{'tau':>8}{'cells':>12}{'steps':>10}{'core-hours':>13}"
          f"{'GPU-hours':>12}{'memory':>12}{'snapshots':>12}")
    for k in (2.0, 4.0, 6.0):
        n_p = n * 10 ** (k * (0.5 + h / 2))
        cells = n_p * n_p * N_THETA_3D
        steps = 2 * n * 10 ** (k * (0.5 + h / 2))
        ch = core_hours(cells, steps, FLOP_FV)
        print(f"{'1e-' + format(k, '.0f'):>8}{cells:>12.1e}{steps:>10.1e}"
              f"{ch:>13.1e}{ch / GPU_PER_CORE:>12.1e}"
              f"{fmt_bytes(cells * BYTES_RESIDENT):>12}"
              f"{fmt_bytes(cells * SNAP_BYTES * N_SNAP):>12}")


def table_D(h=H_PROOF, ds=1.0e-3):
    sep("(D) Similarity variables (X, eta, theta, s = -ln tau), fixed grid,"
        f" ds = {ds:g} [EST]")
    print("The grid is FIXED: the collapse is absorbed by the coordinates, so the")
    print("only cost that grows with depth is the number of s-steps, linearly in")
    print("s = k*ln(10), plus a slow growth of resolution as waves per core rise.")
    print()
    print(f"{'tau':>10}{'s_max':>8}{'waves/core':>12}{'grid':>12}{'steps':>10}"
          f"{'2D core-h':>12}{'3D core-h':>12}{'3D GPU-h':>11}")
    rows = ((6.0, 256, 128), (100.0, 256, 128), (400.0, 1280, 640))
    for k, n_x, n_eta in rows:
        s_max = k * np.log(10.0)
        waves = 10 ** (k * h / 2)
        cells = n_x * n_eta
        steps = s_max / ds
        ch2 = core_hours(cells, steps, FLOP_SPECTRAL)
        ch3 = ch2 * N_THETA_3D
        print(f"{'1e-' + format(k, '.0f'):>10}{s_max:>8.0f}{waves:>12.1f}"
              f"{format(n_x, 'd') + 'x' + format(n_eta, 'd'):>12}{steps:>10.1e}"
              f"{ch2:>12.1e}{ch3:>12.1e}{ch3 / GPU_PER_CORE:>11.1e}")


def inversion(h=H_PROOF, n=N_SPECTRAL, ds=1.0e-3):
    sep("THE COST INVERSION -- read row (D) against row (C)")
    k = 6.0
    n_p = n * 10 ** (k * (0.5 + h / 2))
    cells_c = n_p * n_p * N_THETA_3D
    steps_c = 2 * n * 10 ** (k * (0.5 + h / 2))
    ch_c = core_hours(cells_c, steps_c, FLOP_FV)

    k_d = 400.0
    cells_d = 1280 * 640
    steps_d = k_d * np.log(10.0) / ds
    ch_d = core_hours(cells_d, steps_d, FLOP_SPECTRAL) * N_THETA_3D

    print(f"Physical-space 3D to tau = 1e-6, where the construction's entire")
    print(f"signature is a {100 * (10 ** (k * h) - 1):.1f}% amplitude excess over an ordinary")
    print(f"Leray collapse and Re_theta = {10 ** (k * h):.2f}:")
    print(f"    {ch_c:.1e} core-hours  ({ch_c / GPU_PER_CORE:.1e} GPU-hours)")
    print(f"Similarity-variable 3D to tau = 1e-400 (s = {k_d * np.log(10.0):.0f}), where the")
    print(f"pulses are genuinely sub-core (q^(h/2) = {10 ** (-k_d * h / 2):.2f}) and the")
    print("background expansion parameter is small:")
    print(f"    {ch_d:.1e} core-hours  ({ch_d / GPU_PER_CORE:.1e} GPU-hours)")
    print()
    print(f"Ratio (useless regime) / (meaningful regime) = {ch_c / ch_d:.2f}")
    print("The regime that matters is CHEAPER than the regime that does not.")
    print("Caveat: this compares a physical-space 3D DNS with a similarity-variable")
    print("3D run of the leading-order system on a fixed grid; they are different")
    print("computations, and the comparison is a statement about where an")
    print("allocation should go, not about equal-fidelity cost.")


def self_check():
    """Assert the model reproduces the figures quoted in the assessment's
    section 1.4 tables, to within order-of-magnitude rounding."""
    h, n = H_PROOF, N_SPECTRAL
    sep("SELF-CHECK against the figures quoted in assessment section 1.4")
    results = []

    def add(name, got, want, rel):
        results.append((name, got, want, abs(got - want) <= rel * abs(want)))

    # (A) tau = 1e-6
    n_r = n * 10 ** 3
    n_z = n * 10 ** (6 * (0.5 - h))
    add("A k=6 cells", n_r * n_z, 2.4e8, 0.05)
    add("A k=6 core-hours", core_hours(n_r * n_z, 2 * n * 1e3, FLOP_FV), 6e2, 0.10)
    add("A k=6 memory GB", gb(n_r * n_z * BYTES_RESIDENT), 19.0, 0.10)
    # (B) tau = 1e-6
    n_p = n * 10 ** (6 * (0.5 + h / 2))
    add("B k=6 cells", n_p * n_p, 2.7e8, 0.05)
    add("B k=6 core-hours",
        core_hours(n_p * n_p, 2 * n * 10 ** (6 * (0.5 + h / 2)), FLOP_FV), 8e2, 0.10)
    # (C) tau = 1e-6
    cells_c = n_p * n_p * N_THETA_3D
    steps_c = 2 * n * 10 ** (6 * (0.5 + h / 2))
    add("C k=6 cells", cells_c, 8.8e9, 0.05)
    add("C k=6 core-hours", core_hours(cells_c, steps_c, FLOP_FV), 2.4e4, 0.10)
    add("C k=6 snapshot TB", cells_c * SNAP_BYTES * N_SNAP / 1e12, 70.0, 0.10)
    # (D)
    add("D k=400 s_max", 400 * np.log(10.0), 921.0, 0.01)
    add("D k=400 2D core-h",
        core_hours(1280 * 640, 400 * np.log(10.0) / 1e-3, FLOP_SPECTRAL), 2e2, 0.10)
    add("D k=400 3D core-h",
        core_hours(1280 * 640, 400 * np.log(10.0) / 1e-3, FLOP_SPECTRAL) * N_THETA_3D,
        7e3, 0.10)
    add("D k=100 3D core-h",
        core_hours(256 * 128, 100 * np.log(10.0) / 1e-3, FLOP_SPECTRAL) * N_THETA_3D,
        7e1, 0.10)

    ok = True
    for name, got, want, good in results:
        ok = ok and good
        print(f"  {'PASS' if good else 'FAIL'}  {name:<20} computed {got:>12.3e}"
              f"   quoted {want:>10.3e}")
    print()
    print("ALL CHECKS PASS" if ok else "SOME CHECKS FAILED")
    return ok


def main():
    print("cost_model.py -- order-of-magnitude cost model, h = 1/200.")
    print(f"python {sys.version.split()[0]}  numpy {np.__version__}")
    print("All figures are ESTIMATES from the stated model, not measurements.")
    table_A()
    table_B()
    table_C()
    table_D()
    inversion()
    ok = self_check()
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
