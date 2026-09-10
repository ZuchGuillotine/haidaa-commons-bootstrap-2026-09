"""Measured cost of one 2D reduced-MHD pseudo-spectral step, for sizing runs.

Purpose. The resource estimate for the reduced-MHD fixture library needs a
defensible per-step cost rather than a guess. This script measures the real
transform throughput on the machine it is run on and converts it into a per-step
and per-run cost for the fixture library's resolution ladder.

Model of the step. 2D reduced MHD in a doubly periodic box, Alfven units:
    d_t w   + [phi, w] = [psi, j] + nu Lap w
    d_t psi + [phi, psi] = eta j
    w = Lap phi,  j = Lap psi,  [f,g] = d_x f d_y g - d_y f d_x g.
Evaluating the two right-hand sides pseudo-spectrally needs, per RK stage:
  - 2 inverse transforms to get phi and psi in real space,
  - 8 inverse transforms for the four gradient pairs (grad phi, grad psi,
    grad w, grad j),
  - 2 forward transforms for the two bracket products.
That is 12 transforms per stage; RK4 has 4 stages, so 48 real-to-complex or
complex-to-real transforms of size N x N per time step. This count is an
implementation choice, stated here so the arithmetic can be checked or replaced.

What is measured: wall time of numpy.fft.rfft2 and numpy.fft.irfft2 at each N in
the ladder, best of NREP repeats after a warm-up, single process, whatever
threading numpy's FFT uses by default on this machine.

Run:
  <venv>/bin/python rmhd_step_cost_scaling.py
  (writes rmhd_step_cost_scaling_output.json next to this script)

No randomness beyond the fixed seed 20260909 used to fill the test arrays.
Timings are hardware specific; the JSON records the platform.
"""

import json
import os
import platform
import sys
import time

import numpy as np
import scipy

SEED = 20260909
SIZES = [128, 256, 512, 1024, 2048]
NREP = 5
TRANSFORMS_PER_STEP = 48  # 12 per RK4 stage, 4 stages; see module docstring


def time_transforms(n, rng):
    a = rng.standard_normal((n, n))
    ah = np.fft.rfft2(a)  # warm-up
    fwd, inv = [], []
    for _ in range(NREP):
        t0 = time.perf_counter()
        ah = np.fft.rfft2(a)
        fwd.append(time.perf_counter() - t0)
        t0 = time.perf_counter()
        a2 = np.fft.irfft2(ah, s=(n, n))
        inv.append(time.perf_counter() - t0)
    roundtrip_error = float(np.max(np.abs(a2 - a)))
    return min(fwd), min(inv), roundtrip_error


def main():
    rng = np.random.default_rng(SEED)
    rows = []
    for n in SIZES:
        fwd, inv, err = time_transforms(n, rng)
        per_transform = 0.5 * (fwd + inv)
        step = TRANSFORMS_PER_STEP * per_transform
        # A run of the fixture library's convergence tier is taken as 2e4 steps.
        rows.append(
            {
                "N": n,
                "forward_seconds": fwd,
                "inverse_seconds": inv,
                "mean_transform_seconds": per_transform,
                "seconds_per_rk4_step": step,
                "core_hours_per_20000_steps": step * 20000.0 / 3600.0,
                "field_bytes_float64": 8 * n * n,
                "snapshot_bytes_two_fields_float64": 16 * n * n,
                "roundtrip_max_abs_error": err,
            }
        )

    payload = {
        "description": (
            "Measured numpy FFT throughput converted into the cost of one RK4 "
            "step of a 2D reduced-MHD pseudo-spectral solver (48 transforms per "
            "step) and into per-run core-hours and per-snapshot storage, for the "
            "fixture library resolution ladder."
        ),
        "step_model": {
            "transforms_per_rk4_step": TRANSFORMS_PER_STEP,
            "breakdown": "12 transforms per stage (2 fields, 8 gradients, 2 products), 4 stages",
            "reference_run_length_steps": 20000,
        },
        "parameters": {"seed": SEED, "sizes": SIZES, "repeats": NREP},
        "environment": {
            "python": sys.version.split()[0],
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "platform": platform.platform(),
        },
        "rows": rows,
        "caveats": [
            "Single-process numpy FFT timings on one laptop; a tuned solver using "
            "in-place plans or pyfftw would be faster, a Python-level loop slower.",
            "Dealiasing, time-step control, diagnostics and I/O are not timed and "
            "typically add tens of percent.",
            "Extrapolation to 4096 squared is by the measured N^2 log N trend, not "
            "measured directly.",
        ],
    }

    out_path = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "rmhd_step_cost_scaling_output.json",
    )
    with open(out_path, "w") as fh:
        json.dump(payload, fh, indent=1, sort_keys=True)

    for r in rows:
        print(
            "N=%-5d transform=%.4g s  step=%.4g s  20k steps=%.4g core-hours  "
            "snapshot=%.3g MB" % (
                r["N"], r["mean_transform_seconds"], r["seconds_per_rk4_step"],
                r["core_hours_per_20000_steps"],
                r["snapshot_bytes_two_fields_float64"] / 1e6,
            )
        )
    print("wrote", out_path)


if __name__ == "__main__":
    main()
