"""Exact Reynolds-minus-Maxwell stress cancellation for a shear-Alfven wave.

Fixture for task mhd-task-alfven-stress-cancellation (T3).

Setup. Uniform background field B0 along z, uniform density rho. Work in Alfven
units, so b denotes the magnetic fluctuation divided by sqrt(mu0 rho) and has the
dimensions of a velocity. A shear-Alfven disturbance is transverse and depends on
z only: v = (vx(z), vy(z), 0), b = -(1 + eps) * v, with eps = 0 the exact
equipartitioned Alfven wave (Elsasser field z+ = v + b = 0 identically).

Sign convention (fixed here, to be reconciled in T2): the wave stress is the
z-flux of perpendicular momentum,
    T_ij = <v_i v_j> - <b_i b_j>,
Reynolds part minus Maxwell part, angle brackets = average over one periodic box.

Analytic identity being checked. With b = -(1 + eps) v,
    T_ij = <v_i v_j> * (1 - (1 + eps)^2) = -<v_i v_j> * (2 eps + eps^2),
so T vanishes identically at eps = 0 for every component and every amplitude and
profile, and degrades linearly in eps with a known quadratic correction. The
ideal nonlinearity also vanishes identically for this class: with z+ and z- both
transverse and z-dependent only, (z+ . grad) z- = z+_z d_z z- = 0.

Run:
  <venv>/bin/python alfven_stress_cancellation.py > /dev/null
  (writes alfven_stress_cancellation_output.json next to this script)

No external input. One fixed seed (20260908) for the multi-mode random-phase
profile; single-mode cases are deterministic.
"""

import json
import os
import platform
import sys

import numpy as np
import scipy

SEED = 20260908
NZ = 4096
LZ = 2.0 * np.pi
AMPLITUDE = 0.7          # finite amplitude: the identity is amplitude-independent
POLARIZATION = 1.3       # vy/vx amplitude ratio
PHASE = 0.4              # phase offset, so <vx vy> is nonzero
EPS_LIST = [0.0, 1e-12, 1e-8, 1e-6, 1e-4, 1e-3, 1e-2, 1e-1, 0.5]


def single_mode(z, k=1):
    """Transverse velocity of one finite-amplitude shear-Alfven wave."""
    vx = AMPLITUDE * np.cos(k * z)
    vy = AMPLITUDE * POLARIZATION * np.cos(k * z + PHASE)
    return vx, vy


def multi_mode(z, rng, nmodes=6):
    """Random-phase superposition: still an exact solution of ideal MHD."""
    vx = np.zeros_like(z)
    vy = np.zeros_like(z)
    for k in range(1, nmodes + 1):
        a = AMPLITUDE / k
        vx += a * np.cos(k * z + rng.uniform(0.0, 2.0 * np.pi))
        vy += a * POLARIZATION * np.cos(k * z + rng.uniform(0.0, 2.0 * np.pi))
    return vx, vy


def stresses(vx, vy, eps):
    """Reynolds, Maxwell and net stress components for b = -(1+eps) v."""
    bx, by = -(1.0 + eps) * vx, -(1.0 + eps) * vy
    out = {}
    for name, (a1, a2, m1, m2) in {
        "xy": (vx, vy, bx, by),
        "xx": (vx, vx, bx, bx),
        "yy": (vy, vy, by, by),
    }.items():
        reynolds = float(np.mean(a1 * a2))
        maxwell = float(np.mean(m1 * m2))
        out[name] = {
            "reynolds": reynolds,
            "maxwell": maxwell,
            "net": reynolds - maxwell,
            "predicted_net": -reynolds * (2.0 * eps + eps * eps),
        }
    return out


def nonlinearity_residual(z, vx, vy, eps):
    """Spectral check that the ideal MHD nonlinearity vanishes for this class."""
    kz = np.fft.rfftfreq(z.size, d=(z[1] - z[0])) * 2.0 * np.pi

    def ddz(f):
        return np.fft.irfft(1j * kz * np.fft.rfft(f), n=z.size)

    bx, by = -(1.0 + eps) * vx, -(1.0 + eps) * vy
    zpx, zpy = vx + bx, vy + by
    zmx, zmy = vx - bx, vy - by
    # z-components of both Elsasser fields are identically zero, so the only
    # advective operator available, z^{+/-}_z d_z, annihilates everything.
    zp_z = np.zeros_like(z)
    zm_z = np.zeros_like(z)
    terms = [zp_z * ddz(zmx), zp_z * ddz(zmy), zm_z * ddz(zpx), zm_z * ddz(zpy)]
    return float(max(np.max(np.abs(t)) for t in terms))


def main():
    z = np.linspace(0.0, LZ, NZ, endpoint=False)
    rng = np.random.default_rng(SEED)
    fields = {"single_mode": single_mode(z), "multi_mode": multi_mode(z, rng)}

    results = {}
    for fname, (vx, vy) in fields.items():
        rows = []
        for eps in EPS_LIST:
            s = stresses(vx, vy, eps)
            row = {"epsilon": eps, "components": s}
            row["net_xy_over_epsilon"] = (
                s["xy"]["net"] / eps if eps != 0.0 else None
            )
            row["deviation_from_linear_xy"] = (
                s["xy"]["net"] + 2.0 * eps * s["xy"]["reynolds"]
            )
            row["net_minus_predicted_xy"] = s["xy"]["net"] - s["xy"]["predicted_net"]
            rows.append(row)
        results[fname] = {
            "reynolds_xy_at_eps0": float(np.mean(vx * vy)),
            "reynolds_xx_at_eps0": float(np.mean(vx * vx)),
            "rows": rows,
            "ideal_nonlinearity_max_abs_eps0": nonlinearity_residual(z, vx, vy, 0.0),
            "ideal_nonlinearity_max_abs_eps0p1": nonlinearity_residual(z, vx, vy, 0.1),
        }

    payload = {
        "description": (
            "Net wave stress (Reynolds minus Maxwell) for a finite-amplitude "
            "shear-Alfven wave on a uniform field, versus the equipartition-"
            "breaking parameter epsilon, with b = -(1+epsilon) v in Alfven units."
        ),
        "parameters": {
            "seed": SEED,
            "nz": NZ,
            "Lz": LZ,
            "amplitude": AMPLITUDE,
            "polarization_ratio": POLARIZATION,
            "phase_offset": PHASE,
            "epsilon_list": EPS_LIST,
        },
        "environment": {
            "python": sys.version.split()[0],
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "platform": platform.platform(),
        },
        "results": results,
    }

    out_path = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "alfven_stress_cancellation_output.json",
    )
    with open(out_path, "w") as fh:
        json.dump(payload, fh, indent=1, sort_keys=True)

    for fname, res in results.items():
        print(fname, "R_xy =", res["reynolds_xy_at_eps0"])
        for row in res["rows"]:
            print(
                "  eps=%-8g net_xy=%+.6e ratio=%s dev_lin=%+.3e" % (
                    row["epsilon"],
                    row["components"]["xy"]["net"],
                    ("%+.6f" % row["net_xy_over_epsilon"])
                    if row["net_xy_over_epsilon"] is not None else "n/a",
                    row["deviation_from_linear_xy"],
                )
            )
        print("  ideal nonlinearity max |.| :", res["ideal_nonlinearity_max_abs_eps0"])
    print("wrote", out_path)


if __name__ == "__main__":
    main()
