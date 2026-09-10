"""Stress cone spanned by families of exact Alfvenic disturbances.

Companion fixture to alfven_stress_cancellation.py. That script showed the net
wave stress (Reynolds minus Maxwell) of a single shear-Alfven wave vanishes at
equipartition and degrades linearly in the equipartition-breaking parameter eps.
This script asks the next question, which is the one the 2026 forced
Navier-Stokes construction needs answered in the affirmative: taking several
such wave families at once, with NON-NEGATIVE amplitude weights (the weights are
squares of pulse amplitudes in that construction, so they cannot change sign),
what set of mean stress tensors can they supply? That set is the "admissible
stress cone".

Setup. Uniform mean field B0 along z, uniform density, Alfven units so that b is
the magnetic fluctuation divided by sqrt(mu0 rho). A shear-Alfven disturbance is
transverse and depends on z alone:
    v = (vx(z), vy(z), 0),  b = -(1 + eps) * v,  eps = 0 the exact equipartition.
Sign convention (same as the companion script): the wave stress is
    T_ij = <v_i v_j> - <b_i b_j>,   Reynolds part minus Maxwell part,
angle brackets = average over one periodic box in z.

Two structural facts this script measures rather than assumes.

(1) Parallel-flux components vanish identically. v_z = b_z = 0 for the whole
    class, so T_xz = T_yz = T_zz = 0 for every family, every amplitude and every
    eps. The component that would transport transverse momentum ALONG the mean
    field - the structural analogue of the radial momentum flux the Navier-Stokes
    construction needs - is not merely small, it is exactly zero.

(2) In the transverse block the net stress of family s is
        T^(s) = -(2 eps_s + eps_s^2) R^(s),   R^(s)_ij = <v_i v_j>,
    a scalar multiple of that family's own (positive semi-definite) Reynolds
    tensor. So a family supplies a single RAY in stress space whose direction is
    fixed by its polarization and whose sign is fixed by the sign of eps_s, and
    whose magnitude is O(eps) relative to the fluctuation energy it carries.

Cone tests performed.
 A. Single-sign set (all eps_s > 0): is the conical hull of the generators
    pointed, i.e. contained in some closed half-space? Tested by asking a linear
    program for a separating vector c with c . g_s >= 1 for every generator g_s.
    Feasible => the cone is pointed => whole directions of stress space are
    unreachable with non-negative weights.
 B. Mixed-sign set (both signs of eps allowed, i.e. some families magnetically
    dominated and some kinetically dominated): same test.
 C. Reachability of a target with a parallel-flux component: exact distance from
    the linear span of the generators to the unit vector along T_xz.
 D. Magnitude scaling: Frobenius norm of the net stress divided by the total
    fluctuation energy, versus eps.

Run:
  <venv>/bin/python alfven_two_family_stress_cone.py
  (writes alfven_two_family_stress_cone_output.json next to this script)

Fixed seed 20260909 for the random phases. Runtime a few seconds.

Limitation, stated up front: this is the class of EXACT finite-amplitude Alfvenic
solutions on a uniform field. It is not a sheared-box amplification calculation,
so it does not test hypothesis H1 (stress-cone dimension) for amplified
disturbances in a sheared magnetized background. It bounds only what the exact
Alfvenic class can supply.
"""

import json
import os
import platform
import sys

import numpy as np
import scipy
from scipy.optimize import linprog

SEED = 20260909
NZ = 4096
LZ = 2.0 * np.pi
NMODES = 4

# (label, ax, ay, phase difference between the x and y components, eps)
FAMILIES_SINGLE_SIGN = [
    ("f1", 0.70, 0.90, 0.40, 0.10),
    ("f2", 0.90, 0.30, 1.90, 0.10),
    ("f3", 0.50, 0.50, 0.00, 0.30),
    ("f4", 0.20, 1.00, 2.60, 0.05),
    ("f5", 1.00, 0.60, 1.10, 0.20),
    ("f6", 0.80, 0.80, 2.20, 0.02),
]
FAMILIES_MIXED_SIGN = [
    (lab, ax, ay, dphi, eps if i % 2 == 0 else -eps)
    for i, (lab, ax, ay, dphi, eps) in enumerate(FAMILIES_SINGLE_SIGN)
]
EPS_SCAN = [1e-6, 1e-4, 1e-2, 1e-1, 0.5, 1.0]


def field(z, ax, ay, dphi, rng):
    """Transverse velocity of a finite-amplitude shear-Alfven packet."""
    vx = np.zeros_like(z)
    vy = np.zeros_like(z)
    for k in range(1, NMODES + 1):
        ph = rng.uniform(0.0, 2.0 * np.pi)
        vx += (ax / k) * np.cos(k * z + ph)
        vy += (ay / k) * np.cos(k * z + ph + dphi)
    return vx, vy


def stress_tensor(vx, vy, eps):
    """Full 3x3 net stress T_ij = <v_i v_j> - <b_i b_j>, b = -(1+eps) v."""
    vz = np.zeros_like(vx)
    v = [vx, vy, vz]
    b = [-(1.0 + eps) * c for c in v]
    T = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            T[i, j] = float(np.mean(v[i] * v[j]) - np.mean(b[i] * b[j]))
    return T


def energies(vx, vy, eps):
    ekin = 0.5 * float(np.mean(vx * vx + vy * vy))
    emag = 0.5 * (1.0 + eps) ** 2 * float(np.mean(vx * vx + vy * vy))
    return ekin, emag


def pointedness(generators, tol=1e-9):
    """Is the conical hull of `generators` contained in a closed half-space?

    Feasibility LP: find c with c . g_s >= 1 for all s (g_s normalized).
    Feasible  => pointed cone => some directions are unreachable.
    Infeasible => the generators positively span their linear span.
    """
    G = np.array([g / max(np.linalg.norm(g), 1e-300) for g in generators])
    n = G.shape[1]
    res = linprog(
        c=np.zeros(n),
        A_ub=-G,
        b_ub=-np.ones(G.shape[0]),
        bounds=[(-100.0, 100.0)] * n,
        method="highs",
    )
    out = {
        "separating_vector_found": bool(res.status == 0),
        "lp_status_message": str(res.message),
    }
    if res.status == 0:
        c = np.asarray(res.x, dtype=float)
        margins = (G @ c).tolist()
        out["separating_vector"] = c.tolist()
        out["min_margin"] = float(min(margins))
    return out


def main():
    z = np.linspace(0.0, LZ, NZ, endpoint=False)
    rng = np.random.default_rng(SEED)

    per_family = {}
    gens_single, gens_mixed = [], []
    max_abs_parallel = 0.0

    for (lab, ax, ay, dphi, eps), (_, _, _, _, eps_m) in zip(
        FAMILIES_SINGLE_SIGN, FAMILIES_MIXED_SIGN
    ):
        vx, vy = field(z, ax, ay, dphi, rng)
        T = stress_tensor(vx, vy, eps)
        Tm = stress_tensor(vx, vy, eps_m)
        R = np.array(
            [
                [float(np.mean(vx * vx)), float(np.mean(vx * vy))],
                [float(np.mean(vx * vy)), float(np.mean(vy * vy))],
            ]
        )
        pref = -(2.0 * eps + eps * eps)
        ekin, emag = energies(vx, vy, eps)
        resid = float(np.max(np.abs(T[:2, :2] - pref * R)))
        parallel = [float(T[0, 2]), float(T[1, 2]), float(T[2, 2])]
        max_abs_parallel = max(max_abs_parallel, max(abs(p) for p in parallel))
        per_family[lab] = {
            "params": {"ax": ax, "ay": ay, "dphi": dphi, "eps": eps,
                       "eps_mixed_set": eps_m},
            "reynolds_transverse": R.tolist(),
            "net_stress_transverse": T[:2, :2].tolist(),
            "net_stress_parallel_xz_yz_zz": parallel,
            "prefactor_minus_2eps_minus_eps2": pref,
            "max_abs_deviation_from_prefactor_times_R": resid,
            "energy_kinetic": ekin,
            "energy_magnetic": emag,
            "frobenius_net_over_total_energy": float(
                np.linalg.norm(T) / (ekin + emag)
            ),
        }
        gens_single.append(np.array([T[0, 0], T[1, 1], T[0, 1]]))
        gens_mixed.append(np.array([Tm[0, 0], Tm[1, 1], Tm[0, 1]]))

    def span_report(gens):
        G = np.array(gens)
        s = np.linalg.svd(G, compute_uv=False)
        rank = int(np.sum(s > 1e-12 * max(s[0], 1e-300)))
        return {
            "singular_values": s.tolist(),
            "linear_span_dimension": rank,
            "traces": [float(g[0] + g[1]) for g in gens],
            "cone": pointedness(gens),
        }

    # C. reachability of a target carrying a parallel momentum flux.
    # Generators live in the 3D transverse block of the 6D symmetric-tensor
    # space (xx, yy, zz, xy, xz, yz); the target is the unit vector along xz.
    gen6 = []
    for g in gens_mixed:
        gen6.append(np.array([g[0], g[1], 0.0, g[2], 0.0, 0.0]))
    G6 = np.array(gen6)
    target = np.array([0.0, 0.0, 0.0, 0.0, 1.0, 0.0])
    proj, *_ = np.linalg.lstsq(G6.T, target, rcond=None)
    residual_norm = float(np.linalg.norm(G6.T @ proj - target))

    # D. magnitude scaling with eps for one fixed family.
    rng2 = np.random.default_rng(SEED)
    vx, vy = field(z, 0.7, 0.9, 0.4, rng2)
    scan = []
    for eps in EPS_SCAN:
        T = stress_tensor(vx, vy, eps)
        ekin, emag = energies(vx, vy, eps)
        scan.append(
            {
                "epsilon": eps,
                "frobenius_net_stress": float(np.linalg.norm(T)),
                "total_fluctuation_energy": ekin + emag,
                "ratio": float(np.linalg.norm(T) / (ekin + emag)),
            }
        )

    payload = {
        "description": (
            "Admissible stress cone spanned by families of exact finite-amplitude "
            "shear-Alfven disturbances on a uniform field, with non-negative "
            "amplitude weights. Reports per-family net stress (Reynolds minus "
            "Maxwell), the identically vanishing parallel-flux components, the "
            "linear span dimension and pointedness of the cone for single-sign "
            "and mixed-sign equipartition breaking, the distance from the span "
            "to a target carrying a parallel momentum flux, and the O(eps) "
            "magnitude scaling."
        ),
        "sign_convention": "T_ij = <v_i v_j> - <b_i b_j>, Alfven units, b = -(1+eps) v",
        "parameters": {
            "seed": SEED,
            "nz": NZ,
            "Lz": LZ,
            "nmodes_per_family": NMODES,
            "families_single_sign": FAMILIES_SINGLE_SIGN,
            "families_mixed_sign": FAMILIES_MIXED_SIGN,
            "epsilon_scan": EPS_SCAN,
        },
        "environment": {
            "python": sys.version.split()[0],
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "platform": platform.platform(),
        },
        "per_family": per_family,
        "max_abs_parallel_flux_component_over_all_families": max_abs_parallel,
        "cone_single_sign_eps_positive": span_report(gens_single),
        "cone_mixed_sign_eps": span_report(gens_mixed),
        "parallel_flux_target": {
            "target": "unit symmetric tensor along xz, components (xx,yy,zz,xy,xz,yz)",
            "least_squares_residual_norm": residual_norm,
            "interpretation": (
                "residual 1.0 means the target is exactly orthogonal to the "
                "generators' linear span, so no combination of these families, "
                "with weights of either sign, supplies any parallel momentum flux"
            ),
        },
        "magnitude_scaling_with_epsilon": scan,
    }

    out_path = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "alfven_two_family_stress_cone_output.json",
    )
    with open(out_path, "w") as fh:
        json.dump(payload, fh, indent=1, sort_keys=True)

    print("max |parallel flux component| over all families:", max_abs_parallel)
    for name, key in [
        ("single-sign (all eps > 0)", "cone_single_sign_eps_positive"),
        ("mixed-sign eps", "cone_mixed_sign_eps"),
    ]:
        r = payload[key]
        print(name, "span dim =", r["linear_span_dimension"],
              "pointed =", r["cone"]["separating_vector_found"],
              "traces =", ["%+.4f" % t for t in r["traces"]])
    print("parallel-flux target lstsq residual:", residual_norm)
    for row in scan:
        print("  eps=%-8g |T|_F=%.6e  |T|/E=%.6e" %
              (row["epsilon"], row["frobenius_net_stress"], row["ratio"]))
    print("wrote", out_path)


if __name__ == "__main__":
    main()
