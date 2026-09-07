"""
inv_paramchoice.py -- failure map for Tikhonov regularization parameter-choice
rules on ill-conditioned linear inverse problems.

Rules compared: L-curve corner (maximum curvature of the log-log L-curve, in the
sense of Hansen & O'Leary 1993), generalized cross-validation (Golub, Heath &
Wahba 1979), and the Morozov discrepancy principle with the TRUE noise norm and
with a mis-specified noise level (x0.5 and x2.0).

Baseline: the oracle Tikhonov parameter, i.e. the grid point minimizing the true
relative error ||x_lam - x_true|| / ||x_true||.  A rule "fails" on a given
instance when its relative error exceeds 2x the oracle relative error.

Noise: white Gaussian and AR(1)-correlated (rho = 0.5, 0.9), both scaled to a
prescribed relative noise level ||e|| / ||b_clean||.  20 seeds per cell.

Fixtures: it first tries to import `inv_fixtures` (written by the fixture-builder
role in the same Common).  If that module is unavailable, it falls back to the
self-contained `shaw` and `deconv_gauss` generators defined below; the output
JSON records which source was used in the field "fixture_source".

Usage:
    python inv_paramchoice.py [output.json]

Deterministic: every random draw is seeded (numpy default_rng with an explicit
integer derived from the cell and the seed index).
"""

import hashlib
import json
import os
import platform
import sys

import numpy as np
import scipy

# --------------------------------------------------------------------------
# Fixtures
# --------------------------------------------------------------------------

FIXTURE_SOURCE = None
_fx = None
try:  # prefer the shared fixture library if the fixture-builder has landed it
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import inv_fixtures as _fx  # type: ignore
    if hasattr(_fx, "shaw") and hasattr(_fx, "deconv_gauss"):
        FIXTURE_SOURCE = "inv_fixtures.py (shared fixture library)"
    else:
        _fx = None
except Exception:
    _fx = None
if _fx is None:
    FIXTURE_SOURCE = "local generators inside inv_paramchoice.py (inv_fixtures.py absent)"


def shaw_local(n):
    """Shaw's 1-D image-restoration test problem (standard discretization).

    Kernel K(s,t) = (cos s + cos t)^2 * (sin(u)/u)^2, u = pi (sin s + sin t),
    on [-pi/2, pi/2], midpoint collocation with n points.  This is the standard
    test problem popularized by Hansen's Regularization Tools; the code here is
    an independent NumPy implementation and has NOT been checked for numerical
    identity against Hansen's MATLAB shaw.m.
    """
    h = np.pi / n
    s = -np.pi / 2 + (np.arange(n) + 0.5) * h
    co = np.cos(s)
    u = np.pi * (np.sin(s)[:, None] + np.sin(s)[None, :])
    A = ((co[:, None] + co[None, :]) * np.sinc(u / np.pi)) ** 2 * h
    x = 2.0 * np.exp(-6.0 * (s - 0.8) ** 2) + 1.0 * np.exp(-2.0 * (s + 0.5) ** 2)
    return A, x, A @ x


def deconv_gauss_local(n, sigma=0.05):
    """1-D Gaussian-kernel deconvolution with a piecewise-smooth true signal.

    A_ij = h * exp(-(t_i - t_j)^2 / (2 sigma^2)) / (sqrt(2 pi) sigma), t in [0,1].
    x_true = smooth bump + rectangular plateau + a linear ramp segment.
    """
    h = 1.0 / n
    t = (np.arange(n) + 0.5) * h
    D = t[:, None] - t[None, :]
    A = h * np.exp(-(D ** 2) / (2 * sigma ** 2)) / (np.sqrt(2 * np.pi) * sigma)
    x = np.zeros(n)
    x += 1.0 * np.exp(-((t - 0.2) ** 2) / (2 * 0.03 ** 2))
    x[(t > 0.45) & (t < 0.62)] += 0.8
    ramp = (t > 0.72) & (t < 0.92)
    x[ramp] += 1.2 * (t[ramp] - 0.72) / 0.20
    return A, x, A @ x


def get_fixture(name):
    """name in {shaw, deconv_gauss, deconv_smooth}.

    deconv_smooth reuses the deconv_gauss operator but with a very smooth true
    solution (a single wide Gaussian), the classic setting in which the L-curve
    is reported to lose its distinct corner.
    """
    if name == "deconv_smooth":
        A, _, _ = get_fixture("deconv_gauss")
        n = A.shape[1]
        t = (np.arange(n) + 0.5) / n
        x = np.exp(-((t - 0.5) ** 2) / (2 * 0.15 ** 2))
        return A, x, A @ x
    if _fx is not None:
        if name == "shaw":
            out = _fx.shaw(64)
        else:
            out = _fx.deconv_gauss(128, 0.05)
        # tolerate a mapping with keys A/x_true/b_clean, or a tuple in that order
        if hasattr(out, "keys"):
            A, x_true, b_clean = out["A"], out["x_true"], out["b_clean"]
        else:
            A, x_true, b_clean = out[0], out[1], out[2]
        return np.asarray(A, float), np.asarray(x_true, float), np.asarray(b_clean, float)
    if name == "shaw":
        return shaw_local(64)
    return deconv_gauss_local(128, 0.05)


# --------------------------------------------------------------------------
# Noise
# --------------------------------------------------------------------------

def make_noise(m, kind, rho, rel_level, b_clean, rng):
    """Noise vector scaled so that ||e|| = rel_level * ||b_clean||."""
    if kind == "white":
        e = rng.standard_normal(m)
    elif kind == "ar1":
        w = rng.standard_normal(m)
        e = np.empty(m)
        # stationary start so the process is not transient at i = 0
        e[0] = w[0] / np.sqrt(1.0 - rho ** 2)
        for i in range(1, m):
            e[i] = rho * e[i - 1] + w[i]
    else:
        raise ValueError(kind)
    e = e / np.linalg.norm(e) * rel_level * np.linalg.norm(b_clean)
    return e


# --------------------------------------------------------------------------
# Tikhonov via SVD and the parameter-choice rules
# --------------------------------------------------------------------------

def tikh_curves(U, s, Vt, b, lambdas):
    """Return solutions' norms, residual norms, errors-ready solutions.

    x_lam = V diag(s / (s^2 + lam^2)) U^T b  (square/overdetermined, full rank
    assumed after truncation of exact zeros).
    """
    beta = U.T @ b                       # (k,)
    k = s.size
    f = (s ** 2)[None, :] / ((s ** 2)[None, :] + (lambdas ** 2)[:, None])  # (L,k)
    coef = f * (beta / s)[None, :]
    X = coef @ Vt                         # (L, n)
    res_par = ((1.0 - f) * beta[None, :]) ** 2
    # residual component outside the range of the retained singular vectors
    b_out2 = max(float(b @ b - beta @ beta), 0.0)
    resid = np.sqrt(res_par.sum(axis=1) + b_out2)
    xnorm = np.linalg.norm(X, axis=1)
    trace = float(b.size) - f.sum(axis=1)  # m - sum(filter factors)
    return X, resid, xnorm, trace, k


def pick_gcv(resid, trace, m):
    G = m * resid ** 2 / np.maximum(trace, 1e-300) ** 2
    return int(np.argmin(G)), G


def pick_lcurve(resid, xnorm, lambdas):
    """Maximum-curvature corner of the log-log L-curve.

    Curvature of the parametric curve (xi(t), eta(t)) with t = log(lambda),
    xi = log(resid), eta = log(xnorm), by centred finite differences on the
    (uniform in log lambda) grid.  Boundary points are excluded so an edge of
    the grid cannot be reported as a corner.
    """
    xi = np.log(np.maximum(resid, 1e-300))
    eta = np.log(np.maximum(xnorm, 1e-300))
    t = np.log(lambdas)
    d1x = np.gradient(xi, t)
    d1y = np.gradient(eta, t)
    d2x = np.gradient(d1x, t)
    d2y = np.gradient(d1y, t)
    denom = (d1x ** 2 + d1y ** 2) ** 1.5
    kappa = (d1x * d2y - d2x * d1y) / np.maximum(denom, 1e-300)
    interior = np.zeros_like(kappa, dtype=bool)
    interior[2:-2] = True
    kk = np.where(interior, kappa, -np.inf)
    idx = int(np.argmax(kk))
    # diagnostics: how many interior local maxima of the curvature exist, and
    # how peaked the chosen corner is relative to the curvature spread.
    loc = 0
    for i in range(3, kappa.size - 3):
        if kappa[i] > kappa[i - 1] and kappa[i] >= kappa[i + 1] and kappa[i] > 0:
            loc += 1
    peak = float(kappa[idx])
    diag = {"n_local_curvature_maxima": loc, "corner_curvature": peak}
    return idx, kappa, diag


def pick_discrepancy(resid, lambdas, delta_assumed, tau=1.0):
    """Largest lambda on the grid with ||A x_lam - b|| <= tau * delta_assumed.

    resid is increasing in lambda, so this is the standard Morozov choice.
    Returns (index, status) with status in {ok, all_below, infeasible_target}.
    "infeasible_target" means no lambda on the grid reaches the requested
    discrepancy: the residual has a floor (components of b outside the retained
    numerical range of A cannot be fitted at any lambda).  The fallback is then
    the smallest lambda on the grid, which is what an implementation that simply
    drives lambda down until the discrepancy is met would produce.
    """
    target = tau * delta_assumed
    ok = np.nonzero(resid <= target)[0]
    if ok.size == 0:
        return 0, "infeasible_target"
    if ok.size == resid.size:
        return int(resid.size - 1), "all_below"
    return int(ok[-1]), "ok"


# --------------------------------------------------------------------------
# Experiment
# --------------------------------------------------------------------------

METHODS = ["lcurve", "gcv", "discrepancy_true", "discrepancy_half", "discrepancy_double"]
FAIL_FACTOR = 2.0
N_SEEDS = 20
N_LAMBDA = 300


def run_cell(fixture, noise_kind, rho, rel_level, n_seeds=N_SEEDS):
    A, x_true, b_clean = get_fixture(fixture)
    m, n = A.shape
    U, s, Vt = np.linalg.svd(A, full_matrices=False)
    keep = s > s[0] * 1e-14
    U, s, Vt = U[:, keep], s[keep], Vt[keep, :]
    lambdas = np.logspace(np.log10(s[0]), np.log10(max(s[-1], s[0] * 1e-12)) - 2.0, N_LAMBDA)
    lambdas = np.sort(lambdas)
    xt_norm = np.linalg.norm(x_true)

    per_seed = []
    for si in range(n_seeds):
        key = "%s|%s|%.3f|%.4f|%d" % (fixture, noise_kind, rho, rel_level, si)
        seed = int(hashlib.sha256(key.encode()).hexdigest()[:8], 16)
        rng = np.random.default_rng(seed)
        e = make_noise(m, noise_kind, rho, rel_level, b_clean, rng)
        b = b_clean + e
        delta = float(np.linalg.norm(e))
        X, resid, xnorm, trace, k = tikh_curves(U, s, Vt, b, lambdas)
        errs = np.linalg.norm(X - x_true[None, :], axis=1) / xt_norm
        i_or = int(np.argmin(errs))
        e_or = float(errs[i_or])

        picks = {}
        i_l, _, ldiag = pick_lcurve(resid, xnorm, lambdas)
        picks["lcurve"] = (i_l, "ok")
        i_g, _ = pick_gcv(resid, trace, m)
        picks["gcv"] = (i_g, "ok")
        for tag, fac in (("discrepancy_true", 1.0), ("discrepancy_half", 0.5),
                         ("discrepancy_double", 2.0)):
            picks[tag] = pick_discrepancy(resid, lambdas, fac * delta)

        rec = {"seed": seed, "oracle_relerr": e_or,
               "oracle_lambda": float(lambdas[i_or]), "noise_norm": delta,
               "resid_floor_over_delta": float(resid.min() / delta),
               "lcurve_diagnostics": ldiag, "methods": {}}
        for meth, (idx, status) in picks.items():
            rec["methods"][meth] = {
                "lambda": float(lambdas[idx]),
                "relerr": float(errs[idx]),
                "ratio_to_oracle": float(errs[idx] / e_or),
                "log10_lambda_offset": float(np.log10(lambdas[idx] / lambdas[i_or])),
                "failed": bool(errs[idx] > FAIL_FACTOR * e_or),
                "status": status,
            }
        per_seed.append(rec)

    summary = {}
    for meth in METHODS:
        ratios = np.array([r["methods"][meth]["ratio_to_oracle"] for r in per_seed])
        offs = np.array([r["methods"][meth]["log10_lambda_offset"] for r in per_seed])
        fails = np.array([r["methods"][meth]["failed"] for r in per_seed])
        summary[meth] = {
            "failure_rate": float(fails.mean()),
            "n_failures": int(fails.sum()),
            "median_ratio_to_oracle": float(np.median(ratios)),
            "max_ratio_to_oracle": float(ratios.max()),
            "median_relerr": float(np.median([r["methods"][meth]["relerr"] for r in per_seed])),
            "median_log10_lambda_offset": float(np.median(offs)),
            "status_counts": {st: int(sum(1 for r in per_seed
                                          if r["methods"][meth]["status"] == st))
                              for st in ("ok", "all_below", "infeasible_target")},
        }
        feas = [r for r in per_seed if r["methods"][meth]["status"] == "ok"]
        if len(feas) < len(per_seed):
            summary[meth]["n_feasible"] = len(feas)
            summary[meth]["failure_rate_feasible_only"] = (
                float(np.mean([r["methods"][meth]["failed"] for r in feas])) if feas else None)
            summary[meth]["median_ratio_feasible_only"] = (
                float(np.median([r["methods"][meth]["ratio_to_oracle"] for r in feas]))
                if feas else None)
    summary["lcurve"]["median_n_local_curvature_maxima"] = float(
        np.median([r["lcurve_diagnostics"]["n_local_curvature_maxima"] for r in per_seed]))
    summary["lcurve"]["median_corner_curvature"] = float(
        np.median([r["lcurve_diagnostics"]["corner_curvature"] for r in per_seed]))
    return {
        "fixture": fixture, "m": int(m), "n": int(n),
        "cond": float(s[0] / s[-1]),
        "numerical_rank_at_1e-14": int(s.size),
        "median_resid_floor_over_delta": float(
            np.median([r["resid_floor_over_delta"] for r in per_seed])),
        "noise": noise_kind, "rho": rho, "rel_noise_level": rel_level,
        "n_seeds": n_seeds, "n_lambda": N_LAMBDA,
        "lambda_min": float(lambdas[0]), "lambda_max": float(lambdas[-1]),
        "median_oracle_relerr": float(np.median([r["oracle_relerr"] for r in per_seed])),
        "summary": summary, "per_seed": per_seed,
    }


def main():
    out_path = sys.argv[1] if len(sys.argv) > 1 else "inv_paramchoice_output.json"
    cells = []
    for fixture in ("shaw", "deconv_gauss", "deconv_smooth"):
        for rel_level in (0.01, 0.05):
            for noise_kind, rho in (("white", 0.0), ("ar1", 0.5), ("ar1", 0.9)):
                cells.append(run_cell(fixture, noise_kind, rho, rel_level))

    # overall failure rate per method, and per noise class
    overall = {}
    for meth in METHODS:
        tot = sum(c["summary"][meth]["n_failures"] for c in cells)
        n = sum(c["n_seeds"] for c in cells)
        overall[meth] = {"failure_rate": tot / n, "n_failures": tot, "n_instances": n}
    by_noise = {}
    for c in cells:
        key = c["noise"] if c["noise"] == "white" else "ar1_rho%.1f" % c["rho"]
        d = by_noise.setdefault(key, {m: [0, 0] for m in METHODS})
        for meth in METHODS:
            d[meth][0] += c["summary"][meth]["n_failures"]
            d[meth][1] += c["n_seeds"]
    by_noise = {k: {m: {"failure_rate": v[m][0] / v[m][1], "n_failures": v[m][0],
                        "n_instances": v[m][1]} for m in METHODS}
                for k, v in by_noise.items()}

    by_fixture = {}
    for c in cells:
        d = by_fixture.setdefault(c["fixture"], {m: [0, 0] for m in METHODS})
        for meth in METHODS:
            d[meth][0] += c["summary"][meth]["n_failures"]
            d[meth][1] += c["n_seeds"]
    by_fixture = {k: {m: {"failure_rate": v[m][0] / v[m][1], "n_failures": v[m][0],
                          "n_instances": v[m][1]} for m in METHODS}
                  for k, v in by_fixture.items()}

    doc = {
        "script": "inv_paramchoice.py",
        "fixture_source": FIXTURE_SOURCE,
        "description": ("Tikhonov parameter-choice failure map: L-curve corner, GCV, and "
                        "Morozov discrepancy with true and mis-specified noise level, under "
                        "white and AR(1) noise. Failure = relative error > %.1fx the oracle "
                        "Tikhonov parameter's relative error." % FAIL_FACTOR),
        "fail_factor": FAIL_FACTOR,
        "methods": METHODS,
        "environment": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "platform": platform.platform(),
        },
        "overall_failure_rate": overall,
        "failure_rate_by_noise": by_noise,
        "failure_rate_by_fixture": by_fixture,
        "cells": cells,
    }
    txt = json.dumps(doc, indent=1, sort_keys=True)
    with open(out_path, "w") as fh:
        fh.write(txt)
    print("wrote", out_path, len(txt), "bytes")
    print("sha256", hashlib.sha256(txt.encode()).hexdigest())
    print("fixture_source:", FIXTURE_SOURCE)
    for meth in METHODS:
        print("%-22s overall failure rate %.3f" % (meth, overall[meth]["failure_rate"]))
    for key in sorted(by_noise):
        print("noise", key, {m: round(by_noise[key][m]["failure_rate"], 3) for m in METHODS})
    for key in sorted(by_fixture):
        print("fixture", key,
              {m: round(by_fixture[key][m]["failure_rate"], 3) for m in METHODS})
    for c in cells:
        print("%-13s %-5s rho=%.1f level=%.2f oracle=%.4f  " % (
            c["fixture"], c["noise"], c["rho"], c["rel_noise_level"],
            c["median_oracle_relerr"]),
            {m: (round(c["summary"][m]["failure_rate"], 2),
                 round(c["summary"][m]["median_ratio_to_oracle"], 2)) for m in METHODS},
            "Lcorners=%.1f" % c["summary"]["lcurve"]["median_n_local_curvature_maxima"])


if __name__ == "__main__":
    main()
