"""
CGLS semi-convergence, stopping rules, and float32 vs float64 arithmetic
on two ill-posed 1D Fredholm-first-kind fixtures (shaw, deconv_gauss).

Common C, role: iterative-and-sparse-investigator (Part 1).

Usage:
    python inv_iterative.py            # writes inv_iterative_output.json next to this file

Fixtures: if inv_fixtures.py (from the fixture-builder role) is importable and
exposes shaw()/deconv_gauss(), it is used.  Otherwise the minimal fallback
generators defined in this file are used, and the output JSON records which
path was taken in the field "fixture_source".

All randomness is seeded explicitly (numpy default_rng with integer seeds).
"""

import hashlib
import json
import os
import platform
import sys
import time

import numpy as np
import scipy
from scipy.sparse.linalg import lsqr as scipy_lsqr

HERE = os.path.dirname(os.path.abspath(__file__))

# ----------------------------------------------------------------------------
# Fixtures
# ----------------------------------------------------------------------------

FIXTURE_SOURCE = "builtin_fallback_in_inv_iterative.py"
_ext = None
try:
    sys.path.insert(0, HERE)
    import inv_fixtures as _ext  # type: ignore
    if hasattr(_ext, "shaw") and hasattr(_ext, "deconv_gauss"):
        FIXTURE_SOURCE = "inv_fixtures.py"
    else:
        _ext = None
except Exception:
    _ext = None


def shaw_local(n):
    """Standard 'shaw' 1D image restoration test problem (Fredholm 1st kind).

    Discretization as described for Hansen's Regularization Tools shaw(n)
    (midpoint rule on [-pi/2, pi/2]).  NOT verified to be byte-identical to the
    MATLAB implementation; used here only as a fixed, reproducible ill-posed
    operator with exponentially decaying singular values.
    """
    h = np.pi / n
    idx = np.arange(1, n + 1)
    s = (idx - n / 2.0 - 0.5) * h          # co-latitudes
    A = np.zeros((n, n))
    S, T = np.meshgrid(s, s, indexing="ij")
    u = np.pi * (np.sin(S) + np.sin(T))
    sinc = np.ones_like(u)
    nz = np.abs(u) > 1e-14
    sinc[nz] = np.sin(u[nz]) / u[nz]
    A = (np.cos(S) + np.cos(T)) ** 2 * sinc ** 2 * h
    a1, c1, t1 = 2.0, 6.0, 0.8
    a2, c2, t2 = 1.0, 2.0, -0.5
    x = a1 * np.exp(-c1 * (s - t1) ** 2) + a2 * np.exp(-c2 * (s - t2) ** 2)
    return A, x, A @ x


def deconv_gauss_local(n, sigma=0.05):
    """1D Gaussian-kernel deconvolution on [0,1] with a piecewise-smooth signal.

    A[i,j] = h * exp(-(s_i - t_j)^2 / (2 sigma^2)) / (sqrt(2 pi) sigma).
    x_true = box + smooth bump + linear ramp segment (piecewise smooth).
    """
    h = 1.0 / n
    t = (np.arange(n) + 0.5) * h
    S, T = np.meshgrid(t, t, indexing="ij")
    A = h * np.exp(-((S - T) ** 2) / (2 * sigma ** 2)) / (np.sqrt(2 * np.pi) * sigma)
    x = np.zeros(n)
    x[(t >= 0.10) & (t < 0.25)] = 1.0                       # box
    m = (t >= 0.40) & (t < 0.60)
    x[m] = 1.5 * np.exp(-((t[m] - 0.50) ** 2) / (2 * 0.03 ** 2))  # bump
    m = (t >= 0.72) & (t < 0.92)
    x[m] = 2.0 * (t[m] - 0.72) / 0.20                        # ramp
    return A, x, A @ x


USED_SOURCE = {"shaw": None, "deconv_gauss": None}


def _unpack(out):
    """inv_fixtures.Fixture is a dict subclass with keys A / x_true / b_clean;
    tolerate a plain (A, x_true, b_clean) tuple as well."""
    if isinstance(out, dict):
        return out["A"], out["x_true"], out["b_clean"]
    return out[0], out[1], out[2]


def get_fixture(name, n):
    if _ext is not None:
        try:
            out = _ext.shaw(n) if name == "shaw" else _ext.deconv_gauss(n, 0.05)
            A, x, b = _unpack(out)
            A = np.asarray(A, float)
            if A.shape == (n, n):
                USED_SOURCE[name] = "inv_fixtures.py"
                return A, np.asarray(x, float), np.asarray(b, float)
        except Exception as exc:  # pragma: no cover
            print("inv_fixtures fallback for", name, repr(exc))
    USED_SOURCE[name] = "builtin_fallback_in_inv_iterative.py"
    if name == "shaw":
        return shaw_local(n)
    return deconv_gauss_local(n, 0.05)


# ----------------------------------------------------------------------------
# CGLS (Hestenes-Stiefel CG on the normal equations), dtype-parametric
# ----------------------------------------------------------------------------

def cgls_iterates(A, b, maxit, dtype=np.float64):
    """Return list of iterates x_0..x_maxit and residual norms ||b - A x_k||.

    All arithmetic is performed in `dtype`.
    """
    A = np.asarray(A, dtype)
    b = np.asarray(b, dtype)
    n = A.shape[1]
    x = np.zeros(n, dtype)
    r = b.copy()
    s = A.T @ r
    p = s.copy()
    gamma = float(s @ s)
    xs = [x.copy()]
    res = [float(np.linalg.norm(r))]
    for _ in range(maxit):
        q = A @ p
        qq = float(q @ q)
        if qq <= 0.0 or gamma <= 0.0 or not np.isfinite(qq):
            xs.append(x.copy())
            res.append(float(np.linalg.norm(r)))
            continue
        alpha = dtype(gamma / qq)
        x = x + alpha * p
        r = r - alpha * q
        s = A.T @ r
        gamma_new = float(s @ s)
        beta = dtype(gamma_new / gamma) if gamma > 0 else dtype(0.0)
        p = s + beta * p
        gamma = gamma_new
        xs.append(x.copy())
        res.append(float(np.linalg.norm(r)))
    return xs, np.array(res)


def relerr(x, xt):
    return float(np.linalg.norm(np.asarray(x, float) - xt) / np.linalg.norm(xt))


def first_index_below(res, thresh):
    idx = np.nonzero(res <= thresh)[0]
    return int(idx[0]) if idx.size else None


# ----------------------------------------------------------------------------
# Experiment 1: semi-convergence and stopping rules
# ----------------------------------------------------------------------------

MAXIT = 200
FIXED_KS = [10, 25, 50, 100]
TAU = 1.01
SEEDS = list(range(20))


def run_semiconvergence(fixture, n, noise_level):
    A, x_true, b_clean = get_fixture(fixture, n)
    nb = np.linalg.norm(b_clean)
    rec = {
        "fixture": fixture, "n": n, "noise_level": noise_level,
        "cond2": float(np.linalg.cond(A)),
        "seeds": SEEDS, "maxit": MAXIT, "tau": TAU, "fixed_ks": FIXED_KS,
        "per_seed": [],
    }
    for seed in SEEDS:
        rng = np.random.default_rng(1000 + seed)
        e = rng.standard_normal(n)
        e *= noise_level * nb / np.linalg.norm(e)
        b = b_clean + e
        delta = float(np.linalg.norm(e))
        xs, res = cgls_iterates(A, b, MAXIT, np.float64)
        errs = np.array([relerr(x, x_true) for x in xs])
        kopt = int(np.argmin(errs))
        kdp = first_index_below(res, TAU * delta)
        row = {
            "seed": seed, "delta": delta,
            "k_opt": kopt, "err_opt": float(errs[kopt]),
            "k_dp": kdp,
            "err_dp": float(errs[kdp]) if kdp is not None else None,
            "err_at_fixed": {str(k): float(errs[k]) for k in FIXED_KS},
            "err_last": float(errs[MAXIT]),
            "err_curve_subsample": [float(errs[i]) for i in range(0, MAXIT + 1, 5)],
        }
        rec["per_seed"].append(row)
    return rec


def summarize(rec):
    ps = rec["per_seed"]
    kopt = np.array([r["k_opt"] for r in ps], float)
    eopt = np.array([r["err_opt"] for r in ps], float)
    kdp = np.array([r["k_dp"] if r["k_dp"] is not None else np.nan for r in ps], float)
    edp = np.array([r["err_dp"] if r["err_dp"] is not None else np.nan for r in ps], float)
    elast = np.array([r["err_last"] for r in ps], float)
    out = {
        "k_opt_median": float(np.median(kopt)),
        "k_opt_min": float(np.min(kopt)), "k_opt_max": float(np.max(kopt)),
        "err_opt_median": float(np.median(eopt)),
        "k_dp_median": float(np.nanmedian(kdp)),
        "dp_never_triggered_count": int(np.sum(np.isnan(kdp))),
        "err_dp_median": float(np.nanmedian(edp)),
        "dp_over_oracle_median": float(np.nanmedian(edp / eopt)),
        "dp_failure_rate_2x": float(np.mean(np.nan_to_num(edp / eopt, nan=np.inf) > 2.0)),
        "err_last_median": float(np.median(elast)),
        "last_over_oracle_median": float(np.median(elast / eopt)),
        "fixed_k": {},
    }
    for k in rec["fixed_ks"]:
        ek = np.array([r["err_at_fixed"][str(k)] for r in ps], float)
        out["fixed_k"][str(k)] = {
            "err_median": float(np.median(ek)),
            "over_oracle_median": float(np.median(ek / eopt)),
            "failure_rate_2x": float(np.mean(ek / eopt > 2.0)),
        }
    return out


# ----------------------------------------------------------------------------
# Experiment 2: float32 vs float64 on shaw
# ----------------------------------------------------------------------------

def run_precision(fixture="shaw", ns=(64, 128, 256), noise_level=0.01, seeds=range(20)):
    out = []
    for n in ns:
        A, x_true, b_clean = get_fixture(fixture, n)
        sv = np.linalg.svd(A, compute_uv=False)
        nb = np.linalg.norm(b_clean)
        rows = []
        for seed in seeds:
            rng = np.random.default_rng(2000 + seed)
            e = rng.standard_normal(n)
            e *= noise_level * nb / np.linalg.norm(e)
            b = b_clean + e
            r = {"seed": int(seed)}
            for tag, dt in (("f64", np.float64), ("f32", np.float32)):
                xs, res = cgls_iterates(A, b, MAXIT, dt)
                errs = np.array([relerr(x, x_true) for x in xs])
                errs = np.where(np.isfinite(errs), errs, np.inf)
                k = int(np.argmin(errs))
                r["k_opt_" + tag] = k
                r["err_opt_" + tag] = float(errs[k])
                r["err_last_" + tag] = float(errs[MAXIT]) if np.isfinite(errs[MAXIT]) else None
                r["err_at50_" + tag] = float(errs[50]) if np.isfinite(errs[50]) else None
            rows.append(r)
        k64 = np.array([r["k_opt_f64"] for r in rows], float)
        k32 = np.array([r["k_opt_f32"] for r in rows], float)
        e64 = np.array([r["err_opt_f64"] for r in rows], float)
        e32 = np.array([r["err_opt_f32"] for r in rows], float)
        l64 = np.array([r["err_last_f64"] for r in rows], float)
        l32 = np.array([r["err_last_f32"] for r in rows], float)
        out.append({
            "n": n, "noise_level": noise_level,
            "cond2": float(sv[0] / sv[-1]),
            "smallest_sv": float(sv[-1]), "largest_sv": float(sv[0]),
            "sv_below_f32_eps_rel": int(np.sum(sv / sv[0] < np.finfo(np.float32).eps)),
            "k_opt_median_f64": float(np.median(k64)),
            "k_opt_median_f32": float(np.median(k32)),
            "k_opt_max_abs_diff": float(np.max(np.abs(k64 - k32))),
            "k_opt_identical_fraction": float(np.mean(k64 == k32)),
            "err_opt_median_f64": float(np.median(e64)),
            "err_opt_median_f32": float(np.median(e32)),
            "err_opt_median_rel_diff": float(np.median(np.abs(e32 - e64) / e64)),
            "err_last_median_f64": float(np.median(l64)),
            "err_last_median_f32": float(np.median(l32)),
            "per_seed": rows,
        })
    return out


# ----------------------------------------------------------------------------
# Cross-check: CGLS iterate vs scipy.sparse.linalg.lsqr at the same iteration
# ----------------------------------------------------------------------------

def lsqr_crosscheck(fixture="shaw", n=128, noise_level=0.01, ks=(5, 10, 20, 40)):
    A, x_true, b_clean = get_fixture(fixture, n)
    rng = np.random.default_rng(31415)
    e = rng.standard_normal(n)
    e *= noise_level * np.linalg.norm(b_clean) / np.linalg.norm(e)
    b = b_clean + e
    xs, _ = cgls_iterates(A, b, max(ks), np.float64)
    out = []
    for k in ks:
        xl = scipy_lsqr(A, b, atol=0.0, btol=0.0, conlim=0.0, iter_lim=k)[0]
        d = np.linalg.norm(xl - xs[k]) / max(np.linalg.norm(xs[k]), 1e-300)
        out.append({"k": int(k), "rel_diff_cgls_vs_lsqr": float(d),
                    "relerr_cgls": relerr(xs[k], x_true),
                    "relerr_lsqr": relerr(xl, x_true)})
    return out


def _sha_of(fname):
    p = os.path.join(HERE, fname)
    if not os.path.exists(p):
        return None
    with open(p, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def main():
    t0 = time.time()
    results = {
        "script": "inv_iterative.py",
        "role": "iterative-and-sparse-investigator",
        "common": "C (failure maps for numerical methods on ill-conditioned inverse problems)",
        "fixture_source": FIXTURE_SOURCE,
        "fixture_source_used_per_fixture": USED_SOURCE,
        "inv_fixtures_sha256": _sha_of("inv_fixtures.py"),
        "env": {
            "python": sys.version.split()[0],
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "platform": platform.platform(),
        },
        "semiconvergence": [],
        "precision": None,
        "lsqr_crosscheck": None,
    }
    for fixture, n in (("shaw", 128), ("deconv_gauss", 128)):
        for lvl in (0.01, 0.05):
            rec = run_semiconvergence(fixture, n, lvl)
            rec["summary"] = summarize(rec)
            results["semiconvergence"].append(rec)
            s = rec["summary"]
            print(f"{fixture} n={n} noise={lvl}: k_opt_med={s['k_opt_median']} "
                  f"err_opt_med={s['err_opt_median']:.4f} k_dp_med={s['k_dp_median']} "
                  f"dp/oracle={s['dp_over_oracle_median']:.3f} "
                  f"last/oracle={s['last_over_oracle_median']:.3g}")
    results["precision"] = run_precision()
    for p in results["precision"]:
        print(f"precision shaw n={p['n']}: cond={p['cond2']:.3e} "
              f"k_opt f64={p['k_opt_median_f64']} f32={p['k_opt_median_f32']} "
              f"identical={p['k_opt_identical_fraction']:.2f} "
              f"err_opt f64={p['err_opt_median_f64']:.4f} f32={p['err_opt_median_f32']:.4f} "
              f"err_last f64={p['err_last_median_f64']:.3g} f32={p['err_last_median_f32']:.3g}")
    results["lsqr_crosscheck"] = lsqr_crosscheck()
    print("lsqr crosscheck:", results["lsqr_crosscheck"])
    results["runtime_seconds"] = round(time.time() - t0, 2)

    path = os.path.join(HERE, "inv_iterative_output.json")
    with open(path, "w") as f:
        json.dump(results, f, indent=1, sort_keys=True)
    with open(path, "rb") as f:
        h = hashlib.sha256(f.read()).hexdigest()
    print("wrote", path, "sha256", h, "in", results["runtime_seconds"], "s")


if __name__ == "__main__":
    main()
