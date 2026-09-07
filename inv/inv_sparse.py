"""
ISTA / FISTA support recovery versus mutual coherence on an overcomplete
dictionary fixture.

Common C, role: iterative-and-sparse-investigator (Part 2).

Usage:
    python inv_sparse.py           # writes inv_sparse_output.json next to this file

Fixtures: if inv_fixtures.py (fixture-builder role) exposes
coherent_dictionary(m, n, coherence_target, seed) it is used; otherwise the
fallback generator in this file is used.  The output JSON records which in
"fixture_source".  Two dictionary constructions are exercised:
  "common"  - every column shares a fixed direction u (all pairwise inner
              products approximately equal to the target),
  "pair"    - columns come in strongly correlated pairs, the rest incoherent.

Solvers are written here from scratch (no external sparse-recovery package):
ISTA and FISTA (Beck-Teboulle 2009) for the LASSO
    min_x 0.5 ||A x - b||_2^2 + lam ||x||_1 .
"""

import hashlib
import json
import os
import platform
import sys
import time

import numpy as np
import scipy

HERE = os.path.dirname(os.path.abspath(__file__))

FIXTURE_SOURCE = "builtin_fallback_in_inv_sparse.py"
_ext = None
try:
    sys.path.insert(0, HERE)
    import inv_fixtures as _ext  # type: ignore
    if hasattr(_ext, "coherent_dictionary"):
        FIXTURE_SOURCE = "inv_fixtures.py"
    else:
        _ext = None
except Exception:
    _ext = None


# ----------------------------------------------------------------------------
# Dictionary fixture
# ----------------------------------------------------------------------------

def _normalize_cols(A):
    return A / np.linalg.norm(A, axis=0, keepdims=True)


def coherent_dictionary_local(m, n, mu_target, seed, mode="common"):
    rng = np.random.default_rng(seed)
    if mode == "common":
        u = rng.standard_normal(m)
        u /= np.linalg.norm(u)
        G = _normalize_cols(rng.standard_normal((m, n)))
        t = float(mu_target)
        A = np.sqrt(t) * u[:, None] + np.sqrt(1.0 - t) * G
    elif mode == "pair":
        assert n % 2 == 0
        B = _normalize_cols(rng.standard_normal((m, n // 2)))
        W = _normalize_cols(rng.standard_normal((m, n // 2)))
        c = float(mu_target)
        P = c * B + np.sqrt(max(1.0 - c * c, 0.0)) * W
        A = np.empty((m, n))
        A[:, 0::2] = B
        A[:, 1::2] = P
    else:
        raise ValueError(mode)
    return _normalize_cols(A)


USED_SOURCE = {}
_DICT_CACHE = {}


def get_dictionary(m, n, mu_target, seed, mode):
    """mode 'common' prefers inv_fixtures.coherent_dictionary (its Fixture is a
    dict subclass with key 'A'); 'pair' and 'twin' are local constructions."""
    key = (m, n, mu_target, seed, mode)
    if key in _DICT_CACHE:
        return _DICT_CACHE[key]
    A = None
    if _ext is not None and mode == "common":
        try:
            out = _ext.coherent_dictionary(m, n, mu_target, seed)
            A = np.asarray(out["A"] if isinstance(out, dict) else out[0], float)
            if A.shape != (m, n):
                A = None
            else:
                USED_SOURCE[mode] = "inv_fixtures.py"
        except Exception as exc:
            print("inv_fixtures fallback:", repr(exc))
            A = None
    if A is None:
        USED_SOURCE.setdefault(mode, "builtin_fallback_in_inv_sparse.py")
        A = coherent_dictionary_local(m, n, mu_target, seed, mode)
    A = _normalize_cols(A)
    if len(_DICT_CACHE) < 400:
        _DICT_CACHE[key] = A
    return A


def mutual_coherence(A):
    G = np.abs(A.T @ A)
    np.fill_diagonal(G, 0.0)
    return float(G.max())


def typical_coherence(A):
    G = np.abs(A.T @ A)
    iu = np.triu_indices(A.shape[1], 1)
    return float(np.median(G[iu]))


def make_signal(n, k, rng):
    supp = np.sort(rng.choice(n, size=k, replace=False))
    x = np.zeros(n)
    x[supp] = rng.choice([-1.0, 1.0], size=k) * rng.uniform(1.0, 3.0, size=k)
    return x, set(int(i) for i in supp)


def tropp_erc(A, supp):
    """Tropp 2004 exact recovery coefficient: 1 - max_{j not in S} ||A_S^+ a_j||_1.

    ERC > 0 is a sufficient condition for L1 (basis pursuit) to recover the
    support in the noiseless case.  It is sufficient, not necessary.
    """
    S = sorted(supp)
    comp = [j for j in range(A.shape[1]) if j not in supp]
    As = A[:, S]
    P = np.linalg.pinv(As)
    C = np.abs(P @ A[:, comp]).sum(axis=0)
    return float(1.0 - C.max())


# ----------------------------------------------------------------------------
# ISTA / FISTA (multi-lambda: X has one column per lambda)
# ----------------------------------------------------------------------------

def soft(X, tau):
    return np.sign(X) * np.maximum(np.abs(X) - tau, 0.0)


def fista_multi(A, b, lams, n_iter, L=None, checkpoints=None):
    n = A.shape[1]
    lams = np.asarray(lams, float)
    if L is None:
        L = float(np.linalg.norm(A, 2) ** 2)
    step = 1.0 / L
    tau = lams * step
    X = np.zeros((n, len(lams)))
    Y = X.copy()
    t = 1.0
    snaps = {}
    for it in range(1, n_iter + 1):
        Gr = A.T @ (A @ Y - b[:, None])
        Xn = soft(Y - step * Gr, tau[None, :])
        tn = 0.5 * (1.0 + np.sqrt(1.0 + 4.0 * t * t))
        Y = Xn + ((t - 1.0) / tn) * (Xn - X)
        X, t = Xn, tn
        if checkpoints and it in checkpoints:
            snaps[it] = X.copy()
    return X, snaps, L


def ista_multi(A, b, lams, n_iter, L=None, checkpoints=None):
    n = A.shape[1]
    lams = np.asarray(lams, float)
    if L is None:
        L = float(np.linalg.norm(A, 2) ** 2)
    step = 1.0 / L
    tau = lams * step
    X = np.zeros((n, len(lams)))
    snaps = {}
    for it in range(1, n_iter + 1):
        Gr = A.T @ (A @ X - b[:, None])
        X = soft(X - step * Gr, tau[None, :])
        if checkpoints and it in checkpoints:
            snaps[it] = X.copy()
    return X, snaps, L


def objective(A, b, x, lam):
    r = A @ x - b
    return float(0.5 * r @ r + lam * np.abs(x).sum())


# ----------------------------------------------------------------------------
# Main sweep
# ----------------------------------------------------------------------------

M, N = 100, 200
MUS = [0.1, 0.5, 0.7, 0.9, 0.99]
KS = [3, 5, 8, 12, 20, 30, 40]
NOISES = [0.0, 0.01]
TRIALS = 10
LAM_C = [0.5, 0.2, 0.1, 0.05, 0.02, 0.005]
FIXED_C = 0.05
N_ITER = 2000
SUPP_TOL = 1e-8


def one_cell(mu_t, k, noise, mode):
    rows = []
    for trial in range(TRIALS):
        # dictionary seed depends only on (mu_target, trial) so the same
        # dictionary is reused across k and noise level (cheaper, less variance)
        dseed = 700000 + int(round(mu_t * 1000)) * 100 + trial
        seed = dseed * 10 + k + (7 if noise > 0 else 0)
        A = get_dictionary(M, N, mu_t, dseed, mode)
        mu = mutual_coherence(A)
        rng = np.random.default_rng(seed + 1)
        x_true, supp = make_signal(N, k, rng)
        b_clean = A @ x_true
        if noise > 0:
            e = rng.standard_normal(M)
            e *= noise * np.linalg.norm(b_clean) / np.linalg.norm(e)
        else:
            e = np.zeros(M)
        b = b_clean + e
        lams = [c * float(np.abs(A.T @ b).max()) for c in LAM_C]
        X, _, _ = fista_multi(A, b, lams, N_ITER)
        per_lam = []
        for i, c in enumerate(LAM_C):
            xh = X[:, i]
            s_exact = set(int(j) for j in np.nonzero(np.abs(xh) > SUPP_TOL)[0])
            topk = set(int(j) for j in np.argsort(-np.abs(xh))[:k])
            per_lam.append({
                "c": c,
                "relerr": float(np.linalg.norm(xh - x_true) / np.linalg.norm(x_true)),
                "nnz": len(s_exact),
                "exact_support": s_exact == supp,
                "topk_support": topk == supp,
                "supp_recall": len(topk & supp) / k,
            })
        best = int(np.argmin([p["relerr"] for p in per_lam]))
        fixed = LAM_C.index(FIXED_C)
        rows.append({
            "trial": trial, "seed": seed, "dict_seed": dseed, "mu_measured": mu,
            "mu_typical": typical_coherence(A),
            "erc": tropp_erc(A, supp),
            "donoho_elad_k_star": 0.5 * (1.0 + 1.0 / mu),
            "oracle_lambda_c": per_lam[best]["c"],
            "oracle_relerr": per_lam[best]["relerr"],
            "oracle_exact_support": per_lam[best]["exact_support"],
            "oracle_topk_support": per_lam[best]["topk_support"],
            "oracle_supp_recall": per_lam[best]["supp_recall"],
            "fixed_relerr": per_lam[fixed]["relerr"],
            "fixed_topk_support": per_lam[fixed]["topk_support"],
            "fixed_exact_support": per_lam[fixed]["exact_support"],
            "fixed_nnz": per_lam[fixed]["nnz"],
            "per_lambda": per_lam,
        })
    agg = {
        "mode": mode, "mu_target": mu_t, "k": k, "noise": noise,
        "trials": TRIALS,
        "mu_measured_mean": float(np.mean([r["mu_measured"] for r in rows])),
        "mu_typical_mean": float(np.mean([r["mu_typical"] for r in rows])),
        "erc_mean": float(np.mean([r["erc"] for r in rows])),
        "erc_positive_fraction": float(np.mean([r["erc"] > 0 for r in rows])),
        "donoho_elad_k_star_mean": float(np.mean([r["donoho_elad_k_star"] for r in rows])),
        "oracle_topk_support_rate": float(np.mean([r["oracle_topk_support"] for r in rows])),
        "oracle_exact_support_rate": float(np.mean([r["oracle_exact_support"] for r in rows])),
        "oracle_supp_recall_mean": float(np.mean([r["oracle_supp_recall"] for r in rows])),
        "oracle_relerr_median": float(np.median([r["oracle_relerr"] for r in rows])),
        "fixed_topk_support_rate": float(np.mean([r["fixed_topk_support"] for r in rows])),
        "fixed_relerr_median": float(np.median([r["fixed_relerr"] for r in rows])),
        "fixed_nnz_median": float(np.median([r["fixed_nnz"] for r in rows])),
    }
    return agg, rows


TWIN_CS = [0.90, 0.99, 0.999]
TWIN_TRIALS = 20


def twin_sweep():
    """Adversarial construction (local, not from inv_fixtures): columns come in
    pairs with inner product exactly c.  Two support regimes:
      'twin_off_support' - exactly one member of each of the k pairs is active,
                           so every active atom has a near-duplicate that is
                           NOT in the support (Tropp's ERC is violated by
                           construction as c -> 1);
      'twin_in_support'  - both members of one pair are active.
    Reported: top-k support recovery and the twin-confusion rate, i.e. how often
    the recovered top-k set contains an inactive twin in place of a true atom.
    """
    out, agg = [], []
    k = 5
    for c in TWIN_CS:
        for regime in ("twin_off_support", "twin_in_support"):
            for noise in (0.0, 0.01):
                rows = []
                for trial in range(TWIN_TRIALS):
                    dseed = 500000 + int(c * 1000) * 100 + trial
                    A = get_dictionary(M, N, c, dseed, "pair")
                    rng = np.random.default_rng(dseed * 7 + 3)
                    npair = N // 2
                    if regime == "twin_off_support":
                        pairs = rng.choice(npair, size=k, replace=False)
                        supp_idx = [2 * p + int(rng.integers(2)) for p in pairs]
                    else:
                        p0 = int(rng.integers(npair))
                        others = [p for p in rng.permutation(npair) if p != p0][: k - 2]
                        supp_idx = [2 * p0, 2 * p0 + 1] + [
                            2 * p + int(rng.integers(2)) for p in others]
                    supp = set(int(i) for i in supp_idx)
                    x_true = np.zeros(N)
                    x_true[sorted(supp)] = rng.choice([-1.0, 1.0], size=k) * \
                        rng.uniform(1.0, 3.0, size=k)
                    b_clean = A @ x_true
                    b = b_clean.copy()
                    if noise > 0:
                        e = rng.standard_normal(M)
                        b = b_clean + e * noise * np.linalg.norm(b_clean) / np.linalg.norm(e)
                    lams = [cc * float(np.abs(A.T @ b).max()) for cc in LAM_C]
                    X, _, _ = fista_multi(A, b, lams, N_ITER)
                    per = []
                    for i, cc in enumerate(LAM_C):
                        xh = X[:, i]
                        topk = set(int(j) for j in np.argsort(-np.abs(xh))[:k])
                        twins = set(j ^ 1 for j in supp) - supp
                        per.append({
                            "c": cc,
                            "relerr": float(np.linalg.norm(xh - x_true) /
                                            np.linalg.norm(x_true)),
                            "topk_support": topk == supp,
                            "twin_confusion": len(topk & twins) > 0,
                            "supp_recall": len(topk & supp) / k,
                        })
                    best = int(np.argmin([p["relerr"] for p in per]))
                    fx = LAM_C.index(FIXED_C)
                    rows.append({
                        "trial": trial, "dict_seed": dseed,
                        "mu_measured": mutual_coherence(A),
                        "oracle_topk_support": per[best]["topk_support"],
                        "oracle_twin_confusion": per[best]["twin_confusion"],
                        "oracle_supp_recall": per[best]["supp_recall"],
                        "oracle_relerr": per[best]["relerr"],
                        "fixed_topk_support": per[fx]["topk_support"],
                        "fixed_twin_confusion": per[fx]["twin_confusion"],
                        "fixed_relerr": per[fx]["relerr"],
                    })
                a = {
                    "twin_c": c, "regime": regime, "noise": noise, "k": k,
                    "trials": TWIN_TRIALS,
                    "mu_measured_mean": float(np.mean([r["mu_measured"] for r in rows])),
                    "oracle_topk_support_rate": float(np.mean(
                        [r["oracle_topk_support"] for r in rows])),
                    "oracle_twin_confusion_rate": float(np.mean(
                        [r["oracle_twin_confusion"] for r in rows])),
                    "oracle_supp_recall_mean": float(np.mean(
                        [r["oracle_supp_recall"] for r in rows])),
                    "oracle_relerr_median": float(np.median(
                        [r["oracle_relerr"] for r in rows])),
                    "fixed_topk_support_rate": float(np.mean(
                        [r["fixed_topk_support"] for r in rows])),
                    "fixed_twin_confusion_rate": float(np.mean(
                        [r["fixed_twin_confusion"] for r in rows])),
                    "fixed_relerr_median": float(np.median(
                        [r["fixed_relerr"] for r in rows])),
                }
                agg.append(a)
                out.append({"cell": a, "rows": rows})
                print(f"twin c={c} {regime} noise={noise}: mu={a['mu_measured_mean']:.4f} "
                      f"oracle_topk={a['oracle_topk_support_rate']:.2f} "
                      f"twinconf={a['oracle_twin_confusion_rate']:.2f} "
                      f"fixed_topk={a['fixed_topk_support_rate']:.2f} "
                      f"relerr={a['oracle_relerr_median']:.3f}")
    return {"aggregate": agg, "per_trial": out}


def ista_vs_fista():
    """Is a support-recovery failure statistical (L1 cannot) or algorithmic
    (the solver has not converged)?  Fixed lambda, increasing budgets."""
    out = []
    cps = [100, 1000, 10000]
    for mu_t in (0.1, 0.5, 0.9):
        for trial in range(5):
            seed = 900000 + int(mu_t * 1000) * 10 + trial
            A = get_dictionary(M, N, mu_t, seed, "common")
            rng = np.random.default_rng(seed + 1)
            x_true, supp = make_signal(N, 5, rng)
            b = A @ x_true
            lam = [FIXED_C * float(np.abs(A.T @ b).max())]
            _, sF, L = fista_multi(A, b, lam, 10000, checkpoints=cps)
            _, sI, _ = ista_multi(A, b, lam, 10000, L=L, checkpoints=cps)
            best = min(objective(A, b, sF[10000][:, 0], lam[0]),
                       objective(A, b, sI[10000][:, 0], lam[0]))
            row = {"mu_target": mu_t, "trial": trial, "seed": seed,
                   "mu_measured": mutual_coherence(A), "lam": lam[0]}
            for tag, sn in (("fista", sF), ("ista", sI)):
                for c in cps:
                    x = sn[c][:, 0]
                    topk = set(int(j) for j in np.argsort(-np.abs(x))[:5])
                    row[f"{tag}_obj_gap_{c}"] = objective(A, b, x, lam[0]) - best
                    row[f"{tag}_topk_{c}"] = topk == supp
                    row[f"{tag}_relerr_{c}"] = float(
                        np.linalg.norm(x - x_true) / np.linalg.norm(x_true))
            out.append(row)
    agg = {}
    for mu_t in (0.1, 0.5, 0.9):
        rs = [r for r in out if r["mu_target"] == mu_t]
        agg[str(mu_t)] = {
            f"{tag}_topk_rate_{c}": float(np.mean([r[f"{tag}_topk_{c}"] for r in rs]))
            for tag in ("ista", "fista") for c in cps
        }
        agg[str(mu_t)].update({
            f"{tag}_obj_gap_median_{c}": float(np.median([r[f"{tag}_obj_gap_{c}"] for r in rs]))
            for tag in ("ista", "fista") for c in cps
        })
    return {"per_trial": out, "aggregate": agg, "checkpoints": cps,
            "k": 5, "lambda_rule": f"{FIXED_C} * ||A^T b||_inf", "noise": 0.0}


def _sha_of(fname):
    p = os.path.join(HERE, fname)
    if not os.path.exists(p):
        return None
    with open(p, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def main():
    t0 = time.time()
    res = {
        "script": "inv_sparse.py",
        "role": "iterative-and-sparse-investigator",
        "common": "C (failure maps for numerical methods on ill-conditioned inverse problems)",
        "fixture_source": FIXTURE_SOURCE,
        "env": {"python": sys.version.split()[0], "numpy": np.__version__,
                "scipy": scipy.__version__, "platform": platform.platform()},
        "setup": {"m": M, "n": N, "mus": MUS, "ks": KS, "noises": NOISES,
                  "trials": TRIALS, "lambda_c_grid": LAM_C, "fixed_c": FIXED_C,
                  "fista_iterations": N_ITER, "support_tol": SUPP_TOL,
                  "amplitude_range": [1.0, 3.0]},
        "inv_fixtures_sha256": _sha_of("inv_fixtures.py"),
        "grid": [], "grid_per_trial": [], "twin_sweep": None,
        "ista_vs_fista": None,
    }
    for mode in ("common", "pair"):
        for mu_t in MUS:
            for k in KS:
                for noise in NOISES:
                    if mode == "pair" and noise > 0:
                        continue  # keep the secondary construction cheap
                    agg, rows = one_cell(mu_t, k, noise, mode)
                    res["grid"].append(agg)
                    res["grid_per_trial"].append({"mode": mode, "mu_target": mu_t,
                                                  "k": k, "noise": noise, "rows": rows})
                    print(f"{mode} mu~{agg['mu_measured_mean']:.3f} "
                          f"(typ {agg['mu_typical_mean']:.3f}) k={k} noise={noise}: "
                          f"oracle_topk={agg['oracle_topk_support_rate']:.2f} "
                          f"oracle_exact={agg['oracle_exact_support_rate']:.2f} "
                          f"fixed_topk={agg['fixed_topk_support_rate']:.2f} "
                          f"ERC>0={agg['erc_positive_fraction']:.2f} "
                          f"relerr_med={agg['oracle_relerr_median']:.3f}")
    res["twin_sweep"] = twin_sweep()
    res["fixture_source_used"] = dict(USED_SOURCE)
    res["ista_vs_fista"] = ista_vs_fista()
    print("ista_vs_fista:", json.dumps(res["ista_vs_fista"]["aggregate"], indent=1))
    res["runtime_seconds"] = round(time.time() - t0, 2)

    path = os.path.join(HERE, "inv_sparse_output.json")
    with open(path, "w") as f:
        json.dump(res, f, indent=1, sort_keys=True)
    with open(path, "rb") as f:
        h = hashlib.sha256(f.read()).hexdigest()
    print("wrote", path, "sha256", h, "in", res["runtime_seconds"], "s")


if __name__ == "__main__":
    main()
