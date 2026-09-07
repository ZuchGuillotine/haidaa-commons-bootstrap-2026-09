"""Inert fixture library for ill-conditioned linear (and one mildly nonlinear)
inverse problems.

Purpose
-------
Provide small, fully deterministic, self-contained test problems with documented
conditioning, so that failure claims about numerical methods can be attached to a
named fixture instead of to an ad-hoc script.

Provenance / honesty note
-------------------------
`shaw`, `phillips` and `heat` are *re-implementations* of classical Fredholm
first-kind test problems in the style of P. C. Hansen's Regularization Tools
(Hansen 1994, Numer. Algorithms 6:1-35; Hansen 2007, Numer. Algorithms 46:189-194).
They are written here from the published integral equations and a stated quadrature
rule. They have NOT been compared element-by-element against Hansen's MATLAB code,
so do not assume numerical identity with `shaw.m` / `phillips.m` / `heat.m`.
Singular values and condition numbers reported by this module are the ones measured
for THIS code. Any claim that reuses them should name this file and the bundle hash.

Conventions
-----------
* All linear fixtures return an object with attributes/keys:
      A         (m, n) float64 matrix
      x_true    (n,)   float64 ground truth
      b_clean   (m,)   float64 = A @ x_true  (noise-free right-hand side)
      meta      dict   name, parameters, grids, conditioning descriptors
  The object supports both `f.A` and `f["A"]`.
* Noise is never added here. Use `noise.py`.
* Nothing is random unless a `seed` argument is present; every random fixture takes
  an explicit integer seed and uses numpy.random.default_rng(seed).
* Import has no side effects. Building the shipped bundle is `python inv_fixtures.py`.

Dependencies: numpy, scipy only.
"""

from __future__ import annotations

import hashlib
import json
import os

import numpy as np

__all__ = [
    "Fixture",
    "shaw",
    "phillips",
    "heat",
    "deconv_gauss",
    "coherent_dictionary",
    "expsum",
    "svd_descriptor",
    "mutual_coherence",
    "build_bundle",
    "load_bundle",
    "BUNDLE_NAME",
    "MANIFEST_NAME",
]

BUNDLE_NAME = "inv_fixtures_bundle.npz"
MANIFEST_NAME = "inv_fixtures_manifest.json"


class Fixture(dict):
    """dict with attribute access, so f.A and f['A'] both work."""

    def __getattr__(self, item):
        try:
            return self[item]
        except KeyError as exc:  # pragma: no cover
            raise AttributeError(item) from exc

    def __repr__(self):  # pragma: no cover
        name = self.get("meta", {}).get("name", "fixture")
        A = self.get("A")
        shp = None if A is None else A.shape
        return f"<Fixture {name} A={shp}>"


# ----------------------------------------------------------------------------
# conditioning descriptors
# ----------------------------------------------------------------------------


def svd_descriptor(A, tail_frac=0.75):
    """Measure conditioning and classify singular value decay.

    Fits, over the leading `tail_frac` fraction of the spectrum (indices where the
    singular values are still above the double-precision noise floor):
        exponential model   log s_i = a - b * i
        algebraic  model    log s_i = a - p * log i
    and reports the coefficient of determination of each fit. The class with the
    larger R^2 is reported as `decay_class`; when the two R^2 values are within
    0.02 the class is reported as "mixed". This is a descriptive heuristic, not a
    theorem about the underlying operator.

    Returns a dict of plain Python floats/strings (JSON-serializable).
    """
    A = np.asarray(A, dtype=np.float64)
    s = np.linalg.svd(A, compute_uv=False)
    s = np.maximum(s, 0.0)
    smax = float(s[0])
    smin = float(s[-1])
    cond = float(smax / smin) if smin > 0 else float("inf")

    # numerical rank at the standard tolerance
    tol = max(A.shape) * np.finfo(np.float64).eps * smax
    num_rank = int(np.sum(s > tol))

    # usable part of the spectrum for fitting: above the double precision floor
    floor = smax * 1e-15
    keep = s > floor
    idx = np.arange(1, s.size + 1)[keep]
    sv = s[keep]
    ncut = max(3, int(np.ceil(tail_frac * sv.size)))
    idx = idx[:ncut]
    sv = sv[:ncut]

    y = np.log(sv)

    def _r2(x):
        Xd = np.vstack([x, np.ones_like(x)]).T
        coef, *_ = np.linalg.lstsq(Xd, y, rcond=None)
        pred = Xd @ coef
        ss_res = float(np.sum((y - pred) ** 2))
        ss_tot = float(np.sum((y - y.mean()) ** 2))
        r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
        return float(coef[0]), float(r2)

    b_exp, r2_exp = _r2(idx.astype(np.float64))
    p_alg, r2_alg = _r2(np.log(idx.astype(np.float64)))

    if abs(r2_exp - r2_alg) <= 0.02:
        cls = "mixed"
    elif r2_exp > r2_alg:
        cls = "exponential"
    else:
        cls = "algebraic"

    return {
        "shape": [int(A.shape[0]), int(A.shape[1])],
        "smax": smax,
        "smin": smin,
        "cond2": cond,
        "numerical_rank": num_rank,
        "n_fitted": int(sv.size),
        "exp_rate_per_index": -b_exp,
        "r2_exponential": r2_exp,
        "alg_power": -p_alg,
        "r2_algebraic": r2_alg,
        "decay_class": cls,
        "note": (
            "cond2 above ~1e16 is at or past the double-precision resolution limit; "
            "treat it as 'numerically singular', not as an accurate number."
        ),
    }


def mutual_coherence(A):
    """Max |<a_i, a_j>| over i != j for column-normalized A."""
    A = np.asarray(A, dtype=np.float64)
    norms = np.linalg.norm(A, axis=0)
    Q = A / norms
    G = np.abs(Q.T @ Q)
    np.fill_diagonal(G, 0.0)
    return float(G.max())


# ----------------------------------------------------------------------------
# Fredholm first-kind fixtures
# ----------------------------------------------------------------------------


def shaw(n):
    """Shaw's one-dimensional image restoration model (Shaw 1972).

    Integral equation on s, t in [-pi/2, pi/2]:
        K(s,t) = (cos(s) + cos(t))^2 * ( sin(u)/u )^2,  u = pi*(sin s + sin t)
    discretized by the midpoint rule with n equal cells (h = pi/n), i.e.
        A[i,j] = h * K(s_i, t_j),  s_i = t_i = -pi/2 + (i+0.5)*h.
    The limit u -> 0 is taken as (sin u / u) -> 1.

    Ground truth: the standard two-Gaussian profile
        x(t) = 2*exp(-6*(t-0.8)^2) + exp(-2*(t+0.5)^2).

    n must be even. Severely ill-posed: exponentially decaying singular values.
    """
    n = int(n)
    if n % 2 != 0:
        raise ValueError("shaw: n must be even")
    h = np.pi / n
    t = -np.pi / 2.0 + (np.arange(n) + 0.5) * h

    co = np.cos(t)
    u = np.pi * (np.sin(t)[:, None] + np.sin(t)[None, :])
    sinc = np.ones_like(u)
    nz = u != 0.0
    sinc[nz] = np.sin(u[nz]) / u[nz]
    A = h * ((co[:, None] + co[None, :]) ** 2) * (sinc ** 2)

    x_true = 2.0 * np.exp(-6.0 * (t - 0.8) ** 2) + 1.0 * np.exp(-2.0 * (t + 0.5) ** 2)
    b_clean = A @ x_true

    meta = {
        "name": f"shaw_n{n}",
        "family": "shaw",
        "n": n,
        "quadrature": "midpoint, h = pi/n",
        "domain": [-np.pi / 2, np.pi / 2],
        "x_true_desc": "2*exp(-6*(t-0.8)^2) + exp(-2*(t+0.5)^2)",
        "reference": "Shaw 1972 J. Math. Anal. Appl.; Hansen 1994 Regularization Tools",
        "identity_to_hansen_matlab": "not verified",
    }
    meta["svd"] = svd_descriptor(A)
    return Fixture(A=A, x_true=x_true, b_clean=b_clean, grid_s=t, grid_t=t, meta=meta)


def phillips(n):
    """Phillips' test problem (Phillips 1962).

    phi(x) = 1 + cos(pi*x/3) for |x| < 3, else 0.
    Kernel K(s,t) = phi(s-t) on s, t in [-6, 6]; true solution x(t) = phi(t).
    Midpoint quadrature with h = 12/n:
        A[i,j] = h * phi(s_i - t_j),  s_i = t_i = -6 + (i+0.5)*h.
    b_clean = A @ x_true. The analytic right-hand side
        g(s) = (6-|s|)*(1 + 0.5*cos(pi*s/3)) + (9/(2*pi))*sin(pi*|s|/3),  |s| <= 6
    is also returned as meta-level array `b_analytic` for reference; the fixture's
    b_clean is the discrete A @ x_true so that the discrete system is exactly
    consistent.

    n must be a multiple of 4. Mildly-to-moderately ill-posed.
    """
    n = int(n)
    if n % 4 != 0:
        raise ValueError("phillips: n must be a multiple of 4")
    h = 12.0 / n
    t = -6.0 + (np.arange(n) + 0.5) * h

    def phi(x):
        out = np.zeros_like(x)
        m = np.abs(x) < 3.0
        out[m] = 1.0 + np.cos(np.pi * x[m] / 3.0)
        return out

    A = h * phi(t[:, None] - t[None, :])
    x_true = phi(t)
    b_clean = A @ x_true

    s = t
    b_analytic = (6.0 - np.abs(s)) * (1.0 + 0.5 * np.cos(np.pi * s / 3.0)) + (
        9.0 / (2.0 * np.pi)
    ) * np.sin(np.pi * np.abs(s) / 3.0)
    b_analytic = np.where(np.abs(s) <= 6.0, b_analytic, 0.0)

    meta = {
        "name": f"phillips_n{n}",
        "family": "phillips",
        "n": n,
        "quadrature": "midpoint, h = 12/n",
        "domain": [-6.0, 6.0],
        "x_true_desc": "phi(t) = 1 + cos(pi*t/3) on |t|<3, 0 otherwise",
        "reference": "Phillips 1962 J. ACM 9:84-97; Hansen 1994 Regularization Tools",
        "identity_to_hansen_matlab": "not verified",
    }
    meta["svd"] = svd_descriptor(A)
    return Fixture(
        A=A,
        x_true=x_true,
        b_clean=b_clean,
        b_analytic=b_analytic,
        grid_s=s,
        grid_t=t,
        meta=meta,
    )


def heat(n, kappa=1.0):
    """Inverse heat equation (Volterra convolution equation of the first kind).

    K(s,t) = k(s-t) with
        k(u) = u^{-3/2} / (2*kappa*sqrt(pi)) * exp(-1/(4*kappa^2*u)),  u > 0
    on s, t in (0, 1]. Quadrature: the unknown is sampled at cell midpoints
    t_j = (j+0.5)*h and the data at cell right endpoints s_i = (i+1)*h, h = 1/n, so
        A[i,j] = h * k(s_i - t_j) = h * k((i-j+0.5)*h) for j <= i, else 0,
    a lower-triangular Toeplitz matrix with a strictly positive diagonal. (Sampling
    both grids at midpoints would put k(0)=0 on the diagonal and make A exactly
    singular by construction, which is a discretization artifact, not ill-posedness.)

    Ground truth (documented here, piecewise smooth with a kink; the shape follows
    the profile used in Regularization Tools but the exact breakpoints are as coded
    below and were not checked against the MATLAB source):
        with v = 20*t,
        x = 0.75*v^2/4                 for v < 2
        x = 0.75 + (v-2)*(3-v)         for 2 <= v < 3
        x = 0.75*exp(-2*(v-3))         for 3 <= v <= 10
        x = 0                          for v > 10   (i.e. the second half of the grid)

    kappa controls the severity: smaller kappa is more ill-posed. n must be even.
    """
    n = int(n)
    if n % 2 != 0:
        raise ValueError("heat: n must be even")
    h = 1.0 / n
    t = (np.arange(n) + 0.5) * h
    s = (np.arange(n) + 1.0) * h

    def k(u):
        out = np.zeros_like(u)
        m = u > 0
        uu = u[m]
        out[m] = uu ** (-1.5) / (2.0 * kappa * np.sqrt(np.pi)) * np.exp(
            -1.0 / (4.0 * kappa ** 2 * uu)
        )
        return out

    D = s[:, None] - t[None, :]
    A = np.tril(h * k(D))

    v = 20.0 * t
    x_true = np.zeros(n)
    m1 = v < 2.0
    m2 = (v >= 2.0) & (v < 3.0)
    m3 = (v >= 3.0) & (v <= 10.0)
    x_true[m1] = 0.75 * v[m1] ** 2 / 4.0
    x_true[m2] = 0.75 + (v[m2] - 2.0) * (3.0 - v[m2])
    x_true[m3] = 0.75 * np.exp(-2.0 * (v[m3] - 3.0))
    b_clean = A @ x_true

    meta = {
        "name": f"heat_n{n}_kappa{kappa:g}",
        "family": "heat",
        "n": n,
        "kappa": float(kappa),
        "quadrature": (
            "h = 1/n; unknown at midpoints t_j=(j+0.5)h, data at endpoints s_i=(i+1)h; "
            "lower-triangular Toeplitz"
        ),
        "domain": [0.0, 1.0],
        "x_true_desc": "piecewise smooth ramp-peak-decay, zero on the second half",
        "reference": (
            "inverse heat conduction / Volterra first-kind equation; the discretization "
            "follows Hansen 1994 Regularization Tools (heat.m). Carasso 1982 SIAM J. Appl. "
            "Math. is often cited for this problem but that citation is unverified here."
        ),
        "identity_to_hansen_matlab": "not verified",
    }
    meta["svd"] = svd_descriptor(A)
    return Fixture(A=A, x_true=x_true, b_clean=b_clean, grid_s=s, grid_t=t, meta=meta)


def deconv_gauss(n, sigma=0.03):
    """1D deconvolution with a Gaussian point spread function on [0, 1].

        K(s,t) = exp(-(s-t)^2 / (2*sigma^2)) / (sigma*sqrt(2*pi))
        A[i,j] = h * K(s_i, t_j),  h = 1/n, t_i = (i+0.5)*h

    Ground truth: a deliberately piecewise-smooth signal with three distinct
    features, so that both smoothing bias and ringing are visible:
        1. a rectangular box of height 1 on [0.12, 0.28]
        2. a smooth Gaussian bump  0.8*exp(-((t-0.5)/0.05)^2)
        3. a linear ramp on [0.68, 0.88] rising from 0 to 0.6 with a jump back to 0
    sigma sets the severity: larger sigma -> faster singular value decay.
    """
    n = int(n)
    h = 1.0 / n
    t = (np.arange(n) + 0.5) * h
    D = t[:, None] - t[None, :]
    A = h * np.exp(-(D ** 2) / (2.0 * sigma ** 2)) / (sigma * np.sqrt(2.0 * np.pi))

    x_true = np.zeros(n)
    x_true += np.where((t >= 0.12) & (t <= 0.28), 1.0, 0.0)
    x_true += 0.8 * np.exp(-(((t - 0.5) / 0.05) ** 2))
    ramp = (t - 0.68) / (0.88 - 0.68) * 0.6
    x_true += np.where((t >= 0.68) & (t <= 0.88), ramp, 0.0)

    b_clean = A @ x_true

    meta = {
        "name": f"deconv_gauss_n{n}_sigma{sigma:g}",
        "family": "deconv_gauss",
        "n": n,
        "sigma": float(sigma),
        "quadrature": "midpoint, h = 1/n",
        "domain": [0.0, 1.0],
        "x_true_desc": "box on [0.12,0.28] + Gaussian bump at 0.5 + ramp on [0.68,0.88]",
        "reference": "standard Gaussian-blur deconvolution model; original to this library",
    }
    meta["svd"] = svd_descriptor(A)
    return Fixture(A=A, x_true=x_true, b_clean=b_clean, grid_s=t, grid_t=t, meta=meta)


# ----------------------------------------------------------------------------
# coherent dictionary
# ----------------------------------------------------------------------------


def _dictionary_from_rho(m, n, rho, rng_seed):
    """Columns a_j = sqrt(rho)*u + sqrt(1-rho)*w_j, w_j random unit vectors in the
    orthogonal complement of u. Then <a_i,a_j> = rho + (1-rho)*<w_i,w_j>, so rho
    shifts every pairwise inner product upward by a controlled amount."""
    rng = np.random.default_rng(rng_seed)
    u = rng.standard_normal(m)
    u /= np.linalg.norm(u)
    W = rng.standard_normal((m, n))
    W -= np.outer(u, u @ W)  # project out u
    W /= np.linalg.norm(W, axis=0)
    A = np.sqrt(rho) * u[:, None] + np.sqrt(1.0 - rho) * W
    A /= np.linalg.norm(A, axis=0)
    return A


def coherent_dictionary(m, n, coherence_target, seed, k=5, amp_floor=1.0):
    """Overcomplete dictionary (m < n) with mutual coherence calibrated to a target,
    plus a k-sparse ground truth.

    Construction: every column shares a common direction u with weight sqrt(rho);
    rho is found by bisection so that the *measured* mutual coherence
    max_{i!=j} |<a_i,a_j>| matches `coherence_target` as closely as possible.
    The randomness (u, the complementary directions, the support and the
    amplitudes) is fully determined by `seed`.

    Floor: a dictionary of n random unit columns in R^m already has a mutual
    coherence of roughly sqrt(2*log(n^2)/m); targets below that floor are NOT
    reachable by this construction. In that case rho = 0 is used and the achieved
    coherence (which will exceed the target) is reported in meta["coherence"].
    Always read meta["coherence"], never assume the target was met.

    Ground truth: support drawn uniformly without replacement, amplitudes
    sign * (amp_floor + |N(0,1)|) so no nonzero is negligible.
    """
    m, n, k = int(m), int(n), int(k)
    if m >= n:
        raise ValueError("coherent_dictionary: want an overcomplete dictionary, m < n")

    lo, hi = 0.0, 0.999
    A = _dictionary_from_rho(m, n, lo, seed)
    mu_lo = mutual_coherence(A)
    best = (abs(mu_lo - coherence_target), lo, A, mu_lo)
    if mu_lo < coherence_target:
        for _ in range(40):
            mid = 0.5 * (lo + hi)
            Am = _dictionary_from_rho(m, n, mid, seed)
            mu = mutual_coherence(Am)
            cand = (abs(mu - coherence_target), mid, Am, mu)
            if cand[0] < best[0]:
                best = cand
            if mu < coherence_target:
                lo = mid
            else:
                hi = mid
            if abs(mu - coherence_target) < 1e-4:
                break
    _, rho, A, mu = best

    rng = np.random.default_rng(seed + 10_000)
    support = np.sort(rng.choice(n, size=k, replace=False))
    signs = rng.choice(np.array([-1.0, 1.0]), size=k)
    amps = signs * (amp_floor + np.abs(rng.standard_normal(k)))
    x_true = np.zeros(n)
    x_true[support] = amps
    b_clean = A @ x_true

    G = np.abs(A.T @ A)
    np.fill_diagonal(G, 0.0)

    meta = {
        "name": f"coherent_m{m}_n{n}_mu{coherence_target:g}_k{k}_s{seed}",
        "family": "coherent_dictionary",
        "m": m,
        "n": n,
        "k": k,
        "seed": int(seed),
        "coherence_target": float(coherence_target),
        "coherence": float(mu),
        "coherence_mean_offdiag": float(G.mean()),
        "rho_common_component": float(rho),
        "target_reached": bool(abs(mu - coherence_target) < 0.02),
        "exact_recovery_coefficient_bound": float((2.0 * k - 1.0) * mu),
        "donoho_elad_k_bound": float(0.5 * (1.0 + 1.0 / mu)),
        "support": support.tolist(),
        "reference": (
            "Donoho & Elad 2003 PNAS 100:2197-2202; Tropp 2004 IEEE Trans. Inf. Theory "
            "50:2231-2242 (coherence-based uniqueness/recovery conditions)"
        ),
        "note": (
            "The sufficient condition k < 0.5*(1+1/mu) is only sufficient; failure "
            "above it is expected but not guaranteed, and success above it is common."
        ),
    }
    return Fixture(A=A, x_true=x_true, b_clean=b_clean, support=support, meta=meta)


# ----------------------------------------------------------------------------
# nonlinear: exponential sum
# ----------------------------------------------------------------------------


LANCZOS_RATES = np.array([1.0, 3.0, 5.0])
LANCZOS_AMPS = np.array([0.0951, 0.8607, 1.5576])


def expsum(t, rates=None, amps=None):
    """Exponential-sum (Lanczos-type) model for nonlinear least squares.

        f(t; a, r) = sum_k a_k * exp(-r_k * t)

    Defaults are the amplitudes/rates of the classical Lanczos exponential-fitting
    problem used as the NIST StRD `Lanczos` data generator:
        a = (0.0951, 0.8607, 1.5576), r = (1, 3, 5).
    (Reported from memory of the NIST StRD certified values; verify against
    itl.nist.gov/div898/strd before quoting the constants as authoritative.)

    Returns a Fixture with
        model(theta)     -> (len(t),) prediction; theta = [a_1..a_K, r_1..r_K]
        jac(theta)       -> (len(t), 2K) analytic Jacobian
        residual(theta, b) -> model(theta) - b
        theta_true, x_true (= theta_true), b_clean, t
    Parameter ordering is amplitudes first, then rates. The problem is notoriously
    ill-conditioned in the rates: nearby rates are nearly unidentifiable from
    noisy samples, which is exactly the initialization-dependence regime of
    interest.
    """
    t = np.asarray(t, dtype=np.float64)
    rates = LANCZOS_RATES.copy() if rates is None else np.asarray(rates, np.float64)
    amps = LANCZOS_AMPS.copy() if amps is None else np.asarray(amps, np.float64)
    if rates.shape != amps.shape:
        raise ValueError("expsum: rates and amps must have the same shape")
    K = rates.size
    theta_true = np.concatenate([amps, rates])

    def model(theta):
        theta = np.asarray(theta, dtype=np.float64)
        a, r = theta[:K], theta[K:]
        return np.exp(-np.outer(t, r)) @ a

    def jac(theta):
        theta = np.asarray(theta, dtype=np.float64)
        a, r = theta[:K], theta[K:]
        E = np.exp(-np.outer(t, r))          # (len(t), K)
        J = np.empty((t.size, 2 * K))
        J[:, :K] = E
        J[:, K:] = -E * (t[:, None] * a[None, :])
        return J

    def residual(theta, b):
        return model(theta) - np.asarray(b, dtype=np.float64)

    b_clean = model(theta_true)
    J0 = jac(theta_true)
    sJ = np.linalg.svd(J0, compute_uv=False)

    meta = {
        "name": f"expsum_K{K}_nt{t.size}",
        "family": "expsum",
        "K": int(K),
        "n_samples": int(t.size),
        "t_range": [float(t[0]), float(t[-1])],
        "rates": rates.tolist(),
        "amps": amps.tolist(),
        "param_order": "[a_1..a_K, r_1..r_K]",
        "jacobian_cond_at_truth": float(sJ[0] / sJ[-1]),
        "jacobian_svals_at_truth": sJ.tolist(),
        "reference": (
            "Lanczos 1956 Applied Analysis; NIST StRD nonlinear regression set "
            "(Lanczos1/2/3) - constants reported from memory, verify before quoting"
        ),
    }
    return Fixture(
        model=model,
        jac=jac,
        residual=residual,
        theta_true=theta_true,
        x_true=theta_true,
        b_clean=b_clean,
        t=t,
        meta=meta,
    )


# ----------------------------------------------------------------------------
# bundle
# ----------------------------------------------------------------------------


def _bundle_spec():
    """The exact list of instances shipped in inv_fixtures_bundle.npz."""
    spec = []
    for n in (64, 128, 256):
        spec.append(("shaw", {"n": n}))
    for n in (128, 256):
        spec.append(("phillips", {"n": n}))
    for n in (128, 256):
        spec.append(("heat", {"n": n, "kappa": 1.0}))
    spec.append(("heat", {"n": 128, "kappa": 5.0}))
    for sig in (0.02, 0.05):
        spec.append(("deconv_gauss", {"n": 256, "sigma": sig}))
    for mu in (0.5, 0.7, 0.9):
        for k in (3, 5, 8, 12):
            spec.append(
                (
                    "coherent_dictionary",
                    {"m": 128, "n": 256, "coherence_target": mu, "seed": 20250906, "k": k},
                )
            )
    return spec


def build_bundle(outdir=None):
    """Build the fixed-seed bundle + manifest. Returns the manifest dict."""
    outdir = outdir or os.path.dirname(os.path.abspath(__file__))
    arrays = {}
    entries = []

    builders = {
        "shaw": shaw,
        "phillips": phillips,
        "heat": heat,
        "deconv_gauss": deconv_gauss,
        "coherent_dictionary": coherent_dictionary,
    }

    for family, kw in _bundle_spec():
        f = builders[family](**kw)
        name = f.meta["name"]
        arrays[f"{name}/A"] = f.A
        arrays[f"{name}/x_true"] = f.x_true
        arrays[f"{name}/b_clean"] = f.b_clean
        if "grid_t" in f:
            arrays[f"{name}/grid_t"] = f.grid_t
        entry = {
            "key_prefix": name,
            "family": family,
            "call": f"{family}(" + ", ".join(f"{k}={v!r}" for k, v in kw.items()) + ")",
            "shape_A": [int(f.A.shape[0]), int(f.A.shape[1])],
            "meta": {k: v for k, v in f.meta.items() if k != "svd"},
        }
        if "svd" in f.meta:
            entry["conditioning"] = f.meta["svd"]
        else:
            entry["conditioning"] = svd_descriptor(f.A)
        entries.append(entry)

    # nonlinear fixture: store the sampled data, the model lives in code
    t = np.linspace(0.0, 1.15, 24)
    fe = expsum(t)
    arrays["expsum_K3_nt24/t"] = fe.t
    arrays["expsum_K3_nt24/b_clean"] = fe.b_clean
    arrays["expsum_K3_nt24/theta_true"] = fe.theta_true
    entries.append(
        {
            "key_prefix": "expsum_K3_nt24",
            "family": "expsum",
            "call": "expsum(np.linspace(0.0, 1.15, 24))",
            "shape_A": None,
            "meta": fe.meta,
            "conditioning": {
                "jacobian_cond_at_truth": fe.meta["jacobian_cond_at_truth"],
                "note": "nonlinear model; conditioning is of the Jacobian at the truth",
            },
        }
    )

    npz_path = os.path.join(outdir, BUNDLE_NAME)
    np.savez_compressed(npz_path, **arrays)
    sha = hashlib.sha256(open(npz_path, "rb").read()).hexdigest()

    import scipy  # noqa: F401  (recorded for provenance)
    import sys

    manifest = {
        "bundle": BUNDLE_NAME,
        "sha256": sha,
        "bytes": os.path.getsize(npz_path),
        "generator": "inv_fixtures.py",
        "generator_sha256": hashlib.sha256(
            open(os.path.abspath(__file__), "rb").read()
        ).hexdigest(),
        "build_command": "python inv_fixtures.py",
        "npz_key_convention": "<key_prefix>/<array>, arrays: A, x_true, b_clean, grid_t",
        "global_seed_note": (
            "Only the coherent_dictionary fixtures are random; they use "
            "numpy.random.default_rng(seed) with seed=20250906 (support/amplitudes use "
            "seed+10000). All other fixtures are deterministic functions of n/sigma/kappa."
        ),
        "environment": {
            "python": sys.version.split()[0],
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "platform": sys.platform,
        },
        "n_entries": len(entries),
        "entries": entries,
    }
    with open(os.path.join(outdir, MANIFEST_NAME), "w") as fh:
        json.dump(manifest, fh, indent=2, sort_keys=False)
    return manifest


def load_bundle(path=None):
    """Load the npz bundle. Returns a dict prefix -> {array_name: ndarray}."""
    path = path or os.path.join(os.path.dirname(os.path.abspath(__file__)), BUNDLE_NAME)
    z = np.load(path)
    out = {}
    for key in z.files:
        prefix, _, arr = key.rpartition("/")
        out.setdefault(prefix, {})[arr] = z[key]
    return out


if __name__ == "__main__":  # pragma: no cover
    man = build_bundle()
    print(f"wrote {man['bundle']} sha256={man['sha256']}")
    for e in man["entries"]:
        c = e["conditioning"]
        if "coherence" in e["meta"]:
            print(
                f"  {e['key_prefix']:40s} mu={e['meta']['coherence']:.4f} "
                f"(target {e['meta']['coherence_target']:.2f}, reached="
                f"{e['meta']['target_reached']}) cond2={c['cond2']:.4e}"
            )
        elif "cond2" in c:
            print(
                f"  {e['key_prefix']:36s} cond2={c['cond2']:.4e} "
                f"rank={c['numerical_rank']:4d} decay={c['decay_class']}"
            )
        elif "coherence" in e["meta"]:
            print(f"  {e['key_prefix']:36s} mu={e['meta']['coherence']:.4f}")
        else:
            print(f"  {e['key_prefix']:36s} {c}")
