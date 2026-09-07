"""Noise models for the ill-conditioned inverse-problem fixture library.

Companion to inv_fixtures.py. Every generator is seeded explicitly; nothing uses
global numpy random state.

Level convention (important - quote it whenever you report a noise level)
------------------------------------------------------------------------
`rel_level` is the noise standard deviation expressed as a fraction of the RMS of
the clean right-hand side:

    sigma = rel_level * ||b_clean||_2 / sqrt(m)

so that E||e||_2 / ||b_clean||_2 ~= rel_level. "1 % noise" therefore means
rel_level = 0.01 and an expected relative perturbation of about 1 % in the
2-norm. The AR(1) model is scaled to the SAME marginal standard deviation as the
white model, so white and correlated runs at the same rel_level carry the same
expected noise energy and differ only in the spectrum of the perturbation.

Mis-specification
-----------------
`misspecified_sigma` returns the value a parameter-choice rule is *told*, which is
the true sigma multiplied by a factor (e.g. 0.5 or 2.0). The data itself is always
generated with the true sigma. Discrepancy-principle style rules consume the told
value; the truth is kept for scoring.

Dependencies: numpy only.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "sigma_from_relative_level",
    "white_gaussian",
    "ar1_correlated",
    "ar1_covariance",
    "misspecified_sigma",
    "add_noise",
    "relative_error",
]


def sigma_from_relative_level(b_clean, rel_level):
    """sigma = rel_level * ||b_clean||_2 / sqrt(m). Returns a float."""
    b_clean = np.asarray(b_clean, dtype=np.float64)
    m = b_clean.size
    return float(rel_level * np.linalg.norm(b_clean) / np.sqrt(m))


def white_gaussian(b_clean, rel_level, seed):
    """i.i.d. N(0, sigma^2) perturbation.

    Returns (b_noisy, info) where info holds sigma_true, the realized noise norm,
    and the realized relative error, so a caller can report what actually happened
    rather than the nominal level.
    """
    b_clean = np.asarray(b_clean, dtype=np.float64)
    rng = np.random.default_rng(seed)
    sigma = sigma_from_relative_level(b_clean, rel_level)
    e = sigma * rng.standard_normal(b_clean.size)
    b = b_clean + e
    return b, _info("white", rel_level, sigma, b_clean, e, seed, rho=0.0)


def ar1_correlated(b_clean, rel_level, rho, seed, burn_in=200):
    """Stationary AR(1) perturbation: e_t = rho*e_{t-1} + sqrt(1-rho^2)*sigma*z_t.

    The process is started from its stationary distribution and a burn-in is run
    anyway, so the realization is stationary along the whole index range. Marginal
    std is sigma, identical to `white_gaussian` at the same rel_level; the
    difference is entirely in the correlation structure. rho in (-1, 1).

    Why it matters here: a large positive rho concentrates the noise power at low
    frequencies, i.e. in the directions of the LARGE singular vectors of a smoothing
    operator, which is exactly where rules that assume white noise (GCV, and the
    plain discrepancy principle with a scalar sigma) misjudge the trade-off.
    """
    if not -1.0 < rho < 1.0:
        raise ValueError("ar1_correlated: need -1 < rho < 1")
    b_clean = np.asarray(b_clean, dtype=np.float64)
    m = b_clean.size
    rng = np.random.default_rng(seed)
    sigma = sigma_from_relative_level(b_clean, rel_level)

    z = rng.standard_normal(m + burn_in)
    e = np.empty(m + burn_in)
    e[0] = z[0]                      # stationary start: Var = 1
    c = np.sqrt(1.0 - rho ** 2)
    for i in range(1, m + burn_in):
        e[i] = rho * e[i - 1] + c * z[i]
    e = sigma * e[burn_in:]
    b = b_clean + e
    return b, _info("ar1", rel_level, sigma, b_clean, e, seed, rho=rho)


def ar1_covariance(m, rho, sigma=1.0):
    """The exact stationary AR(1) covariance matrix sigma^2 * rho^{|i-j|}.

    Useful for whitening (Cholesky) or for building the generalized discrepancy
    functional when a method is allowed to know the true noise covariance.
    """
    idx = np.arange(m)
    return (sigma ** 2) * rho ** np.abs(idx[:, None] - idx[None, :])


def misspecified_sigma(sigma_true, factor):
    """The sigma a rule is told, = factor * sigma_true. factor 0.5 = told half the
    true level (rule will under-regularize); factor 2.0 = told double (will
    over-regularize). The data is never regenerated."""
    return float(factor) * float(sigma_true)


def add_noise(b_clean, kind, rel_level, seed, rho=0.0):
    """Dispatcher. kind in {'white', 'ar1'}. Returns (b_noisy, info)."""
    if kind == "white":
        return white_gaussian(b_clean, rel_level, seed)
    if kind == "ar1":
        return ar1_correlated(b_clean, rel_level, rho, seed)
    raise ValueError(f"add_noise: unknown kind {kind!r}")


def relative_error(x_hat, x_true):
    """||x_hat - x_true||_2 / ||x_true||_2."""
    x_hat = np.asarray(x_hat, dtype=np.float64)
    x_true = np.asarray(x_true, dtype=np.float64)
    return float(np.linalg.norm(x_hat - x_true) / np.linalg.norm(x_true))


def _info(kind, rel_level, sigma, b_clean, e, seed, rho):
    return {
        "kind": kind,
        "rel_level_nominal": float(rel_level),
        "sigma_true": float(sigma),
        "rho": float(rho),
        "seed": int(seed),
        "noise_norm": float(np.linalg.norm(e)),
        "rel_error_realized": float(np.linalg.norm(e) / np.linalg.norm(b_clean)),
        "discrepancy_target": float(sigma * np.sqrt(e.size)),
        "level_convention": "sigma = rel_level * ||b_clean||_2 / sqrt(m)",
    }


if __name__ == "__main__":  # pragma: no cover
    # tiny self-check: marginal std and lag-1 correlation of the AR(1) model
    b = np.ones(4000)
    for rho in (0.0, 0.5, 0.9):
        e_b, info = ar1_correlated(b, 0.01, rho, seed=1) if rho else white_gaussian(
            b, 0.01, seed=1
        )
        e = e_b - b
        r1 = float(np.corrcoef(e[:-1], e[1:])[0, 1])
        print(
            f"rho={rho:.1f} std={e.std():.6f} target={info['sigma_true']:.6f} "
            f"lag1_corr={r1:.3f} rel_err={info['rel_error_realized']:.5f}"
        )
