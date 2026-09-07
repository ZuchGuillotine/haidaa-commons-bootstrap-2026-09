"""
inv_critique_probe.py -- skeptical-reviewer probe of inv_paramchoice.py.

Answers four questions about the parameter-choice failure map, by re-running the
same cells with the knobs that the original script fixed:

A) Lambda-grid floor and SVD truncation.  The original uses
   lambdas = logspace(log10(s1), log10(s_min) - 2, 300) after discarding
   singular values below 1e-14*s1.  How much of the reported catastrophe
   MAGNITUDE (median error ratios of 1e9-1e12) is set by those two constants?
   Variants: (trunc, extra decades below s_min).

B) L-curve corner detector.  The original maximizes a centred-finite-difference
   curvature of the discrete log-log L-curve.  Hansen's reference
   implementations use the analytic curvature (closed form in lambda) or a
   spline-smoothed curve.  Does the analytic curvature change the verdict on
   deconv_smooth?

C) Failure criterion.  "relerr > 2 * oracle relerr" is scale-free, so a cell
   with a tiny oracle error (deconv_smooth: 0.0146) marks a 0.09 reconstruction
   as a failure.  Report absolute errors and an alternative criterion.

D) Discrepancy principle tau.  The original fixes tau = 1 and varies the
   assumed noise level by 0.5/1/2.  Because the target is tau*fac*delta, the
   "discrepancy_double" column is algebraically identical to the textbook
   tau = 2 rule with a CORRECT noise level.  Sweep tau explicitly.

Usage:
    python inv_critique_probe.py [output.json]
"""

import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import inv_paramchoice as pc  # noqa: E402


CELLS = [(f, nk, rho, lv)
         for f in ("shaw", "deconv_gauss", "deconv_smooth")
         for lv in (0.01, 0.05)
         for nk, rho in (("white", 0.0), ("ar1", 0.5), ("ar1", 0.9))]


def analytic_lcurve_corner(s, beta, lambdas, b_out2=0.0):
    """Exact curvature of the log-log L-curve, differentiated analytically.

    Same geometric object the record's detector approximates by centred finite
    differences: the curve (xi, eta_l) = (log ||A x_lam - b||, log ||x_lam||)
    parameterized by t = log lambda.  Here d/dt is exact.

    With f_i = s_i^2/(s_i^2+lam^2), df/dt = -2 f (1-f),
      R  = sum (1-f)^2 beta^2 + b_out2 ,  E = sum f^2 (beta/s)^2
      R' = 4 sum f (1-f)^2 beta^2       ,  E' = -4 sum f^2 (1-f) (beta/s)^2
      R''= -8 sum f (1-f)^2 (1-3f) beta^2, E''= 8 sum f^2 (1-f)(2-3f)(beta/s)^2
    and xi = 0.5 log R, eta_l = 0.5 log E.
    """
    f = (s ** 2)[None, :] / ((s ** 2)[None, :] + (lambdas ** 2)[:, None])
    b2 = (beta ** 2)[None, :]
    q2 = ((beta / s) ** 2)[None, :]
    R = ((1.0 - f) ** 2 * b2).sum(axis=1) + b_out2
    E = (f ** 2 * q2).sum(axis=1)
    Rp = 4.0 * (f * (1.0 - f) ** 2 * b2).sum(axis=1)
    Ep = -4.0 * (f ** 2 * (1.0 - f) * q2).sum(axis=1)
    Rpp = -8.0 * (f * (1.0 - f) ** 2 * (1.0 - 3.0 * f) * b2).sum(axis=1)
    Epp = 8.0 * (f ** 2 * (1.0 - f) * (2.0 - 3.0 * f) * q2).sum(axis=1)
    R = np.maximum(R, 1e-300)
    E = np.maximum(E, 1e-300)
    d1x = 0.5 * Rp / R
    d1y = 0.5 * Ep / E
    d2x = 0.5 * (Rpp * R - Rp ** 2) / R ** 2
    d2y = 0.5 * (Epp * E - Ep ** 2) / E ** 2
    denom = np.maximum((d1x ** 2 + d1y ** 2) ** 1.5, 1e-300)
    return (d1x * d2y - d2x * d1y) / denom


def run(fixture, noise_kind, rho, rel_level, trunc, extra_decades,
        n_seeds=20, n_lambda=300):
    A, x_true, b_clean = pc.get_fixture(fixture)
    m, n = A.shape
    U, s, Vt = np.linalg.svd(A, full_matrices=False)
    keep = s > s[0] * trunc
    U, s, Vt = U[:, keep], s[keep], Vt[keep, :]
    lambdas = np.sort(np.logspace(np.log10(s[0]),
                                  np.log10(max(s[-1], s[0] * 1e-12)) - extra_decades,
                                  n_lambda))
    xt = np.linalg.norm(x_true)
    acc = {k: [] for k in ("gcv_ratio", "gcv_fail", "gcv_relerr",
                           "dhalf_ratio", "dhalf_fail",
                           "lc_fd_ratio", "lc_fd_fail", "lc_fd_relerr",
                           "lc_an_ratio", "lc_an_fail", "lc_an_relerr",
                           "oracle_relerr", "lc_fd_nmax", "lc_an_nmax")}
    for tau in (1.0, 1.01, 1.1, 1.5, 2.0):
        acc["tau%g_ratio" % tau] = []
        acc["tau%g_fail" % tau] = []
    for si in range(n_seeds):
        import hashlib
        key = "%s|%s|%.3f|%.4f|%d" % (fixture, noise_kind, rho, rel_level, si)
        seed = int(hashlib.sha256(key.encode()).hexdigest()[:8], 16)
        rng = np.random.default_rng(seed)
        e = pc.make_noise(m, noise_kind, rho, rel_level, b_clean, rng)
        b = b_clean + e
        delta = float(np.linalg.norm(e))
        X, resid, xnorm, trace, _ = pc.tikh_curves(U, s, Vt, b, lambdas)
        errs = np.linalg.norm(X - x_true[None, :], axis=1) / xt
        e_or = float(errs.min())
        acc["oracle_relerr"].append(e_or)

        ig, _ = pc.pick_gcv(resid, trace, m)
        acc["gcv_ratio"].append(errs[ig] / e_or)
        acc["gcv_fail"].append(errs[ig] > 2 * e_or)
        acc["gcv_relerr"].append(float(errs[ig]))

        ih, _ = pc.pick_discrepancy(resid, lambdas, 0.5 * delta)
        acc["dhalf_ratio"].append(errs[ih] / e_or)
        acc["dhalf_fail"].append(errs[ih] > 2 * e_or)

        il, kap_fd, diag = pc.pick_lcurve(resid, xnorm, lambdas)
        acc["lc_fd_ratio"].append(errs[il] / e_or)
        acc["lc_fd_fail"].append(errs[il] > 2 * e_or)
        acc["lc_fd_relerr"].append(float(errs[il]))
        acc["lc_fd_nmax"].append(diag["n_local_curvature_maxima"])

        beta = U.T @ b
        b_out2 = max(float(b @ b - beta @ beta), 0.0)
        kap = analytic_lcurve_corner(s, beta, lambdas, b_out2)
        kk = np.full_like(kap, -np.inf)
        kk[2:-2] = kap[2:-2]
        ia = int(np.argmax(kk))
        nmax = int(sum(1 for i in range(3, kap.size - 3)
                       if kap[i] > kap[i - 1] and kap[i] >= kap[i + 1] and kap[i] > 0))
        acc["lc_an_ratio"].append(errs[ia] / e_or)
        acc["lc_an_fail"].append(errs[ia] > 2 * e_or)
        acc["lc_an_relerr"].append(float(errs[ia]))
        acc["lc_an_nmax"].append(nmax)

        for tau in (1.0, 1.01, 1.1, 1.5, 2.0):
            it, _ = pc.pick_discrepancy(resid, lambdas, delta, tau=tau)
            acc["tau%g_ratio" % tau].append(errs[it] / e_or)
            acc["tau%g_fail" % tau].append(errs[it] > 2 * e_or)

    out = {"fixture": fixture, "noise": noise_kind, "rho": rho,
           "rel_noise_level": rel_level, "svd_trunc": trunc,
           "extra_decades": extra_decades,
           "lambda_min": float(lambdas[0]), "lambda_max": float(lambdas[-1]),
           "numerical_rank": int(s.size), "cond_kept": float(s[0] / s[-1])}
    for k, v in acc.items():
        v = np.asarray(v, float)
        if k.endswith("_fail"):
            out[k + "_rate"] = float(v.mean())
        else:
            out["median_" + k] = float(np.median(v))
            if k.endswith("_ratio"):
                out["max_" + k] = float(v.max())
    return out


def main():
    outp = sys.argv[1] if len(sys.argv) > 1 else "inv_critique_probe_output.json"
    variants = [(1e-14, 2.0), (1e-14, 0.0), (1e-10, 2.0), (1e-8, 0.0)]
    res = []
    for trunc, dec in variants:
        for (f, nk, rho, lv) in CELLS:
            r = run(f, nk, rho, lv, trunc, dec)
            res.append(r)
    doc = {"script": "inv_critique_probe.py",
           "role": "skeptical-reviewer",
           "variants_trunc_extradecades": variants,
           "baseline_variant": [1e-14, 2.0],
           "cells": res}
    with open(outp, "w") as fh:
        json.dump(doc, fh, indent=1, sort_keys=True)

    # ---- printed summaries -------------------------------------------------
    def sel(trunc, dec):
        return [r for r in res if r["svd_trunc"] == trunc and r["extra_decades"] == dec]

    print("== A) grid floor / SVD truncation sensitivity (pooled over 18 cells x 20 seeds)")
    print("trunc   dec  gcv_failrate gcv_med_ratio(worst cell)  dhalf_failrate dhalf_med_ratio(worst)")
    for trunc, dec in variants:
        rs = sel(trunc, dec)
        g = np.mean([r["gcv_fail_rate"] for r in rs])
        h = np.mean([r["dhalf_fail_rate"] for r in rs])
        gm = max(r["median_gcv_ratio"] for r in rs)
        hm = max(r["median_dhalf_ratio"] for r in rs)
        print(f"{trunc:.0e} {dec:4.1f} {g:12.4f} {gm:24.3e}  {h:14.4f} {hm:20.3e}")

    print("\n== A2) per-cell GCV median ratio, baseline vs trunc=1e-8/floor=s_min")
    base = {(r['fixture'], r['noise'], r['rho'], r['rel_noise_level']): r for r in sel(1e-14, 2.0)}
    alt = {(r['fixture'], r['noise'], r['rho'], r['rel_noise_level']): r for r in sel(1e-8, 0.0)}
    for k in base:
        b_, a_ = base[k], alt[k]
        print(f"{k[0]:13s} {k[1]}{k[2]} lv={k[3]:.2f} | gcv med ratio {b_['median_gcv_ratio']:.3e}"
              f" -> {a_['median_gcv_ratio']:.3e} | failrate {b_['gcv_fail_rate']:.2f} -> {a_['gcv_fail_rate']:.2f}"
              f" | dhalf med {b_['median_dhalf_ratio']:.2e} -> {a_['median_dhalf_ratio']:.2e}")

    print("\n== B) L-curve corner: finite-difference (record) vs analytic curvature (baseline variant)")
    for r in sel(1e-14, 2.0):
        print(f"{r['fixture']:13s} {r['noise']}{r['rho']} lv={r['rel_noise_level']:.2f} |"
              f" FD fail {r['lc_fd_fail_rate']:.2f} med {r['median_lc_fd_ratio']:6.2f} nmax {r['median_lc_fd_nmax']:.0f} |"
              f" AN fail {r['lc_an_fail_rate']:.2f} med {r['median_lc_an_ratio']:6.2f} nmax {r['median_lc_an_nmax']:.0f}")

    print("\n== C) absolute errors behind the '2x oracle' failures (baseline variant)")
    for r in sel(1e-14, 2.0):
        print(f"{r['fixture']:13s} {r['noise']}{r['rho']} lv={r['rel_noise_level']:.2f} |"
              f" oracle relerr {r['median_oracle_relerr']:.4f} |"
              f" Lcurve relerr {r['median_lc_fd_relerr']:.4f} (fail {r['lc_fd_fail_rate']:.2f}) |"
              f" GCV relerr {r['median_gcv_relerr']:.4g}")

    print("\n== D) discrepancy principle: tau sweep with the CORRECT noise norm (baseline variant)")
    rs = sel(1e-14, 2.0)
    for tau in (1.0, 1.01, 1.1, 1.5, 2.0):
        fr = np.mean([r["tau%g_fail_rate" % tau] for r in rs])
        med = np.median([r["median_tau%g_ratio" % tau] for r in rs])
        worst = max(r["max_tau%g_ratio" % tau] for r in rs)
        print(f"tau={tau:<5} pooled failure rate {fr:.4f}  median ratio {med:.3f}  worst instance {worst:.2f}")
    print("\nwrote", outp)


if __name__ == "__main__":
    main()
