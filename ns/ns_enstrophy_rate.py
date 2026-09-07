#!/usr/bin/env python3
"""
ns_enstrophy_rate.py
--------------------
Reproducible fixture: instantaneous enstrophy growth rate for divergence-free
velocity fields on the periodic box [0, 2*pi]^3, compared with the
Lu-Doering-type a-priori bound  dE/dt <= C * E^3 / nu^3.

Conventions (IMPORTANT -- the numbers depend on them):
  omega = curl u
  E     = (1/2) * INTEGRAL_{[0,2pi]^3} |omega|^2 dV          (enstrophy)
  P     = INTEGRAL |curl omega|^2 dV                          (palinstrophy)
  S     = INTEGRAL omega . (omega . grad) u dV                (vortex stretching)
  dE/dt = -nu * P + S
Volume INTEGRALS, not means. E has dimensions L^3 T^-2, so
  R := (dE/dt) * nu^3 / E^3     and     R_S := S * nu^3 / E^3
are dimensionless and invariant under the Navier-Stokes scaling group.

Two further reduced observables are used throughout, with L := 2*pi:
  Q     := S * L^(3/2) / E^(3/2)      dimensionless "stretching shape number",
                                      invariant under u -> lambda*u and under nu
  Re_E  := sqrt(E * L) / nu           enstrophy-based Reynolds number
and they satisfy the exact algebraic identity
  R_S = Q * Re_E^(-3).
So the distance from the bound dE/dt <= C E^3/nu^3 is controlled entirely by the
shape number Q and the Reynolds number; a given field family can only approach
the bound at Re_E = O(1).  This is an identity, not an empirical fit.
If you instead define Etilde = INTEGRAL |omega|^2 (no 1/2), then
Rtilde = R / 4.  If you use domain means instead of integrals on [0,2pi]^3,
multiply E by (2*pi)^-3 and R changes by (2*pi)^6.

The constant C of Lu & Doering (2008) is NOT hard-coded here: it was not
verified against the primary source while this fixture was written.  R and R_S
are reported as raw observables.

Numerics
  - Fourier pseudo-spectral, real FFTs (numpy.fft.rfftn / irfftn).
  - Every test field is band-limited to |k_i| <= KMAX with 3*KMAX < N_eval,
    so the grid quadrature of the CUBIC integrand of S is exact (no aliasing
    error in S).  Static cases evaluate S on a padded grid of size PAD*N.
  - Quadratic integrals (E, P, energy) are exact for 2*KMAX < N.
  - All random fields use numpy.random.default_rng(seed) with fixed seeds.

Usage:
  python ns_enstrophy_rate.py                 # full fixture, prints JSON
  python ns_enstrophy_rate.py --quick         # skip the time integration
  python ns_enstrophy_rate.py --out FILE.json # also write JSON to FILE
"""

import argparse
import json
import platform
import sys
import time

import numpy as np
import scipy

# ----------------------------------------------------------------------------
# spectral machinery
# ----------------------------------------------------------------------------


def wavenumbers(N):
    """Integer wavenumber arrays for a real FFT of an N^3 field on [0,2pi]^3."""
    k1 = np.fft.fftfreq(N, d=1.0 / N)
    k3 = np.fft.rfftfreq(N, d=1.0 / N)
    KX = k1[:, None, None] * np.ones((1, N, 1)) * np.ones((1, 1, k3.size))
    KY = k1[None, :, None] * np.ones((N, 1, 1)) * np.ones((1, 1, k3.size))
    KZ = k3[None, None, :] * np.ones((N, 1, 1)) * np.ones((1, N, 1))
    K = np.stack([KX, KY, KZ]).astype(np.float64)
    K2 = (K ** 2).sum(axis=0)
    K2inv = np.where(K2 > 0, 1.0 / np.where(K2 > 0, K2, 1.0), 0.0)
    return K, K2, K2inv


def fft3(u):
    return np.fft.rfftn(u, axes=(0, 1, 2))


def ifft3(uh, N):
    return np.fft.irfftn(uh, s=(N, N, N), axes=(0, 1, 2))


def vfft3(u):
    return np.stack([fft3(u[i]) for i in range(3)])


def vifft3(uh, N):
    return np.stack([ifft3(uh[i], N) for i in range(3)])


def leray(uh, K, K2inv):
    """Project onto divergence-free fields:  uh - k (k.uh)/|k|^2."""
    kdotu = (K * uh).sum(axis=0)
    return uh - K * (kdotu * K2inv)


def curl_hat(uh, K):
    return 1j * np.stack([
        K[1] * uh[2] - K[2] * uh[1],
        K[2] * uh[0] - K[0] * uh[2],
        K[0] * uh[1] - K[1] * uh[0],
    ])


def bandlimit_mask(K, kmax):
    return (np.abs(K[0]) <= kmax) & (np.abs(K[1]) <= kmax) & (np.abs(K[2]) <= kmax)


def dV(N):
    return (2.0 * np.pi / N) ** 3


def pad_hat(uh, N, Np):
    """Zero-pad a vector spectral field from grid N to grid Np (Np > N).

    The physical field is unchanged; only the sampling grid is refined.
    """
    kn = N // 2
    out = np.zeros((3, Np, Np, Np // 2 + 1), dtype=complex)
    f1 = np.fft.fftfreq(N, d=1.0 / N).astype(int)
    f1p = np.fft.fftfreq(Np, d=1.0 / Np).astype(int)
    idx = np.array([np.where(f1p == kk)[0][0] for kk in f1])
    nz = min(N // 2 + 1, Np // 2 + 1)
    sub = uh[:, :, :, :nz]
    out[np.ix_(np.arange(3), idx, idx, np.arange(nz))] = sub
    # irfftn on the padded grid needs the amplitude convention rescaled
    return out * (float(Np) / float(N)) ** 3, kn


# ----------------------------------------------------------------------------
# diagnostics
# ----------------------------------------------------------------------------


def diagnostics(uh, N, nu, pad=2):
    """Return E, P, S, dE/dt, R, R_S and auxiliary quantities for velocity uh."""
    K, K2, K2inv = wavenumbers(N)
    wh = curl_hat(uh, K)
    ch = curl_hat(wh, K)  # curl omega

    u = vifft3(uh, N)
    w = vifft3(wh, N)
    c = vifft3(ch, N)
    dv = dV(N)

    energy = 0.5 * float((u ** 2).sum()) * dv
    E = 0.5 * float((w ** 2).sum()) * dv
    P = float((c ** 2).sum()) * dv

    # divergence check
    divh = 1j * (K * uh).sum(axis=0)
    div = ifft3(divh, N)
    div_l2 = float(np.sqrt((div ** 2).sum() * dv))
    u_l2 = float(np.sqrt((u ** 2).sum() * dv))

    # vortex stretching on a padded grid so the cubic quadrature is exact
    Np = int(round(pad * N))
    if Np % 2:
        Np += 1
    uhp, _ = pad_hat(uh, N, Np)
    Kp, _, _ = wavenumbers(Np)
    whp = curl_hat(uhp, Kp)
    wp = vifft3(whp, Np)
    S = 0.0
    for i in range(3):
        for j in range(3):
            dj_ui = ifft3(1j * Kp[j] * uhp[i], Np)
            S += float((wp[i] * wp[j] * dj_ui).sum())
    S *= dV(Np)

    dEdt = -nu * P + S

    # resolution indicator: fraction of enstrophy above 0.8*kmax_resolved
    kmag = np.sqrt(K2)
    wgt = np.ones_like(kmag)
    wgt[..., 1:] = 2.0  # rfft: modes with kz>0 counted once, need conjugates
    if N % 2 == 0:
        wgt[..., -1] = 1.0
    ens_spec = 0.5 * wgt * (np.abs(wh) ** 2).sum(axis=0) / (N ** 6) * (2 * np.pi) ** 3
    tot = float(ens_spec.sum())
    khi = 0.8 * (N // 3)
    tail = float(ens_spec[kmag > khi].sum())
    tail_frac = tail / tot if tot > 0 else 0.0

    urms = float(np.sqrt((u ** 2).sum() / (3.0 * N ** 3)))
    L = 2.0 * np.pi
    Re = urms * L / nu
    Q = S * L ** 1.5 / E ** 1.5 if E > 0 else float("nan")
    Re_E = np.sqrt(E * L) / nu if E > 0 else float("nan")

    return dict(
        N=N, nu=nu, energy=energy, enstrophy_E=E, palinstrophy_P=P,
        viscous_term=-nu * P, stretching_S=S, dEdt=dEdt,
        R_total=dEdt * nu ** 3 / E ** 3 if E > 0 else float("nan"),
        R_stretching=S * nu ** 3 / E ** 3 if E > 0 else float("nan"),
        u_rms=urms, Re_box=Re, Q_shape=Q, Re_E=float(Re_E),
        identity_residual=abs(Q * Re_E ** -3.0 - S * nu ** 3 / E ** 3) /
        max(abs(S * nu ** 3 / E ** 3), 1e-300) if E > 0 else float("nan"),
        div_rel_l2=div_l2 / u_l2 if u_l2 > 0 else 0.0,
        enstrophy_tail_frac=tail_frac,
        spectral_E_check=tot / E - 1.0 if E > 0 else float("nan"),
    )


# ----------------------------------------------------------------------------
# test fields
# ----------------------------------------------------------------------------


def taylor_green(N, A=1.0):
    x = np.linspace(0, 2 * np.pi, N, endpoint=False)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    u = np.stack([
        A * np.sin(X) * np.cos(Y) * np.cos(Z),
        -A * np.cos(X) * np.sin(Y) * np.cos(Z),
        np.zeros_like(X),
    ])
    return vfft3(u)


def antiparallel_tubes(N, kmax, Gamma=1.0, sep=0.9, core=0.35,
                       wy=0.15, wz=0.0, phase=0.0, asym=1.0):
    """Two antiparallel Gaussian vortex tubes along x, sinusoidally perturbed.

    The idealized vorticity (purely x-aligned Gaussian tubes with wavy centre
    lines) is Leray-projected and band-limited, so the field actually used is
    divergence-free to machine precision and spectrally smooth, and the
    velocity follows by Biot-Savart, u_hat = i k x omega_hat / |k|^2.

    Defaults (wz=0, phase=0, asym=1) give the mirror-symmetric configuration:
    omega_1, omega_2 even and omega_3 odd about the plane z = pi, which forces
    the stretching integral S to vanish identically.  Set wz > 0 with
    phase != 0 and asym != 1 to break that symmetry.
    """
    x = np.linspace(0, 2 * np.pi, N, endpoint=False)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    yc, zc = np.pi, np.pi
    y1 = yc + sep + wy * np.cos(X)
    z1 = zc + wz * np.sin(X)
    y2 = yc - sep - asym * wy * np.cos(X + phase)
    z2 = zc - asym * wz * np.sin(X + phase)

    def g(dy, dz):
        # periodized Gaussian (nearest image is enough for core << pi)
        dy = (dy + np.pi) % (2 * np.pi) - np.pi
        dz = (dz + np.pi) % (2 * np.pi) - np.pi
        return np.exp(-(dy ** 2 + dz ** 2) / (2 * core ** 2))

    wx = Gamma * (g(Y - y1, Z - z1) - g(Y - y2, Z - z2))
    w = np.stack([wx, np.zeros_like(wx), np.zeros_like(wx)])
    wh = vfft3(w)
    K, K2, K2inv = wavenumbers(N)
    wh = leray(wh, K, K2inv)
    wh *= bandlimit_mask(K, kmax)
    wh[:, 0, 0, 0] = 0.0
    uh = curl_hat(wh, K) * K2inv  # u_hat = i k x omega_hat / |k|^2
    return uh


def random_field(N, kmax, slope, seed, kf_lo=1.0):
    """Divergence-free random field with 3D energy spectrum ~ k^{-slope}."""
    rng = np.random.default_rng(seed)
    noise = rng.standard_normal((3, N, N, N))
    uh = vfft3(noise)
    K, K2, K2inv = wavenumbers(N)
    kmag = np.sqrt(K2)
    amp = np.zeros_like(kmag)
    good = kmag >= kf_lo
    amp[good] = kmag[good] ** (-(slope + 2.0) / 2.0)
    uh = uh * amp
    uh = leray(uh, K, K2inv)
    uh *= bandlimit_mask(K, kmax)
    uh[:, 0, 0, 0] = 0.0
    return uh


def rescale_to_enstrophy(uh, N, target_E):
    K, _, _ = wavenumbers(N)
    w = vifft3(curl_hat(uh, K), N)
    E = 0.5 * float((w ** 2).sum()) * dV(N)
    return uh * np.sqrt(target_E / E)


# ----------------------------------------------------------------------------
# time integration (RK4, 2/3 dealiasing, integrating factor for viscosity)
# ----------------------------------------------------------------------------


def rhs(uh, N, nu, K, K2, K2inv, mask):
    u = vifft3(uh, N)
    wh = curl_hat(uh, K)
    w = vifft3(wh, N)
    # rotation form: u x omega
    nl = np.stack([
        u[1] * w[2] - u[2] * w[1],
        u[2] * w[0] - u[0] * w[2],
        u[0] * w[1] - u[1] * w[0],
    ])
    nlh = vfft3(nl) * mask
    nlh = leray(nlh, K, K2inv)
    return nlh - nu * K2 * uh


def integrate(N, nu, uh0, t_end, dt, sample_every, pad=2):
    K, K2, K2inv = wavenumbers(N)
    mask = bandlimit_mask(K, N // 3)
    uh = uh0 * mask
    uh = leray(uh, K, K2inv)
    t = 0.0
    out = []
    step = 0
    nsteps = int(round(t_end / dt))
    while step <= nsteps:
        if step % sample_every == 0:
            d = diagnostics(uh, N, nu, pad=pad)
            out.append(dict(t=round(t, 6), E=d["enstrophy_E"], energy=d["energy"],
                            S=d["stretching_S"], dEdt=d["dEdt"],
                            Q_shape=d["Q_shape"], Re_E=d["Re_E"],
                            R_stretching=d["R_stretching"], R_total=d["R_total"],
                            tail_frac=d["enstrophy_tail_frac"]))
        k1 = rhs(uh, N, nu, K, K2, K2inv, mask)
        k2 = rhs(uh + 0.5 * dt * k1, N, nu, K, K2, K2inv, mask)
        k3 = rhs(uh + 0.5 * dt * k2, N, nu, K, K2, K2inv, mask)
        k4 = rhs(uh + dt * k3, N, nu, K, K2, K2inv, mask)
        uh = uh + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        uh *= mask
        t += dt
        step += 1
    return out


# ----------------------------------------------------------------------------
# main
# ----------------------------------------------------------------------------


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    t0 = time.time()
    res = {
        "fixture": "ns_enstrophy_rate",
        "version": "1.0",
        "conventions": {
            "domain": "[0,2pi]^3 periodic",
            "E": "0.5 * volume-integral of |curl u|^2",
            "P": "volume-integral of |curl curl u|^2",
            "S": "volume-integral of omega_i omega_j d_j u_i",
            "dEdt": "-nu*P + S",
            "R_total": "dEdt * nu^3 / E^3 (dimensionless)",
            "R_stretching": "S * nu^3 / E^3 (dimensionless)",
            "Q_shape": "S * L^1.5 / E^1.5 with L=2pi; amplitude- and nu-independent",
            "Re_E": "sqrt(E*L)/nu with L=2pi",
            "identity": "R_stretching = Q_shape * Re_E^-3 exactly (algebra, not a fit)",
            "Re_box": "u_rms * 2pi / nu",
            "bound": "dE/dt <= C * E^3 / nu^3 (Lu & Doering 2008); C NOT verified here",
        },
        "environment": {
            "python": sys.version.split()[0],
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "platform": platform.platform(),
            "machine": platform.machine(),
        },
        "cases": {},
    }

    N = 48
    KMAX = N // 3 - 1  # 15: guarantees 3*KMAX < 2*N so cubic quadrature is exact
    res["grid"] = {"N": N, "kmax_bandlimit": KMAX, "pad_factor_for_cubic": 2}

    # ---- 1. analytic verification on the Taylor-Green vortex -----------------
    # For u = A(sin x cos y cos z, -cos x sin y cos z, 0):
    #   E(0)      = 3 A^2 pi^3
    #   P(0)      = 6 E(0)      (every Fourier mode of omega has |k|^2 = 3)
    #   S(0)      = 0           (every term carries an odd power of cos x)
    ver = {}
    for A in (0.5, 1.0, 2.0):
        uh = taylor_green(N, A)
        d = diagnostics(uh, N, nu=0.01, pad=2)
        E_exact = 3.0 * A ** 2 * np.pi ** 3
        P_exact = 6.0 * E_exact
        ver["A=%g" % A] = dict(
            E=d["enstrophy_E"], E_exact=E_exact,
            E_rel_err=abs(d["enstrophy_E"] - E_exact) / E_exact,
            P=d["palinstrophy_P"], P_exact=P_exact,
            P_rel_err=abs(d["palinstrophy_P"] - P_exact) / P_exact,
            S=d["stretching_S"], S_exact=0.0,
            S_abs_err_scaled=abs(d["stretching_S"]) / (E_exact ** 1.5),
            div_rel_l2=d["div_rel_l2"],
        )
    res["cases"]["analytic_check_taylor_green_t0"] = ver

    # ---- 2. Taylor-Green at several amplitudes and viscosities ---------------
    tg = []
    for A in (0.5, 1.0, 2.0):
        for nu in (0.1, 0.02, 0.005):
            uh = taylor_green(N, A)
            d = diagnostics(uh, N, nu, pad=2)
            d["case"] = "taylor_green_t0"
            d["A"] = A
            tg.append(d)
    res["cases"]["taylor_green_t0"] = tg

    # ---- 3. antiparallel vortex tubes ---------------------------------------
    # 3a. mirror-symmetric pair: S vanishes identically by the z-reflection
    #     symmetry (omega_1, omega_2 even and omega_3 odd about z = pi).
    #     This is a further exact check on the stretching operator.
    tubes = []
    uh_sym = antiparallel_tubes(N, KMAX, wy=0.15, wz=0.0, phase=0.0, asym=1.0)
    d_sym = diagnostics(uh_sym, N, 0.02, pad=2)
    d_sym["case"] = "antiparallel_tubes_mirror_symmetric"
    res["cases"]["antiparallel_tubes_mirror_symmetric"] = d_sym

    # 3b. symmetry-broken pair: the two tubes are displaced in z with a phase
    #     offset and unequal amplitudes, so no discrete symmetry forces S = 0.
    uh_t = antiparallel_tubes(N, KMAX, wy=0.15, wz=0.15, phase=0.7, asym=0.6)
    for nu in (0.1, 0.02, 0.005):
        d = diagnostics(uh_t, N, nu, pad=2)
        d["case"] = "antiparallel_tubes_perturbed"
        tubes.append(d)
    # resolution check on the same physical field at N=64
    uh_t64 = antiparallel_tubes(64, KMAX, wy=0.15, wz=0.15, phase=0.7, asym=0.6)
    d64 = diagnostics(uh_t64, 64, 0.02, pad=2)
    res["cases"]["antiparallel_tubes_perturbed"] = tubes
    res["cases"]["antiparallel_tubes_resolution_check_N64"] = d64
    res["cases"]["antiparallel_tubes_N48_vs_N64"] = {
        "Q_N48": tubes[1]["Q_shape"], "Q_N64": d64["Q_shape"],
        "rel_diff": abs(tubes[1]["Q_shape"] - d64["Q_shape"]) /
        abs(d64["Q_shape"]),
    }

    # ---- 4. random divergence-free fields with prescribed spectra -----------
    rnd = []
    for slope in (5.0 / 3.0, 3.0, 1.0):
        for seed in (11, 12, 13):
            uh = random_field(N, KMAX, slope, seed)
            uh = rescale_to_enstrophy(uh, N, target_E=1.0)
            for nu in (0.1, 0.02, 0.005):
                d = diagnostics(uh, N, nu, pad=2)
                d["case"] = "random_spectrum"
                d["slope"] = slope
                d["seed"] = seed
                rnd.append(d)
    res["cases"]["random_fields"] = rnd

    # ---- 5. invariance / identity sweep -------------------------------------
    # R_S = Q * Re_E^-3 is an identity, NOT an empirical scaling law.  What is
    # worth checking numerically is that Q is genuinely invariant under
    # amplitude rescaling and under changes of nu, and that the identity holds
    # to round-off in the code.  Then the whole nu- and amplitude-dependence of
    # the distance to the bound is carried by Re_E alone.
    uh_ref = rescale_to_enstrophy(random_field(N, KMAX, 5.0 / 3.0, 101), N, 1.0)
    sweep = []
    for lam in (0.25, 1.0, 4.0):
        for nu in (1.0, 0.1, 0.01, 0.001):
            d = diagnostics(uh_ref * lam, N, nu, pad=2)
            sweep.append(dict(amplitude_lambda=lam, nu=nu, E=d["enstrophy_E"],
                              S=d["stretching_S"], Q_shape=d["Q_shape"],
                              Re_E=d["Re_E"], R_stretching=d["R_stretching"],
                              R_total=d["R_total"],
                              identity_residual=d["identity_residual"]))
    Qs = np.array([s["Q_shape"] for s in sweep])
    res["cases"]["invariance_sweep"] = {
        "field": "random_spectrum slope=5/3 seed=101, enstrophy normalized to E=1",
        "points": sweep,
        "Q_mean": float(Qs.mean()),
        "Q_max_rel_spread": float(np.abs(Qs / Qs.mean() - 1.0).max()),
        "max_identity_residual": max(s["identity_residual"] for s in sweep),
        "note": "Q is invariant under u -> lambda u and under nu; R_S = Q*Re_E^-3",
    }

    # ---- 6. short time integration of the Taylor-Green vortex ---------------
    # The t=0 fields above are all kinematic.  Real vortex stretching is
    # produced dynamically (strain-vorticity alignment), so the largest Q in
    # this fixture should come from an evolved field, not from a chosen one.
    tgi = None
    if not args.quick:
        nu_t = 0.025
        uh0 = taylor_green(N, 1.0)
        traj = integrate(N, nu_t, uh0, t_end=10.0, dt=0.02, sample_every=10)
        Qs = [p["Q_shape"] for p in traj]
        tgi = {
            "N": N, "nu": nu_t, "dt": 0.02, "t_end": 10.0,
            "scheme": "RK4, rotation-form nonlinearity, 2/3 dealiasing",
            "trajectory": traj,
            "max_Q_shape": max(Qs),
            "t_at_max_Q_shape": traj[int(np.argmax(Qs))]["t"],
            "Re_E_at_max_Q": traj[int(np.argmax(Qs))]["Re_E"],
            "R_stretching_at_max_Q": traj[int(np.argmax(Qs))]["R_stretching"],
            "max_R_stretching": max(p["R_stretching"] for p in traj),
            "max_E": max(p["E"] for p in traj),
            "t_at_max_E": traj[int(np.argmax([p["E"] for p in traj]))]["t"],
            "max_tail_frac": max(p["tail_frac"] for p in traj),
            "energy_t0": traj[0]["energy"], "energy_tend": traj[-1]["energy"],
        }
        res["cases"]["taylor_green_time_integration"] = tgi

        # 6b. same vortex at four times the Reynolds number on a 64^3 grid,
        #     integrated only to t = 6 (before the dissipation peak), which is
        #     where the flow actually shows NET enstrophy growth, dE/dt > 0.
        N2, nu2 = 64, 0.01
        traj2 = integrate(N2, nu2, taylor_green(N2, 1.0), t_end=6.0, dt=0.02,
                          sample_every=15, pad=1.5)
        growth = [p for p in traj2 if p["dEdt"] > 0]
        tgi2 = {
            "N": N2, "nu": nu2, "dt": 0.02, "t_end": 6.0, "pad": 1.5,
            "trajectory": traj2,
            "max_E": max(p["E"] for p in traj2),
            "t_at_max_E": traj2[int(np.argmax([p["E"] for p in traj2]))]["t"],
            "E_t0": traj2[0]["E"],
            "n_samples_with_positive_dEdt": len(growth),
            "max_dEdt": max(p["dEdt"] for p in traj2),
            "t_at_max_dEdt": traj2[int(np.argmax([p["dEdt"] for p in traj2]))]["t"],
            "R_total_at_max_dEdt":
                traj2[int(np.argmax([p["dEdt"] for p in traj2]))]["R_total"],
            "max_R_total": max(p["R_total"] for p in traj2),
            "max_Q_shape": max(p["Q_shape"] for p in traj2),
            "max_tail_frac": max(p["tail_frac"] for p in traj2),
        }
        res["cases"]["taylor_green_time_integration_higher_Re"] = tgi2

    # ---- summary -------------------------------------------------------------
    all_static = tg + tubes + rnd
    qmax = max(d["Q_shape"] for d in all_static)
    qabs = max(abs(d["Q_shape"]) for d in all_static)
    res["summary"] = {
        "max_Q_shape_static_t0": qmax,
        "max_abs_Q_shape_static_t0": qabs,
        "max_Q_shape_evolved_taylor_green": tgi["max_Q_shape"] if tgi else None,
        "t_at_max_Q_evolved": tgi["t_at_max_Q_shape"] if tgi else None,
        "max_R_total_positive_growth_case":
            res["cases"]["taylor_green_time_integration_higher_Re"]["max_R_total"]
            if tgi else None,
        "max_R_stretching_static": max(d["R_stretching"] for d in all_static),
        "min_R_stretching_static": min(d["R_stretching"] for d in all_static),
        "max_identity_residual": max(s["identity_residual"] for s in sweep),
        "max_div_rel_l2": max(d["div_rel_l2"] for d in all_static),
        "max_enstrophy_tail_frac_static": max(
            d["enstrophy_tail_frac"] for d in all_static),
        "runtime_seconds": round(time.time() - t0, 2),
    }

    txt = json.dumps(res, indent=2, sort_keys=False)
    print(txt)
    if args.out:
        with open(args.out, "w") as f:
            f.write(txt + "\n")


if __name__ == "__main__":
    main()
