ns_enstrophy_rate - reproducible enstrophy-growth-rate fixture
=============================================================

Purpose
-------
Measure the instantaneous enstrophy growth rate dE/dt of divergence-free
velocity fields on the periodic box [0,2pi]^3 and compare it with the
a-priori bound of the form dE/dt <= C * E^3 / nu^3 (Lu and Doering 2008).
The constant C is NOT hard-coded: it was not verified against the primary
source, so the fixture reports the raw dimensionless observables instead.

Files
-----
ns_enstrophy_rate.py           the fixture (single file, numpy + scipy only)
ns_enstrophy_rate_output.json  the actual output of the run described below
ns_enstrophy_rate_README.txt   this file

Definitions (the numbers depend on these; read before comparing)
----------------------------------------------------------------
omega = curl u
E     = (1/2) * volume-integral over [0,2pi]^3 of |omega|^2      (enstrophy)
P     = volume-integral of |curl omega|^2                        (palinstrophy)
S     = volume-integral of omega_i omega_j d_j u_i               (stretching)
dE/dt = -nu * P + S
These are volume INTEGRALS, not domain means. E has dimensions L^3 T^-2, so
R     = (dE/dt) * nu^3 / E^3      is dimensionless
R_S   = S * nu^3 / E^3            is dimensionless
With L := 2*pi the fixture also reports
Q     = S * L^(3/2) / E^(3/2)     dimensionless shape number, invariant under
                                  u -> lambda*u and independent of nu
Re_E  = sqrt(E * L) / nu          enstrophy-based Reynolds number
and these satisfy the exact identity  R_S = Q * Re_E^(-3).  That is algebra,
not an empirical scaling law; the fixture only checks that the code respects
it (max relative residual 1.6e-16 over 12 (amplitude, nu) combinations).

Convention traps: if you use Etilde = integral |omega|^2 (no 1/2) then
Rtilde = R/4.  If you use domain means instead of integrals, E is multiplied
by (2*pi)^-3 and R changes by (2*pi)^6 = 6.1e4.

Numerics
--------
Fourier pseudo-spectral with real FFTs (numpy.fft.rfftn / irfftn).
Static test fields are band-limited to |k_i| <= 15 on a 48^3 grid.  The
stretching integrand is CUBIC, so its grid quadrature is exact only when
3*kmax < N_eval; the fixture therefore evaluates S on a zero-padded grid
(96^3 for the 48^3 cases, 96^3 for the 64^3 time integration, pad=1.5).
Quadratic integrals (energy, E, P) are exact for 2*kmax < N.
Random fields use numpy.random.default_rng(seed) with seeds 11, 12, 13, 101.
Time integration: RK4, rotation-form nonlinearity, 2/3 dealiasing.

Correctness checks built into the run
-------------------------------------
1. Taylor-Green vortex u = A(sin x cos y cos z, -cos x sin y cos z, 0) at
   t = 0 has, by hand, E = 3 A^2 pi^3, P = 6 E (every Fourier mode of omega
   has |k|^2 = 3), and S = 0 (every term in the stretching integrand carries
   an odd power of cos x).  Observed relative errors in E and P: 6.1e-16 at
   A = 0.5, 1, 2.  Observed |S| <= 4.6e-17 against an exact 0.
2. Divergence: max relative L2 divergence over all static fields 4.2e-15.
3. Mirror-symmetric antiparallel tube pair: omega_1, omega_2 even and
   omega_3 odd about z = pi forces S = 0 identically.  Observed
   S = 1.7e-19, Q = 7.3e-19.
4. Grid independence: the symmetry-broken tube pair gives identical Q on the
   48^3 and 64^3 grids to a relative difference of 2.0e-15.
5. Amplitude and viscosity invariance of Q, and the identity
   R_S = Q * Re_E^(-3): residual <= 1.6e-16.

Environment used for ns_enstrophy_rate_output.json
--------------------------------------------------
python 3.14.4
numpy  2.5.3
scipy  1.18.1   (imported only to record its version)
platform macOS-26.4 arm64 (Apple M1), Darwin 25.4.0
No BLAS-sensitive operations; results are FFT-only and should reproduce to
near machine precision on any IEEE-754 double platform, and to at least
10 significant figures across FFT library versions.

Exact commands
--------------
Interpreter:
  VENV=/private/tmp/claude-501/-Users-benjamincox-dsm-haidaa/\
da69bd54-2977-4882-88ec-9babfab4e89d/scratchpad/venv/bin/python

Full run (writes the JSON artifact):
  cd <artifacts dir>
  $VENV ns_enstrophy_rate.py --out ns_enstrophy_rate_output.json

Fast subset (skips both time integrations, ~8 s):
  $VENV ns_enstrophy_rate.py --quick

Runtime
-------
Full run: 105 s wall clock (self-reported summary.runtime_seconds = 105.31)
on an Apple M1 laptop, single process.  --quick: about 8 s.

Headline numbers from ns_enstrophy_rate_output.json
---------------------------------------------------
Kinematic (t = 0) fields, shape number Q:
  Taylor-Green, A = 0.5, 1, 2                 Q = 0 exactly (symmetry)
  antiparallel tubes, mirror-symmetric        Q = 7.3e-19 (0 by symmetry)
  antiparallel tubes, symmetry-broken         Q = -3.8365e-4
  random fields, E(k) ~ k^-p, p = 1, 5/3, 3,
    seeds 11/12/13, enstrophy normalized to 1 Q in [-1.82e-2, +3.69e-2]
Dynamically evolved Taylor-Green:
  N = 48, nu = 0.025, t in [0,10]: Q rises to 0.5204 at t = 3.4; E decays
    monotonically from 93.019; no net enstrophy growth at this viscosity.
  N = 64, nu = 0.01, t in [0,6]: E grows from 93.019 to 160.88 at t = 4.8
    (+73%); dE/dt > 0 at 15 of 21 samples; max dE/dt = 25.026 at t = 2.4;
    max Q = 0.6216; Re_E in [2418, 3179]; max R = 1.682e-11 at t = 1.5.

Interpretation in one line: with Q_max ~ 0.62, the identity R = Q Re_E^-3
says that a field of this shape could only reach R = 1 at Re_E ~ 0.85, i.e.
in a regime with essentially no inertia.  Nothing here bears on global
regularity.

Known limitations
-----------------
- The "enstrophy_tail_frac" field is a resolution indicator only for the
  time integrations.  For the exactly band-limited static fields it is large
  for shallow spectra (up to 0.556 for E(k) ~ k^-1) purely because the
  spectrum is flat out to the cutoff; the quadrature there is still exact.
- The 64^3, nu = 0.01 run is well resolved up to about t = 4
  (tail fraction < 1.5e-3) and only marginally resolved at t = 6
  (3.5e-3).  Treat t > 5 as indicative.
- Only three field families are tested.  No optimization over initial data
  was performed, so nothing here is an upper bound on attainable Q.
