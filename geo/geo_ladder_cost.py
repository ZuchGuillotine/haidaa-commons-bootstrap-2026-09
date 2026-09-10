"""
geo_ladder_cost.py

Explicit cost arithmetic for a four-rung resolution ladder (Delta = 100/50/25/12.5 m)
of an idealized closed-domain tornado-like vortex in an open cloud model, plus the
storage arithmetic for the mean/eddy angular-momentum budget diagnostics.

Cost model, stated so a reader can rescale every number:

    N_points = (L_x/Delta) * (L_y/Delta) * (L_z/Delta)
    dt       = CFL * Delta / u_max                      (advective CFL; the split-explicit
                                                         acoustic substepping is folded into
                                                         the per-point-step cost below)
    N_steps  = T_sim / dt
    core_h   = N_points * N_steps * tau_pt / 3600

tau_pt is the sustained wall time per grid point per time step on one modern server
core, in core-seconds. It is the one empirical input and it is an ESTIMATE: it is
memory-bandwidth bound for a stencil code of this class. Everything scales linearly
in tau_pt, so a reader with a measured value can rescale in one multiplication.

Storage model:

    bytes_per_dump_3d = N_points * n_vars_3d * 4          (single precision)
    bytes_per_dump_2d = (L_r/Delta)*(L_z/Delta) * n_vars_2d * 4

The 2D figure is for ON-LINE azimuthal-mean and eddy-covariance accumulation in the
model, reduced to an (r,z) plane; the 3D figure is for OFF-LINE diagnosis from full
snapshots. The difference between them is the single largest lever on the project's
storage requirement.

No fitting, no random numbers, no seeds required.

Run:
    ./venv/bin/python commons/geo/artifacts/geo_ladder_cost.py
"""

import sys
import numpy as np

LX = LY = 8000.0     # m, horizontal extent of the idealized closed domain
LZ = 4000.0          # m, vertical extent
LR = 4000.0          # m, radial extent of the (r,z) diagnostic plane
DELTAS = (100.0, 50.0, 25.0, 12.5)
CFL = 0.5
U_MAX = 100.0        # m/s, peak wind that sets the advective CFL
T_SIM = 3600.0       # s of simulated time per member (spin-up plus contraction window)
TAU_PT = 4.0e-6      # core-seconds per grid point per time step  [ESTIMATE]

N_CLOSURES = 3       # Smagorinsky, TKE-1.5, fixed-viscosity control
N_SWIRL = 2          # two swirl-ratio settings, so the verdict does not rest on one regime
OVERHEAD = 1.20      # restarts, failed runs, spin-up reruns

N_VARS_3D = 10       # u, v, w, p, theta, qv, tke, plus 3 work/diagnostic fields
N_VARS_2D = 30       # mean fields plus all eddy covariances and budget terms in (r,z)
CADENCE = 1.0        # s, output cadence required to close the budget
T_WINDOW = 600.0     # s, the contraction window over which the budget must close
N_SNAPSHOTS_APRIORI = 50   # full 3D snapshots kept for the a-priori sub-grid test

GB = 1024.0 ** 3
TB = 1024.0 ** 4


def main():
    print("geo_ladder_cost.py")
    print("environment: python %s, numpy %s" % (sys.version.split()[0], np.__version__))
    print("command: ./venv/bin/python commons/geo/artifacts/geo_ladder_cost.py")
    print("no random numbers used; no seed required")
    print("")
    print("configuration: domain %.0f x %.0f x %.0f m, T_sim = %.0f s, CFL = %.2f,"
          % (LX, LY, LZ, T_SIM, CFL))
    print("               u_max = %.0f m/s, tau_pt = %.1e core-s per point-step [ESTIMATE]"
          % (U_MAX, TAU_PT))
    print("")

    hdr = ("%10s %14s %10s %10s %14s %14s"
           % ("Delta [m]", "N_points", "dt [s]", "N_steps", "point-steps", "core-hours"))
    print(hdr)
    print("-" * len(hdr))
    per_rung = []
    for d in DELTAS:
        npts = (LX / d) * (LY / d) * (LZ / d)
        dt = CFL * d / U_MAX
        nsteps = T_SIM / dt
        ps = npts * nsteps
        ch = ps * TAU_PT / 3600.0
        per_rung.append((d, npts, dt, nsteps, ps, ch))
        print("%10.1f %14.4g %10.4f %10.0f %14.4g %14.1f"
              % (d, npts, dt, nsteps, ps, ch))
    ladder_ch = sum(r[5] for r in per_rung)
    print("-" * len(hdr))
    print("%10s %14s %10s %10s %14s %14.1f" % ("TOTAL", "", "", "", "", ladder_ch))
    print("")
    print("the finest rung carries %.1f%% of the ladder cost; cost scales as Delta^-4"
          % (100.0 * per_rung[-1][5] / ladder_ch))
    print("(three spatial dimensions plus the CFL-shortened time step), so each halving")
    print("of Delta multiplies the cost by %.0f. Measured ratio 100 m -> 12.5 m: %.0f."
          % (16, per_rung[-1][5] / per_rung[0][5]))
    print("")

    print("=== full experiment ===")
    exp = ladder_ch * N_CLOSURES * N_SWIRL
    print("ladder x %d closures x %d swirl-ratio settings = %.0f core-hours"
          % (N_CLOSURES, N_SWIRL, exp))
    exp_oh = exp * OVERHEAD
    print("x %.2f for restarts, failed runs and spin-up reruns   = %.0f core-hours"
          % (OVERHEAD, exp_oh))
    analysis = 2000.0
    print("plus a-priori sub-grid filtering and budget analysis   = %.0f core-hours"
          % analysis)
    total = exp_oh + analysis
    print("TOTAL                                                  = %.0f core-hours" % total)
    print("")
    posted = 80000.0
    print("the project's posted pooled request is %.0f core-hours, leaving %.0f core-hours"
          % (posted, posted - total))
    print("of headroom (%.0f%%). The estimate is consistent with the posted request."
          % (100.0 * (posted - total) / posted))
    print("")
    print("what would break this estimate, in order of leverage:")
    print("  1. tau_pt: linear. A measured tau_pt of 1e-5 instead of %.0e would give"
          % TAU_PT)
    print("     %.0f core-hours and would exceed the posted request." % (total * 1e-5 / TAU_PT))
    print("  2. adding a fifth rung at Delta = 6.25 m: adds %.0f core-hours per"
          % (per_rung[-1][5] * 16.0))
    print("     closure-swirl member, %.0f in total. Not affordable inside the request."
          % (per_rung[-1][5] * 16.0 * N_CLOSURES * N_SWIRL * OVERHEAD))
    print("  3. a supercell-embedded member rather than a closed idealized domain:")
    print("     the domain volume grows by roughly 10^2-10^3, which is a separate request.")
    print("")

    print("=== storage ===")
    print("%10s %16s %16s %16s" % ("Delta [m]", "3D dump [GB]", "2D dump [MB]", "ratio"))
    n_dumps = T_WINDOW / CADENCE
    for d, npts, dt, nsteps, ps, ch in per_rung:
        b3 = npts * N_VARS_3D * 4.0
        n2 = (LR / d) * (LZ / d)
        b2 = n2 * N_VARS_2D * 4.0
        print("%10.1f %16.3f %16.2f %16.0f"
              % (d, b3 / GB, b2 / (1024.0 ** 2), b3 / b2))
    print("")
    d_fine, npts_fine = per_rung[-1][0], per_rung[-1][1]
    b3_fine = npts_fine * N_VARS_3D * 4.0
    n2_fine = (LR / d_fine) * (LZ / d_fine)
    b2_fine = n2_fine * N_VARS_2D * 4.0
    n_members = N_CLOSURES * N_SWIRL
    # coarse rungs add sum over (1/8, 1/64, 1/512) of the fine-rung volume in 3D
    coarse_factor_3d = sum((d_fine / d) ** 3 for d, *_ in per_rung) / 1.0
    coarse_factor_2d = sum((d_fine / d) ** 2 for d, *_ in per_rung) / 1.0

    offline = b3_fine * n_dumps * coarse_factor_3d * n_members
    online = b2_fine * n_dumps * coarse_factor_2d * n_members
    apriori = b3_fine * N_SNAPSHOTS_APRIORI * n_members

    print("plan A, OFF-LINE diagnosis: keep %.0f full 3D dumps at %.0f s cadence over the"
          % (n_dumps, CADENCE))
    print("        %.0f s contraction window, all rungs, all %d members:"
          % (T_WINDOW, n_members))
    print("        %.1f TB" % (offline / TB))
    print("plan B, ON-LINE azimuthal-mean and eddy accumulation reduced to an (r,z) plane,")
    print("        same cadence, plus %d full 3D snapshots per member kept for the"
          % N_SNAPSHOTS_APRIORI)
    print("        a-priori sub-grid test:")
    print("        %.3f TB (reduced) + %.1f TB (snapshots) = %.1f TB"
          % (online / TB, apriori / TB, (online + apriori) / TB))
    print("")
    print("posted pooled storage request: 60 TB.")
    print("plan A uses %.0f%% of it; plan B uses %.0f%%."
          % (100.0 * offline / (60 * TB), 100.0 * (online + apriori) / (60 * TB)))
    print("Implementing the azimuthal-mean/eddy accumulation INSIDE the model, rather than")
    print("diagnosing it from archived snapshots, reduces the storage requirement by a")
    print("factor of %.0f. That software task is therefore the highest-leverage single"
          % (offline / (online + apriori)))
    print("item in the project, and it is a laptop-scale job, not a compute-block job.")
    print("")
    print("=== cheap precursor: 2D axisymmetric mixing-length sweep ===")
    nr, nz = LR / 12.5, LZ / 12.5
    npts_ax = nr * nz
    dt_ax = CFL * 12.5 / U_MAX
    nsteps_ax = T_SIM / dt_ax
    ch_ax = npts_ax * nsteps_ax * TAU_PT / 3600.0
    print("an axisymmetric (r,z) run at Delta = 12.5 m has %.4g points and costs"
          % npts_ax)
    print("%.2f core-hours per member, i.e. %.0f times cheaper than the 3D fine rung."
          % (ch_ax, per_rung[-1][5] / ch_ax))
    print("A 24-member sweep over mixing length and swirl ratio costs %.0f core-hours"
          % (ch_ax * 24))
    print("and can be run on a workstation before any pooled allocation is requested.")


if __name__ == "__main__":
    main()
