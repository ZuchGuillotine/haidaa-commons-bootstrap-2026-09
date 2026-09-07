#!/usr/bin/env python3
"""dist_checker.py -- a tiny explicit-state breadth-first model checker for
toy distributed-protocol models, plus three models:

  1. lease_mutex             -- lease-based mutual exclusion with clock drift
  2. two_phase_commit        -- 2PC with a crash-stop coordinator and an
                                optional participant timeout-abort
  3. paxos_nondurable_promise-- single-decree Paxos, 3 acceptors, promise
                                optionally lost on crash-recovery

Design
------
A model is an object with:
    name          : str
    params        : dict of explicit parameters (fault model knobs)
    initial()     : the initial state (a hashable tuple)
    actions(s)    : iterable of (label:str, s2:state) successors
    invariant(s)  : True if the safety property holds in s
    render(s)     : dict, a human-readable snapshot of s

bfs_check() does a level-by-level breadth-first search from the initial
state. Because BFS visits states in non-decreasing distance from the
initial state, the FIRST state found that violates the invariant is at
minimal distance; the reconstructed schedule is therefore a violating
schedule of minimal LENGTH IN ACTIONS OF THIS MODEL'S ALPHABET. Minimality
is relative to the action granularity chosen by the model, not an absolute
property of the protocol.

Everything here is a toy, inert model. Nothing targets any live system.
Absence of a counterexample below the stated depth/state bounds is NOT a
proof of correctness; it is a bounded-exploration result.

Determinism: no randomness is used anywhere. Successors are generated in a
fixed order, so the reported trace is reproducible byte-for-byte.

Usage:
    python dist_checker.py --outdir .
"""

import argparse
import collections
import json
import os
import platform
import sys
import time

# ----------------------------------------------------------------------------
# The checker
# ----------------------------------------------------------------------------


def bfs_check(model, max_depth=64, max_states=4_000_000, stop_on_violation=True,
              find_deadlocks=False):
    """Breadth-first explicit-state search.

    Returns a dict with state counts, depth reached, runtime, and (if found)
    the minimal-length violating schedule with per-step state snapshots.
    """
    t0 = time.time()
    s0 = model.initial()
    parent = {s0: (None, None)}          # state -> (predecessor, action label)
    frontier = [s0]
    depth = 0
    explored = 1
    transitions = 0
    violation = None
    deadlocks = []
    truncated_depth = False
    truncated_states = False

    def trace_of(s):
        chain = []
        cur = s
        while True:
            pred, label = parent[cur]
            if pred is None:
                break
            chain.append((label, cur))
            cur = pred
        chain.reverse()
        steps = [{"step": 0, "action": "<initial>", "state": model.render(s0)}]
        for i, (label, st) in enumerate(chain, start=1):
            steps.append({"step": i, "action": label, "state": model.render(st)})
        return steps

    if not model.invariant(s0):
        violation = {"depth": 0, "trace": trace_of(s0)}
        frontier = []

    while frontier and violation is None:
        if depth >= max_depth:
            truncated_depth = True
            break
        nxt = []
        for s in frontier:
            succs = list(model.actions(s))
            transitions += len(succs)
            if find_deadlocks and not succs and model.is_blocked(s):
                deadlocks.append({"depth": depth, "trace": trace_of(s)})
            for label, s2 in succs:
                if s2 in parent:
                    continue
                parent[s2] = (s, label)
                explored += 1
                if not model.invariant(s2):
                    violation = {"depth": depth + 1, "trace": trace_of(s2)}
                    if stop_on_violation:
                        nxt = []
                        break
                nxt.append(s2)
                if explored >= max_states:
                    truncated_states = True
                    break
            if violation is not None and stop_on_violation:
                break
            if truncated_states:
                break
        if truncated_states:
            break
        frontier = nxt
        depth += 1

    elapsed = time.time() - t0
    exhaustive = (violation is None and not truncated_depth
                  and not truncated_states and not frontier)
    return {
        "model": model.name,
        "params": dict(model.params),
        "states_explored": explored,
        "transitions_examined": transitions,
        "max_depth_reached": depth,
        "max_depth_bound": max_depth,
        "max_states_bound": max_states,
        "runtime_seconds": round(elapsed, 3),
        "violation_found": violation is not None,
        "violation_depth": violation["depth"] if violation else None,
        "violation_trace": violation["trace"] if violation else None,
        "deadlock_count": len(deadlocks),
        "deadlocks": deadlocks[:3],
        "truncated_by_depth_bound": truncated_depth,
        "truncated_by_state_bound": truncated_states,
        "state_space_exhausted": exhaustive,
    }


# ----------------------------------------------------------------------------
# Model 1: lease_mutex
# ----------------------------------------------------------------------------

class LeaseMutex(object):
    """Lease-based mutual exclusion with a bounded-lag client clock.

    Two clients, one lock server. Discrete time. server_clock = t.
    Client i's local clock is (t - lag[i]); lag[i] only grows (the client's
    clock ticks no faster than the server's) and is bounded by MAX_LAG.
    This is a bounded-SKEW model of drift: over the whole run the client
    clock may fall behind the server clock by at most MAX_LAG ticks.

    Protocol (the standard grant / expire lease):
      - The server grants a lease of length L to client i when the previous
        lease has expired by the SERVER's clock (t >= lease_end).
      - Client i records its own local clock at grant time, lg[i], and
        believes it holds the lock while
              local_clock - lg[i] < L - GUARD
        where GUARD is the client-side safety margin the protocol uses to
        absorb clock error.

    Safety: at most one client believes it holds the lock at any server time.

    The protocol is safe exactly when GUARD >= MAX_LAG. Grant delivery is
    instantaneous (no message delay is modeled) so the only source of
    disagreement is clock error.
    """

    NONE = -1

    def __init__(self, L=3, guard=1, max_lag=2, horizon=8, name=None):
        self.L = L
        self.guard = guard
        self.max_lag = max_lag
        self.horizon = horizon
        self.name = name or "lease_mutex"
        self.params = {
            "clients": 2,
            "lease_length_L_server_ticks": L,
            "client_guard_margin_ticks": guard,
            "max_client_clock_lag_ticks": max_lag,
            "server_time_horizon": horizon,
            "grant_delivery": "instantaneous",
            "clock_model": "client clock ticks 0 or 1 per server tick; lag "
                           "monotone non-decreasing, capped at max_lag",
            "protocol_is_safe_iff": "guard >= max_lag",
        }

    def initial(self):
        # (t, lag0, lag1, lease_end, lg0, lg1)
        return (0, 0, 0, 0, self.NONE, self.NONE)

    def believes(self, s, i):
        t, lag0, lag1, _le, lg0, lg1 = s
        lag = (lag0, lag1)[i]
        lg = (lg0, lg1)[i]
        if lg == self.NONE:
            return False
        return (t - lag) - lg < self.L - self.guard

    def invariant(self, s):
        return not (self.believes(s, 0) and self.believes(s, 1))

    def actions(self, s):
        t, lag0, lag1, le, lg0, lg1 = s
        out = []
        # server grants the lease to client i once the old lease expired
        if t >= le:
            for i in (0, 1):
                lag = (lag0, lag1)[i]
                nlg = [lg0, lg1]
                nlg[i] = t - lag
                out.append(("grant(server->c%d)@t=%d" % (i, t),
                            (t, lag0, lag1, t + self.L, nlg[0], nlg[1])))
        # a tick of the server clock; each client's clock independently
        # either keeps up (d=0) or loses one tick (d=1, bounded by max_lag)
        if t < self.horizon:
            for d0 in (0, 1):
                for d1 in (0, 1):
                    n0, n1 = lag0 + d0, lag1 + d1
                    if n0 > self.max_lag or n1 > self.max_lag:
                        continue
                    lbl = "tick(server->%d, c0_%s, c1_%s)" % (
                        t + 1,
                        "lags" if d0 else "keeps_up",
                        "lags" if d1 else "keeps_up")
                    out.append((lbl, (t + 1, n0, n1, le, lg0, lg1)))
        return out

    def is_blocked(self, s):
        return False

    def render(self, s):
        t, lag0, lag1, le, lg0, lg1 = s
        return {
            "server_clock": t,
            "lease_end_server_clock": le,
            "client0": {"local_clock": t - lag0, "lag": lag0,
                        "lease_start_local": None if lg0 == self.NONE else lg0,
                        "believes_it_holds_lock": self.believes(s, 0)},
            "client1": {"local_clock": t - lag1, "lag": lag1,
                        "lease_start_local": None if lg1 == self.NONE else lg1,
                        "believes_it_holds_lock": self.believes(s, 1)},
        }


# ----------------------------------------------------------------------------
# Model 2: two_phase_commit
# ----------------------------------------------------------------------------

INIT, YES, COMMIT, ABORT = 0, 1, 2, 3
PSTATE_NAME = {INIT: "init", YES: "voted_yes_uncertain",
               COMMIT: "decided_commit", ABORT: "decided_abort"}
C_INIT, C_COMMIT, C_ABORT = 0, 1, 2
CSTATE_NAME = {C_INIT: "collecting", C_COMMIT: "decided_commit",
               C_ABORT: "decided_abort"}
VNONE, VYES, VNO = 0, 1, 2
VOTE_NAME = {VNONE: None, VYES: "yes", VNO: "no"}
DNONE, DCOMMIT, DABORT = 0, 1, 2
DEC_NAME = {DNONE: None, DCOMMIT: "commit", DABORT: "abort"}


class TwoPhaseCommit(object):
    """2PC: one coordinator, N participants, crash-stop coordinator.

    Fault parameters:
      timeout_abort : if True, a participant that has voted YES and is
                      uncertain may unilaterally decide ABORT on a timeout.
                      This is the "timeout-abort optimization".
      coordinator_may_crash : crash-stop; once crashed the coordinator takes
                      no further steps and never recovers (no durable log).

    Channels are modeled as one in-flight slot per (link, direction); no
    loss, no duplication; ordering between different links is free.

    Safety: no two participants decide differently.
    Liveness (checked separately): a reachable state with no enabled action
    in which some participant is still uncertain is a BLOCKING state.
    """

    def __init__(self, n=2, timeout_abort=True, coordinator_may_crash=True,
                 name=None):
        self.n = n
        self.timeout_abort = timeout_abort
        self.coord_may_crash = coordinator_may_crash
        self.name = name or "two_phase_commit"
        self.params = {
            "participants": n,
            "timeout_abort_optimization": timeout_abort,
            "coordinator_may_crash": coordinator_may_crash,
            "coordinator_crash_model": "crash-stop, no recovery, no durable log",
            "channels": "one in-flight slot per link, reliable, no duplication",
            "votes_allowed": ["yes", "no"],
        }

    def initial(self):
        n = self.n
        # (pstate..., vote_inflight..., vote_recv..., dec_inflight...,
        #  dec_sent..., coord, crashed)
        z = tuple([0] * n)
        return (z, z, z, z, z, C_INIT, 0)

    def invariant(self, s):
        p = s[0]
        return not (COMMIT in p and ABORT in p)

    def is_blocked(self, s):
        p = s[0]
        return any(x in (INIT, YES) for x in p)

    def actions(self, s):
        p, vi, vr, di, ds, coord, crashed = s
        n = self.n
        out = []

        def rep(tup, i, v):
            l = list(tup)
            l[i] = v
            return tuple(l)

        # participants vote
        for i in range(n):
            if p[i] == INIT:
                out.append(("p%d_votes_yes" % i,
                            (rep(p, i, YES), rep(vi, i, VYES), vr, di, ds,
                             coord, crashed)))
                # a NO voter aborts unilaterally (always safe)
                out.append(("p%d_votes_no_and_aborts" % i,
                            (rep(p, i, ABORT), rep(vi, i, VNO), vr, di, ds,
                             coord, crashed)))
        # coordinator receives a vote
        if not crashed:
            for i in range(n):
                if vi[i] != VNONE:
                    out.append(("coord_receives_vote_from_p%d(%s)"
                                % (i, VOTE_NAME[vi[i]]),
                                (p, rep(vi, i, VNONE), rep(vr, i, vi[i]), di,
                                 ds, coord, crashed)))
            # coordinator decides
            if coord == C_INIT:
                if all(v == VYES for v in vr):
                    out.append(("coord_decides_COMMIT",
                                (p, vi, vr, di, ds, C_COMMIT, crashed)))
                elif any(v == VNO for v in vr):
                    out.append(("coord_decides_ABORT",
                                (p, vi, vr, di, ds, C_ABORT, crashed)))
            else:
                dv = DCOMMIT if coord == C_COMMIT else DABORT
                for i in range(n):
                    if ds[i] == 0:
                        out.append(("coord_sends_%s_to_p%d"
                                    % (DEC_NAME[dv], i),
                                    (p, vi, vr, rep(di, i, dv), rep(ds, i, 1),
                                     coord, crashed)))
            if self.coord_may_crash:
                out.append(("coord_CRASHES",
                            (p, vi, vr, di, ds, coord, 1)))
        # decision delivery
        for i in range(n):
            if di[i] != DNONE:
                if p[i] in (INIT, YES):
                    np_ = COMMIT if di[i] == DCOMMIT else ABORT
                    out.append(("p%d_receives_%s" % (i, DEC_NAME[di[i]]),
                                (rep(p, i, np_), vi, vr, rep(di, i, DNONE),
                                 ds, coord, crashed)))
                else:
                    out.append(("p%d_ignores_late_%s (already decided %s)"
                                % (i, DEC_NAME[di[i]], PSTATE_NAME[p[i]]),
                                (p, vi, vr, rep(di, i, DNONE), ds, coord,
                                 crashed)))
        # timeout-abort optimization
        if self.timeout_abort:
            for i in range(n):
                if p[i] == YES:
                    out.append(("p%d_TIMES_OUT_and_unilaterally_ABORTS" % i,
                                (rep(p, i, ABORT), vi, vr, di, ds, coord,
                                 crashed)))
        return out

    def render(self, s):
        p, vi, vr, di, ds, coord, crashed = s
        return {
            "coordinator": {"state": CSTATE_NAME[coord],
                            "crashed": bool(crashed),
                            "votes_received": [VOTE_NAME[v] for v in vr]},
            "participants": [
                {"id": i, "state": PSTATE_NAME[p[i]],
                 "vote_in_flight": VOTE_NAME[vi[i]],
                 "decision_in_flight": DEC_NAME[di[i]]}
                for i in range(self.n)],
        }


# ----------------------------------------------------------------------------
# Model 3: paxos_nondurable_promise
# ----------------------------------------------------------------------------

VAL_NAME = {0: None, 1: "X", 2: "Y"}


class PaxosNondurablePromise(object):
    """Single-decree Paxos, 3 acceptors, 2 proposers, ballots 1 and 2.

    Proposer 1 uses ballot 1 and prefers value X; proposer 2 uses ballot 2
    and prefers value Y (fixing the preferences is a symmetry reduction).

    Fault parameter:
      durable_promise : if False, an acceptor that crashes and recovers
                        loses its highest-promised ballot (reset to 0) while
                        its accepted (ballot, value) pair survives. This is
                        the "promise written to volatile memory only" bug.

    Simplification: prepare and accept are atomic request/response steps
    (send + acceptor handling + response delivery in one action). Message
    reordering and loss are therefore NOT explored; this keeps the state
    space small and isolates durability as the cause. Rejections are
    stutter steps and are not generated.

    "Chosen" is tracked as a monotone ghost variable: a value v enters
    `chosen` as soon as a majority of acceptors simultaneously hold the same
    accepted (ballot, v), and it never leaves.

    Safety: at most one value is ever chosen.
    """

    QUORUM = 2  # majority of 3

    def __init__(self, durable_promise=False, prefs=(1, 2), name=None):
        self.durable = durable_promise
        self.prefs = tuple(prefs)   # prefs[i] = value preferred at ballot i+1
        self.nprop = len(self.prefs)
        self.name = name or "paxos_nondurable_promise"
        self.params = {
            "acceptors": 3,
            "proposers": self.nprop,
            "ballots": list(range(1, self.nprop + 1)),
            "proposer_preferred_values": {
                "ballot%d" % (i + 1): VAL_NAME[v]
                for i, v in enumerate(self.prefs)},
            "quorum_size": self.QUORUM,
            "durable_promise": durable_promise,
            "crash_recovery_effect": ("resets highest-promised ballot to 0; "
                                      "accepted (ballot,value) survives"
                                      if not durable_promise else "none"),
            "message_model": "atomic RPC; no reordering, loss or duplication",
        }

    def initial(self):
        # acceptors: (promised, acc_ballot, acc_value) x 3
        # proposers: (phase, promise_mask, best_ballot_seen, best_value_seen)
        #            phase 1 = preparing, 2 = accepting (value fixed)
        acc = ((0, 0, 0), (0, 0, 0), (0, 0, 0))
        props = tuple((1, 0, 0, 0) for _ in range(self.nprop))
        return (acc, props, ())  # third slot = sorted tuple of chosen values

    def invariant(self, s):
        return len(s[2]) <= 1

    def is_blocked(self, s):
        return False

    def _chosen_after(self, acc, chosen):
        counts = {}
        for (pr, ab, av) in acc:
            if ab > 0:
                counts[(ab, av)] = counts.get((ab, av), 0) + 1
        new = set(chosen)
        for (ab, av), c in counts.items():
            if c >= self.QUORUM:
                new.add(av)
        return tuple(sorted(new))

    def actions(self, s):
        acc, props, chosen = s
        out = []
        for pi in range(self.nprop):
            ballot = pi + 1
            pref = self.prefs[pi]
            phase, mask, bb, bv = props[pi]
            if phase == 1:
                for a in range(3):
                    if mask & (1 << a):
                        continue
                    pr, ab, av = acc[a]
                    if pr >= ballot:
                        continue  # rejected: stutter step, not generated
                    nacc = list(acc)
                    nacc[a] = (ballot, ab, av)
                    nbb, nbv = (ab, av) if ab > bb else (bb, bv)
                    nprops = list(props)
                    nprops[pi] = (1, mask | (1 << a), nbb, nbv)
                    out.append(("b%d_prepare_acceptor%s (promise recorded)"
                                % (ballot, "ABC"[a]),
                                (tuple(nacc), tuple(nprops), chosen)))
                if bin(mask).count("1") >= self.QUORUM:
                    val = bv if bb > 0 else pref
                    nprops = list(props)
                    nprops[pi] = (2, mask, 0, val)
                    out.append(("b%d_picks_value_%s_from_promise_quorum"
                                % (ballot, VAL_NAME[val]),
                                (acc, tuple(nprops), chosen)))
            else:
                val = bv
                for a in range(3):
                    pr, ab, av = acc[a]
                    if pr > ballot:
                        continue  # rejected
                    if (pr, ab, av) == (ballot, ballot, val):
                        continue  # no state change
                    nacc = list(acc)
                    nacc[a] = (ballot, ballot, val)
                    nacc = tuple(nacc)
                    out.append(("b%d_accept_%s_at_acceptor%s"
                                % (ballot, VAL_NAME[val], "ABC"[a]),
                                (nacc, props, self._chosen_after(nacc, chosen))))
        if not self.durable:
            for a in range(3):
                pr, ab, av = acc[a]
                if pr == 0:
                    continue  # nothing to lose
                nacc = list(acc)
                nacc[a] = (0, ab, av)
                out.append(("acceptor%s_CRASHES_and_RECOVERS "
                            "(promise b%d lost, accepted state kept)"
                            % ("ABC"[a], pr),
                            (tuple(nacc), props, chosen)))
        return out

    def render(self, s):
        acc, props, chosen = s
        return {
            "acceptors": [
                {"id": "ABC"[a], "highest_promised_ballot": acc[a][0],
                 "accepted_ballot": acc[a][1],
                 "accepted_value": VAL_NAME[acc[a][2]]}
                for a in range(3)],
            "proposers": [
                {"ballot": p + 1,
                 "phase": "prepare" if props[p][0] == 1 else "accept",
                 "promises_from": [c for c, a in zip("ABC", range(3))
                                   if props[p][1] & (1 << a)],
                 "value": (VAL_NAME[props[p][3]] if props[p][0] == 2
                           else None)}
                for p in range(self.nprop)],
            "chosen_values": [VAL_NAME[v] for v in chosen],
        }


# ----------------------------------------------------------------------------
# Runner
# ----------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()
    out = os.path.abspath(args.outdir)
    os.makedirs(out, exist_ok=True)

    runs = []

    def run(model, tag, note, **kw):
        sys.stderr.write("running %-46s ... " % tag)
        sys.stderr.flush()
        r = bfs_check(model, **kw)
        r["run_tag"] = tag
        r["note"] = note
        sys.stderr.write("%8d states  depth %2d  %6.2fs  violation=%s\n"
                         % (r["states_explored"], r["max_depth_reached"],
                            r["runtime_seconds"], r["violation_found"]))
        runs.append(r)
        return r

    # --- lease_mutex ------------------------------------------------------
    m = LeaseMutex(L=3, guard=1, max_lag=2, horizon=8)
    r_lease_bad = run(m, "lease_mutex/drift_exceeds_guard",
                      "guard=1 < max_lag=2: protocol margin too small",
                      max_depth=32)
    run(LeaseMutex(L=3, guard=2, max_lag=2, horizon=12),
        "lease_mutex/guard_equals_max_lag",
        "guard=2 >= max_lag=2: expected safe within bound", max_depth=64)
    run(LeaseMutex(L=3, guard=1, max_lag=1, horizon=12),
        "lease_mutex/guard_equals_max_lag_small",
        "guard=1 >= max_lag=1: expected safe within bound", max_depth=64)
    run(LeaseMutex(L=4, guard=3, max_lag=4, horizon=12),
        "lease_mutex/large_drift_small_guard",
        "guard=3 < max_lag=4: expected unsafe", max_depth=64)
    run(LeaseMutex(L=4, guard=3, max_lag=3, horizon=20),
        "lease_mutex/guard_equals_max_lag_wide",
        "guard=3 >= max_lag=3, longer horizon: expected safe within bound",
        max_depth=128)

    # --- two_phase_commit -------------------------------------------------
    r_2pc_bad = run(TwoPhaseCommit(n=2, timeout_abort=True),
                    "two_phase_commit/n2_timeout_abort_on",
                    "timeout-abort optimization enabled: expected unsafe",
                    max_depth=32)
    run(TwoPhaseCommit(n=3, timeout_abort=True),
        "two_phase_commit/n3_timeout_abort_on",
        "3 participants, timeout-abort enabled", max_depth=32)
    r_2pc_ok = run(TwoPhaseCommit(n=2, timeout_abort=False),
                   "two_phase_commit/n2_timeout_abort_off",
                   "no optimization: expected safe, but blocking states exist",
                   max_depth=64, find_deadlocks=True)
    run(TwoPhaseCommit(n=3, timeout_abort=False),
        "two_phase_commit/n3_timeout_abort_off",
        "3 participants, no optimization", max_depth=64, find_deadlocks=True)
    run(TwoPhaseCommit(n=2, timeout_abort=False, coordinator_may_crash=False),
        "two_phase_commit/n2_no_opt_no_crash",
        "no optimization, no coordinator crash: expected safe and non-blocking",
        max_depth=64, find_deadlocks=True)

    # --- paxos ------------------------------------------------------------
    r_paxos_bad = run(PaxosNondurablePromise(durable_promise=False),
                      "paxos/nondurable_promise",
                      "promise lost on crash-recovery: expected unsafe",
                      max_depth=24)
    run(PaxosNondurablePromise(durable_promise=True),
        "paxos/durable_promise",
        "promise persisted: expected safe within bound", max_depth=24)
    run(PaxosNondurablePromise(durable_promise=True, prefs=(1, 2, 1)),
        "paxos/durable_promise_3_ballots",
        "3 ballots (values X,Y,X), promise persisted: expected safe",
        max_depth=40)
    run(PaxosNondurablePromise(durable_promise=False, prefs=(1, 2, 1)),
        "paxos/nondurable_promise_3_ballots",
        "3 ballots, promise lost on recovery: expected unsafe", max_depth=40)

    # --- write counterexample artifacts -----------------------------------
    ce_files = {
        "lease_mutex": r_lease_bad,
        "two_phase_commit": r_2pc_bad,
        "paxos_nondurable_promise": r_paxos_bad,
    }
    for tag, r in ce_files.items():
        path = os.path.join(out, "dist_counterexample_%s.json" % tag)
        with open(path, "w") as f:
            json.dump({
                "model": r["model"],
                "run_tag": r["run_tag"],
                "params": r["params"],
                "invariant_violated": INVARIANT_TEXT[tag],
                "minimal_violating_schedule_length": r["violation_depth"],
                "minimality": ("minimal in the number of actions of this "
                               "model's action alphabet; found by BFS"),
                "states_explored_before_violation": r["states_explored"],
                "runtime_seconds": r["runtime_seconds"],
                "steps": r["violation_trace"],
            }, f, indent=2)
        sys.stderr.write("wrote %s\n" % path)

    # 2PC blocking schedule gets its own artifact section
    blocking = r_2pc_ok["deadlocks"]
    summary = {
        "generated_by": "dist_checker.py",
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "randomness": "none; fully deterministic",
        "invariants": INVARIANT_TEXT,
        "runs": runs,
        "two_phase_commit_blocking_example": {
            "note": ("terminal state with no enabled action in which a "
                     "participant is still uncertain (blocking = liveness "
                     "failure, not a safety violation)"),
            "run_tag": r_2pc_ok["run_tag"],
            "terminal_blocking_states_found": r_2pc_ok["deadlock_count"],
            "shortest_examples": blocking,
        },
        "caveats": [
            "Bounded exploration. 'No violation' means no violation within "
            "the stated depth, state and parameter bounds; it is not a proof.",
            "Counterexamples are counterexamples for these toy models under "
            "these explicit fault parameters, not statements about any "
            "real implementation.",
            "Minimal length is relative to each model's action granularity.",
        ],
    }
    path = os.path.join(out, "dist_checker_output.json")
    with open(path, "w") as f:
        json.dump(summary, f, indent=2)
    sys.stderr.write("wrote %s\n" % path)


INVARIANT_TEXT = {
    "lease_mutex": ("MutualExclusion: at no server time do both clients "
                    "believe they hold the lease"),
    "two_phase_commit": ("AgreementOnDecision: no two participants decide "
                         "differently"),
    "paxos_nondurable_promise": ("AtMostOneChosen: at most one value is ever "
                                 "chosen by a quorum"),
}


if __name__ == "__main__":
    main()
