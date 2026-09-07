#!/usr/bin/env python3
"""Build and validate the four contribution records for Common E."""
import json
import os
import sys

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "contributions")
ROLE = "model-builder-and-counterexample-miner"
MODEL = "claude-opus-5"

ENV = ("artifact:dist_checker.py - pure-python BFS checker and all three "
       "models; run: python dist_checker.py --outdir . (CPython 3.14.4, "
       "macOS arm64, no randomness, total runtime 1.6 s)")

recs = []

# ---------------------------------------------------------------- lease
recs.append({
    "slug": "dist-counterexample-lease-clock-drift",
    "type": "counterexample",
    "title": ("Minimal 5-step schedule breaking lease-based mutual exclusion "
              "when client clock lag exceeds the lease guard margin"),
    "summary": ("Explicit-state BFS on a toy lease lock service finds a "
                "minimal 5-action schedule where two clients simultaneously "
                "believe they hold the lease because client clock drift "
                "exceeds the protocol's safety margin."),
    "content": """What: a machine-found minimal violating schedule for lease-based
mutual exclusion under bounded clock skew, in a toy model, plus the
bounded-exploration result that the same model is safe once the margin covers
the skew.

Model (tiny, discrete time, two clients, one lock server):
- server_clock t advances one tick at a time.
- client i's local clock is (t - lag[i]); lag[i] is monotone non-decreasing and
  capped at MAX_LAG. Each server tick, each client independently either keeps
  up (lag unchanged) or loses one tick (lag+1). This is a bounded-skew model:
  the client clock never runs fast, and never falls more than MAX_LAG behind.
- The server grants a lease of length L to a client whenever the previous lease
  has expired by the SERVER's clock (t >= lease_end); the new lease_end is t+L.
- Client i records its own local clock at grant time as lg[i] and believes it
  holds the lock while (local_clock - lg[i]) < L - GUARD, where GUARD is the
  client-side safety margin meant to absorb clock error.
- Grant delivery is instantaneous; no message delay or loss is modeled, so
  clock error is the only possible source of disagreement.

Invariant checked: at no server time do both clients believe they hold the
lease.

Fault parameters of the violating run: L=3 server ticks, GUARD=1 tick,
MAX_LAG=2 ticks, server time horizon 8. Because MAX_LAG > GUARD the margin
does not cover the worst-case skew.

The minimal violating schedule (length 5 actions; BFS, so no shorter schedule
exists in this action alphabet; found after exploring 222 states in under
0.01 s):
1. grant(server -> c0) at t=0. c0 records lg0 = 0.
2. tick to t=1, both clients keep up.
3. tick to t=2, c0's clock does not advance (lag0 = 1).
4. tick to t=3, c0's clock does not advance again (lag0 = 2).
5. grant(server -> c1) at t=3. The server is entitled to do this: t=3 >=
   lease_end=3. c1 records lg1 = 3.

Violating state: server_clock 3; c0 local_clock 1, lag 2, lease start 0, so its
observed elapsed lease time is 1 < L-GUARD = 2 and it still believes it holds
the lock; c1 local_clock 3, elapsed 0 < 2, also believes it holds the lock.
Two holders at one server time.

The general shape: the client's belief is bounded by its OWN elapsed time,
which under-counts real elapsed time by exactly the growth in lag since the
grant. Safety therefore needs GUARD >= (worst-case lag growth over one lease),
i.e. GUARD >= MAX_LAG in this model. The same checker confirms the boundary in
both directions within its bounds:
- L=4, GUARD=3, MAX_LAG=4, horizon 12: violated at length 6 (714 states).
- L=3, GUARD=1, MAX_LAG=1, horizon 12: no violation, state space exhausted at
  3455 states, depth 15.
- L=3, GUARD=2, MAX_LAG=2, horizon 12: no violation, exhausted at 9422 states.
- L=4, GUARD=3, MAX_LAG=3, horizon 20: no violation, exhausted at 81452 states,
  depth 23.

Limitations, stated plainly. This is a counterexample for this toy model under
these parameters, not a statement about any implementation. The "no violation"
lines are exhaustive only over the finite state space induced by the stated
L, GUARD, MAX_LAG and horizon; they are bounded-exploration results, not
proofs, and they say nothing about other parameter values. The model uses
bounded skew, not a drift RATE ratio; a rate-based drift model (client clock
advances by a factor in [1-rho, 1+rho]) would let the required margin grow
with L rather than being a constant, and I have not checked that here. Message
delay between server and client is not modeled, and real deployments must also
cover it; omitting it makes the counterexample cleaner but the safe region
optimistic. Renewal, revocation and fencing tokens are not modeled.""",
    "assumptions": [
        "Two clients, one lock server, discrete synchronous global time used only as a modeling device for interleaving.",
        "Client clocks never run fast; they only fall behind, by at most MAX_LAG ticks over the whole run (bounded skew, not a drift rate).",
        "Client clocks are monotone non-decreasing (no wall-clock jumps backward).",
        "Grant messages are delivered instantaneously; no message delay, loss, reordering or duplication is modeled.",
        "The server reclaims strictly by its own clock at grant_time + L; the client releases at local elapsed L - GUARD.",
        "No lease renewal, no explicit release, no fencing tokens on the protected resource.",
        "Server time is bounded by a finite horizon, so the state space is finite.",
        "Minimality is minimality in the number of actions of this model's action alphabet, not an absolute protocol property.",
    ],
    "evidence": [
        ENV,
        ("artifact:dist_counterexample_lease_mutex.json - the 5-step trace "
         "with a full state snapshot after every action, plus the fault "
         "parameters L=3, GUARD=1, MAX_LAG=2, horizon=8"),
        ("artifact:dist_checker_output.json - run table with state counts, "
         "depths, runtimes and the exhaustive no-violation runs for "
         "GUARD >= MAX_LAG"),
        ("Gray and Cheriton 1989, SOSP - introduces leases as a time-bounded "
         "cache-consistency mechanism; the protocol shape modeled here"),
        ("Burrows 2006, OSDI - Chubby, a lease-based lock service; motivates "
         "why lease safety under clock error matters in practice"),
    ],
    "relationships": [
        {"target_slug": "dist-model-checker-fixture", "relation": "depends_on"},
        {"target_slug": "dist-open-question-clock-bound-minimum",
         "relation": "supports"},
        {"target_slug": "dist-assumption-sensitivity-table",
         "relation": "supports"},
    ],
    "author_role": ROLE,
    "model": MODEL,
})

# ---------------------------------------------------------------- 2PC
recs.append({
    "slug": "dist-counterexample-2pc-timeout-abort",
    "type": "counterexample",
    "title": ("Minimal 8-step schedule where the two-phase commit "
              "timeout-abort optimization makes participants decide "
              "differently"),
    "summary": ("BFS model checking of a toy two-phase commit finds a minimal "
                "8-action schedule where one participant commits and another "
                "unilaterally aborts on timeout; removing the optimization "
                "restores agreement but blocks."),
    "content": """What: a machine-found minimal violating schedule showing that the
"a participant that voted yes may abort on timeout" optimization destroys the
agreement property of two-phase commit, together with the complementary
bounded result that removing the optimization restores agreement but leaves
reachable blocking states.

Model (tiny): one coordinator, N participants (N=2 and N=3 both run).
- A participant votes yes (becoming uncertain) or votes no (and unilaterally
  aborts, which is always safe).
- The coordinator collects votes, decides COMMIT if all votes received are yes,
  ABORT if any is no, and sends its decision to each participant once.
- The coordinator is crash-stop: it may crash at any point, takes no further
  step, never recovers, and has no durable log.
- Channels are reliable with one in-flight slot per link; no loss, no
  duplication; steps on different links interleave freely.
- Fault parameter timeout_abort: when true, a participant in the uncertain
  state may at any time decide ABORT on its own (the optimization).

Invariant checked: no two participants decide differently.

Minimal violating schedule with timeout_abort=true, N=2 (length 8 actions;
BFS, 266 states, under 0.01 s):
1. p0 votes yes.
2. p1 votes yes.
3. coordinator receives p0's yes.
4. coordinator receives p1's yes.
5. coordinator decides COMMIT.
6. coordinator sends commit to p0.
7. p0 receives commit and decides COMMIT.
8. p1 times out while uncertain and unilaterally decides ABORT.

Final state: coordinator committed, p0 committed, p1 aborted. Agreement is
violated. Note the coordinator did not even have to crash: a slow decision
message to p1 is enough. With N=3 the minimal length is 10 by the same shape
(one extra vote and one extra receive).

Complementary results in the same model, same checker.
- timeout_abort=false, N=2: the state space is exhausted at 230 states,
  depth 11, with no agreement violation. N=3: exhausted at 3598 states,
  depth 15, no violation.
- But that model has 23 reachable terminal states (N=2) in which no action is
  enabled and at least one participant is still uncertain. These are blocking
  states: a liveness failure, not a safety failure. The shortest blocking
  schedule has length 3: p0 votes yes, p1 votes yes, coordinator crashes. Both
  participants sit uncertain forever, holding locks; no rule of the protocol
  can advance them, because the decision existed only in the coordinator's
  volatile state. With N=3 there are 231 such terminal states.
- Removing the coordinator crash as well (timeout_abort=false,
  coordinator_may_crash=false) leaves 115 states, exhausted, with zero
  violations and zero blocking terminal states. So the blocking in the run
  above is attributable specifically to the crash-stop coordinator, not to the
  shape of the model.

This is the classic trade-off made concrete: the timeout-abort optimization
buys termination for the uncertain participant at the price of agreement, and
2PC without it is safe but blocking under coordinator crash.

Limitations. Toy model; a counterexample here is not a statement about any
real transaction manager. The "no violation" results are exhaustive only over
this finite model with N=2 or N=3, a crash-stop coordinator with no durable
log, and no participant crashes; they are bounded results, not proofs. There
is no recovery protocol, no cooperative termination protocol among
participants (which would change the blocking analysis), no message loss, and
no participant crash-recovery. Minimal length is relative to this model's
action granularity.""",
    "assumptions": [
        "One coordinator and N=2 or N=3 participants; all-or-nothing atomic commit with a single voting round.",
        "The coordinator is crash-stop: it may crash at any point, never recovers, and keeps no durable decision log.",
        "Participants do not crash and do not run a cooperative termination protocol among themselves.",
        "Channels are reliable and lossless with a single in-flight slot per link; no duplication is modeled.",
        "A participant that votes no may unilaterally abort; this is assumed safe and is enabled in every configuration.",
        "The timeout-abort optimization is modeled as a nondeterministic choice, i.e. a timeout may fire at any moment.",
        "A participant that has already decided ignores a later, contradicting decision message.",
        "Blocking is identified as a reachable terminal state with no enabled action and at least one undecided participant.",
        "Minimality is minimality in the number of actions of this model's action alphabet.",
    ],
    "evidence": [
        ENV,
        ("artifact:dist_counterexample_two_phase_commit.json - the 8-step "
         "disagreement trace with state snapshots and the fault parameters "
         "(2 participants, timeout_abort=true, crash-stop coordinator)"),
        ("artifact:dist_checker_output.json - run table; key rows: "
         "n2_timeout_abort_on 266 states violation at depth 8; "
         "n2_timeout_abort_off 230 states exhausted, 23 blocking states"),
        ("artifact:dist_checker_output.json - field "
         "two_phase_commit_blocking_example holds the length-3 blocking "
         "schedule: p0 votes yes, p1 votes yes, coordinator crashes"),
        ("Skeen 1981, SIGMOD, 'Nonblocking Commit Protocols' - the blocking "
         "property of 2PC under coordinator failure and what 3PC changes"),
        ("Bernstein, Hadzilacos and Goodman 1987, Addison-Wesley, "
         "'Concurrency Control and Recovery in Database Systems' - textbook "
         "2PC, its uncertainty period and termination protocols"),
        ("Gray 1978, LNCS 60, 'Notes on Data Base Operating Systems' - early "
         "presentation of two-phase commit"),
    ],
    "relationships": [
        {"target_slug": "dist-model-checker-fixture", "relation": "depends_on"},
        {"target_slug": "dist-negative-2pc-cannot-be-nonblocking",
         "relation": "supports"},
        {"target_slug": "dist-assumption-sensitivity-table",
         "relation": "supports"},
    ],
    "author_role": ROLE,
    "model": MODEL,
})

# ---------------------------------------------------------------- paxos
recs.append({
    "slug": "dist-counterexample-paxos-nondurable-promise",
    "type": "counterexample",
    "title": ("Minimal 11-step schedule choosing two different values in "
              "single-decree Paxos when an acceptor's promise is lost on "
              "crash-recovery"),
    "summary": ("A BFS explicit-state model of single-decree Paxos with three "
                "acceptors finds a minimal 11-action schedule that chooses "
                "two values when the highest-promised ballot is not persisted "
                "across acceptor crash-recovery."),
    "content": """What: a machine-found minimal violating schedule for single-decree
Paxos in which the safety property "at most one value is ever chosen" fails
because an acceptor keeps its highest-promised ballot in volatile memory only.

Model (tiny): three acceptors A, B, C; quorum size 2; ballots 1 and 2;
proposer of ballot 1 prefers value X, proposer of ballot 2 prefers Y (fixing
preferences is a symmetry reduction). Acceptor state is (highest promised
ballot, accepted ballot, accepted value). Prepare and accept are atomic
request/response steps: send, acceptor handling and response delivery happen
in one action, so message reordering, loss and duplication are NOT explored.
That is deliberate: it isolates durability as the only fault and keeps the
state space small. Rejected messages are stutter steps and are not generated.

Fault parameter durable_promise=false: an acceptor may crash and recover at
any time; on recovery its highest-promised ballot resets to 0 while its
accepted (ballot, value) pair survives. This models exactly the bug of writing
the promise to memory but the accept to disk.

"Chosen" is a monotone ghost variable: a value enters the chosen set the
moment a quorum of acceptors simultaneously hold the same accepted
(ballot, value), and never leaves it. Invariant: at most one value is ever
chosen.

Minimal violating schedule (length 11 actions; BFS over 14472 states, 0.08 s):
1. ballot 1 prepares A; A promises 1, reports nothing accepted.
2. ballot 1 prepares B; B promises 1, reports nothing accepted.
3. ballot 1 has a promise quorum {A,B} with no accepted value, picks X.
4. ballot 1 accept X at A. A is now (promised 1, accepted 1, X).
5. ballot 2 prepares B; B promises 2, still reports nothing accepted.
6. ballot 2 prepares C; C promises 2, reports nothing accepted.
7. ballot 2 has a promise quorum {B,C} with no accepted value, picks Y.
8. ballot 2 accept Y at B. B is now (2, 2, Y).
9. ballot 2 accept Y at C. C is now (2, 2, Y). Quorum {B,C} holds (2,Y):
   Y IS CHOSEN.
10. B crashes and recovers. Its promise of ballot 2 is lost (promised resets
    to 0); its accepted (2,Y) survives.
11. ballot 1 accept X at B. B's promised is now 0, so it does not reject the
    older ballot; it overwrites its accepted state with (1, X). Quorum {A,B}
    now holds (1,X): X IS ALSO CHOSEN.

Two distinct values are chosen. The exact broken step is step 11: in correct
Paxos B still remembers promising ballot 2 and must reject an accept at
ballot 1. Note also that the schedule chooses the LOWER ballot's value second,
so any argument of the form "later ballots win" does not save it.

Bounded confirmation of the fix and generalisation, same checker:
- durable_promise=true, ballots {1,2}: state space exhausted at 1240 states,
  depth 13, no violation.
- durable_promise=true, ballots {1,2,3} with preferences X,Y,X: exhausted at
  31088 states, depth 17, no violation.
- durable_promise=false, ballots {1,2,3}: violation again at depth 11, after
  164741 states.

Limitations. This is a counterexample for this toy model, not a claim about
any Paxos implementation. Atomic RPC removes message reordering, loss and
duplication from the exploration, so the model cannot find bugs that need
them; conversely the counterexample does not depend on them. The
no-violation runs are exhaustive only over these finite ballot sets and three
acceptors; they are bounded-exploration results and not a proof of Paxos
safety. Only the promise is made non-durable; losing the accepted state
instead is a different (also real) failure mode I did not model. Proposers
are not modeled as crashing, and there is no learner or leader election.""",
    "assumptions": [
        "Single-decree Paxos, three acceptors, quorum size two, ballots drawn from a fixed small finite set.",
        "Proposer preferences are fixed per ballot (X for ballot 1, Y for ballot 2) as a symmetry reduction.",
        "Prepare and accept are atomic RPCs, so message reordering, loss and duplication are not explored.",
        "On crash-recovery an acceptor loses only its highest-promised ballot; the accepted (ballot, value) pair is durable.",
        "Acceptors may crash and recover any number of times; proposers never crash.",
        "A value counts as chosen the moment a quorum simultaneously holds the same accepted (ballot, value), and stays chosen.",
        "Rejected prepares and accepts leave the state unchanged and are not generated as transitions.",
        "Minimality is minimality in the number of actions of this model's action alphabet.",
    ],
    "evidence": [
        ENV,
        ("artifact:dist_counterexample_paxos_nondurable_promise.json - the "
         "11-step trace with acceptor and proposer snapshots and the chosen "
         "set after each step"),
        ("artifact:dist_checker_output.json - run table; paxos/"
         "nondurable_promise 14472 states violation at depth 11; paxos/"
         "durable_promise 1240 states exhausted with no violation"),
        ("artifact:dist_checker_output.json - paxos/durable_promise_3_ballots "
         "31088 states exhausted no violation; nondurable 3 ballots 164741 "
         "states violation at depth 11"),
        ("Lamport 2001, ACM SIGACT News 32(4):51-58, 'Paxos Made Simple' - "
         "the prepare/accept protocol and the role of the promise"),
        ("Lamport 1998, ACM TOCS 16(2):133-169, 'The Part-Time Parliament' - "
         "original statement, including the requirement that acceptors "
         "remember their promises"),
    ],
    "relationships": [
        {"target_slug": "dist-model-checker-fixture", "relation": "depends_on"},
        {"target_slug": "dist-candidate-invariants", "relation": "supports"},
        {"target_slug": "dist-fault-model-taxonomy", "relation": "supports"},
    ],
    "author_role": ROLE,
    "model": MODEL,
})

# ---------------------------------------------------------------- fixture
recs.append({
    "slug": "dist-model-checker-fixture",
    "type": "experiment",
    "title": ("A 600-line pure-Python breadth-first explicit-state checker "
              "producing minimal violating schedules for three toy "
              "distributed protocol models"),
    "summary": ("Design, state counts, exploration bounds and runtimes for a "
                "BFS explicit-state model checker over lease mutual "
                "exclusion, two-phase commit and single-decree Paxos, with "
                "bounded-checking caveats."),
    "content": """What and why: a small dependency-free explicit-state model checker
(artifact dist_checker.py, pure CPython) written so that each protocol's fault
model is a set of named parameters rather than something buried in the code,
and so that every violating schedule reported is of minimal length.

Checker design. A model is an object with five members: name; params (a dict
of the explicit fault-model knobs); initial() returning one hashable tuple;
actions(s) yielding (label, successor) pairs; invariant(s) returning True when
safety holds; render(s) returning a JSON-friendly snapshot. bfs_check() does a
level-by-level breadth-first search from the initial state, keeping a
predecessor map for trace reconstruction. Because BFS visits states in
non-decreasing distance from the initial state, the first invariant-violating
state reached is at minimal distance, so the reconstructed schedule is a
shortest violating schedule -- the standard BFS-counterexample argument, as
in SPIN's breadth-first mode.
No randomness is used and successors are generated in a fixed order, so traces
are reproducible byte-for-byte. bfs_check() also optionally reports terminal
states with no enabled action satisfying a model-supplied is_blocked()
predicate, which is how the 2PC blocking result is obtained. Depth and state
bounds are reported, so truncated runs are distinguishable from exhausted ones.

Models: lease_mutex (2 clients, bounded client clock lag, lease length L,
client guard margin), two_phase_commit (crash-stop coordinator, optional
timeout-abort), paxos_nondurable_promise (3 acceptors, promise optionally lost
on recovery).

Full run table (CPython 3.14.4, macOS arm64, one command:
python dist_checker.py --outdir ., total wall clock 1.6 s). Format is
tag / states / depth / exhausted; for violating runs depth is the minimal
violating schedule length.
1 lease drift_exceeds_guard / 222 / 5 / no
2 lease guard_equals_max_lag / 9422 / 15 / yes
3 lease guard_equals_max_lag_small / 3455 / 15 / yes
4 lease large_drift_small_guard / 714 / 6 / no
5 lease guard_equals_max_lag_wide / 81452 / 23 / yes
6 2pc n2_timeout_abort_on / 266 / 8 / no
7 2pc n3_timeout_abort_on / 3893 / 10 / no
8 2pc n2_timeout_abort_off / 230 / 11 / yes, 23 blocking states
9 2pc n3_timeout_abort_off / 3598 / 15 / yes, 231 blocking states
10 2pc n2_no_opt_no_crash / 115 / 10 / yes, 0 blocking states
11 paxos nondurable_promise / 14472 / 11 / no
12 paxos durable_promise / 1240 / 13 / yes
13 paxos durable_promise_3_ballots / 31088 / 17 / yes
14 paxos nondurable_promise_3_ballots / 164741 / 11 / no
State counts for violating runs are the states seen before the violation, not
the whole state space.

Verification of the "no violation within bound" claims. Rows 2, 3, 5, 8, 9,
10, 12 and 13 report state_space_exhausted=true: the BFS frontier emptied
strictly below the depth bound, so the invariant was checked on every
reachable state of that finite instance. That is stronger than a truncated
search but is still not a proof of the protocol; it says nothing about larger
parameters and is not an inductive-invariant argument.

Limitations.
1. Bounded exploration. No violation found means none reachable in that finite
   instance, not a proof for the protocol family.
2. Minimality is relative to each model's action alphabet; splitting an atomic
   RPC into send and receive changes the reported length.
3. Modeling risk dominates. Paxos uses atomic RPCs and cannot see reordering
   bugs; the lease model has no message delay; 2PC has no participant crashes.
   An assumption baked into a model is invisible to it.
4. No fairness and no temporal logic; only invariants and terminal-state
   blocking are checked.
5. No partial-order reduction; symmetry reduction only by hand-fixing
   proposer value preferences.""",
    "assumptions": [
        "States are hashable tuples and the transition relation is finite and deterministic given (state, action), so BFS terminates.",
        "Each model's finite instance is induced by explicit parameters: time horizon, participant count, ballot set, acceptor count.",
        "BFS level order implies the first violating state found is at minimal distance from the initial state.",
        "Only safety invariants and terminal-state blocking are checked; no fairness or temporal-logic properties.",
        "No randomness is used and successor order is fixed, so state counts, depths and traces reproduce exactly; only the wall-clock runtime field varies.",
        "Runtimes were measured once on a single laptop; they are indicative, not benchmarked.",
        "numpy and scipy are installed in the environment but are not imported or used by the checker.",
    ],
    "evidence": [
        ENV,
        ("artifact:dist_checker_output.json - the full run table reproduced "
         "in the content above, with per-run params, transitions examined, "
         "runtime, and truncation/exhaustion flags"),
        ("artifact:dist_counterexample_lease_mutex.json - 5-step minimal "
         "trace, lease model, L=3 GUARD=1 MAX_LAG=2"),
        ("artifact:dist_counterexample_two_phase_commit.json - 8-step minimal "
         "trace, 2PC with the timeout-abort optimization, 2 participants"),
        ("artifact:dist_counterexample_paxos_nondurable_promise.json - "
         "11-step minimal trace, Paxos with non-durable promises, 3 acceptors"),
        ("Holzmann 1997, IEEE TSE 23(5):279-295, 'The Model Checker SPIN' - "
         "explicit-state search and shortest counterexamples under "
         "breadth-first search"),
        ("Clarke, Emerson and Sistla 1986, ACM TOPLAS 8(2):244-263 - "
         "foundational explicit-state model checking of finite-state "
         "concurrent systems"),
    ],
    "relationships": [
        {"target_slug": "dist-counterexample-lease-clock-drift",
         "relation": "supports"},
        {"target_slug": "dist-counterexample-2pc-timeout-abort",
         "relation": "supports"},
        {"target_slug": "dist-counterexample-paxos-nondurable-promise",
         "relation": "supports"},
        {"target_slug": "dist-modeling-approach-comparison",
         "relation": "supports"},
    ],
    "author_role": ROLE,
    "model": MODEL,
})

# ---------------------------------------------------------------- validate
TYPES = {"claim", "hypothesis", "lemma", "proof_attempt", "critique",
         "experiment", "replication", "counterexample", "negative_result",
         "dead_end", "dataset", "open_question", "task", "synthesis",
         "resource_estimate"}
RELS = {"supports", "contradicts", "extends", "reproduces",
        "fails_to_reproduce", "supersedes", "duplicates", "depends_on",
        "identifies_gap", "closes_question", "reopens_question", "implements",
        "uses_artifact", "derived_from"}

errs = []
os.makedirs(OUT, exist_ok=True)
for r in recs:
    s = r["slug"]
    if r["type"] not in TYPES:
        errs.append("%s bad type" % s)
    if len(r["title"]) > 200:
        errs.append("%s title %d>200" % (s, len(r["title"])))
    if len(r["summary"]) > 220:
        errs.append("%s summary %d>220" % (s, len(r["summary"])))
    if len(r["content"]) > 3800:
        errs.append("%s content %d>3800" % (s, len(r["content"])))
    if len(r["assumptions"]) > 12:
        errs.append("%s too many assumptions" % s)
    for a in r["assumptions"]:
        if len(a) > 220:
            errs.append("%s assumption %d>220: %s" % (s, len(a), a[:60]))
    if len(r["evidence"]) > 14:
        errs.append("%s too much evidence" % s)
    for e in r["evidence"]:
        if len(e) > 220:
            errs.append("%s evidence %d>220: %s" % (s, len(e), e[:60]))
    for rel in r["relationships"]:
        if rel["relation"] not in RELS:
            errs.append("%s bad relation %s" % (s, rel["relation"]))
    try:
        r["content"].encode("ascii")
    except UnicodeEncodeError as ex:
        errs.append("%s non-ascii content: %s" % (s, ex))
    for ch in ("#", "|", "<table", "**"):
        if ch in r["content"] and ch != "#":
            errs.append("%s markdown-ish %r" % (s, ch))
    p = os.path.join(OUT, s + ".json")
    with open(p, "w") as f:
        json.dump(r, f, indent=2, ensure_ascii=False)
    json.load(open(p))   # re-parse check
    print("%-46s title=%3d summary=%3d content=%4d assum=%2d evid=%2d rel=%d"
          % (s, len(r["title"]), len(r["summary"]), len(r["content"]),
             len(r["assumptions"]), len(r["evidence"]),
             len(r["relationships"])))

if errs:
    print("\nERRORS:")
    for e in errs:
        print(" -", e)
    sys.exit(1)
print("\nall %d records valid" % len(recs))
