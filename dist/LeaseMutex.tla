---- MODULE LeaseMutex ----
(* UNTESTED. Never run through TLC or the TLA+ parser (SANY). Written by hand as a
   text artifact for a modeling-approach comparison; it is a plausible-looking
   translation of the Python lease_mutex explicit-state model (server grants a
   lease of length L; client believes it holds the lock while its own local
   clock is below grant_time + L; server reclaims once its own clock reaches
   grant_time + L) into TLA+. Treat every line below as a hypothesis, not a
   verified fact. Known-likely problems, left in deliberately rather than
   silently fixed, are listed at the bottom instead of in-line, so a reader
   can compare "what I intended" against "what is probably wrong". This is
   exactly the artifact that dist-task-run-tlc-lease is supposed to run
   through SANY/TLC and correct. *)

EXTENDS Integers, FiniteSets

CONSTANTS
    Clients,        \* nonempty finite set of client ids
    L,               \* lease length, a positive integer in server-clock units
    MaxDrift,        \* max extra ticks a client clock may advance per server tick (>= 0)
    MaxSkew,         \* max |client_clock - server_clock| allowed at Init
    MaxTime          \* bound on server_clock, for finite-state exploration only;
                     \* not part of the protocol, purely a TLC search bound

ASSUME /\ L \in Nat \ {0}
       /\ MaxDrift \in Nat
       /\ MaxSkew \in Nat
       /\ MaxTime \in Nat
       /\ Clients # {}

VARIABLES
    server_clock,     \* Nat: server's own monotonic clock
    client_clock,      \* [Clients -> Int]: each client's local clock
    holder,            \* server's view: Clients \cup {"none"}, who currently holds the lease
    grant_time,        \* Int: server_clock at the most recent grant; -1 if none granted
    believes_holds      \* [Clients -> BOOLEAN]: each client's own belief that it holds the lock

vars == <<server_clock, client_clock, holder, grant_time, believes_holds>>

TypeOK ==
    /\ server_clock \in 0..MaxTime
    /\ client_clock \in [Clients -> (-MaxSkew)..(MaxTime + MaxSkew + MaxDrift)]
    /\ holder \in Clients \cup {"none"}
    /\ grant_time \in (-1)..MaxTime
    /\ believes_holds \in [Clients -> BOOLEAN]

Init ==
    /\ server_clock = 0
    /\ client_clock \in [Clients -> (-MaxSkew)..MaxSkew]
    /\ holder = "none"
    /\ grant_time = -1
    /\ believes_holds = [c \in Clients |-> FALSE]

(* Server grants the lease to an idle client. Simplification: the grant
   message is instantaneous and grant_time is shared state legible to both
   server and client immediately -- i.e. this model does NOT represent
   message delay on the grant path itself, only clock drift afterward. That
   is itself a modeling assumption worth flagging, not a fact about real
   lease services. *)
Grant(c) ==
    /\ holder = "none"
    /\ holder' = c
    /\ grant_time' = server_clock
    /\ believes_holds' = [believes_holds EXCEPT ![c] = TRUE]
    /\ UNCHANGED <<server_clock, client_clock>>

(* Server reclaims once ITS clock says the lease has expired. This does not
   by itself change any client's belief -- that only happens via
   ClientExpire below, and only when the CLIENT's own (possibly drifted)
   clock agrees. The gap between Reclaim and the corresponding ClientExpire
   is exactly where the safety violation is expected to live when a client
   clock runs slow relative to the server. *)
Reclaim ==
    /\ holder # "none"
    /\ server_clock >= grant_time + L
    /\ holder' = "none"
    /\ grant_time' = -1
    /\ UNCHANGED <<server_clock, client_clock, believes_holds>>

ClientExpire(c) ==
    /\ believes_holds[c] = TRUE
    /\ grant_time >= 0
    /\ client_clock[c] >= grant_time + L
    /\ believes_holds' = [believes_holds EXCEPT ![c] = FALSE]
    /\ UNCHANGED <<server_clock, client_clock, holder, grant_time>>

ServerTick ==
    /\ server_clock < MaxTime
    /\ server_clock' = server_clock + 1
    /\ UNCHANGED <<client_clock, holder, grant_time, believes_holds>>

(* Client clock advances by at least 1 and at most 1 + MaxDrift per step, so
   MaxDrift bounds how much FASTER than the server a client can run. Slow
   clients are represented implicitly: because ClientTick is not required to
   fire on every ServerTick (this is an interleaving semantics, not a
   synchronous product), a client that is simply not scheduled falls behind
   the server by an unbounded number of ticks. That asymmetry -- explicit
   bound on fast drift, only an implicit/unbounded bound on slow drift via
   scheduling starvation -- is almost certainly not what the charter's
   "maximum drift ratio" parameter was meant to capture, and is flagged here
   as a probable modeling bug rather than corrected silently. A fix likely
   needs a symmetric per-tick bound, e.g. requiring client and server ticks
   to be interleaved via a joint Tick action with delta in
   -MaxDrift..MaxDrift, or an explicit fairness/weak-fairness constraint so
   TLC's bounded search cannot hide behind starvation alone. *)
ClientTick(c) ==
    /\ client_clock[c] < MaxTime + MaxSkew + MaxDrift
    /\ \E delta \in 0..MaxDrift :
          client_clock' = [client_clock EXCEPT ![c] = @ + 1 + delta]
    /\ UNCHANGED <<server_clock, holder, grant_time, believes_holds>>

Next ==
    \/ \E c \in Clients : Grant(c)
    \/ Reclaim
    \/ \E c \in Clients : ClientExpire(c)
    \/ ServerTick
    \/ \E c \in Clients : ClientTick(c)

Spec == Init /\ [][Next]_vars

(* Safety property this whole model exists to check: at most one client
   simultaneously believes it holds the lease. *)
AtMostOneHolder ==
    Cardinality({c \in Clients : believes_holds[c] = TRUE}) <= 1

====

(* Known-likely problems (untested; listed instead of silently fixed):
   1. Drift is modeled asymmetrically -- MaxDrift bounds only how much a
      client can run FAST relative to the server (via the delta in
      ClientTick); "slow" drift is represented only implicitly through
      TLC/the interleaving scheduler choosing not to fire ClientTick(c),
      which is unbounded, not the bounded skew the charter parameter
      describes. This likely needs a joint tick action or a fairness
      constraint (or both) before the model matches the intended semantics,
      and may currently either fail to reproduce the Python drift
      counterexample or produce a different, weaker one via starvation.
   2. No .cfg file is included. TLC needs CONSTANTS bound to concrete
      values (e.g. Clients = {c1, c2}, L = 3, MaxDrift = 1, MaxSkew = 1,
      MaxTime = 10) and AtMostOneHolder listed as an INVARIANT; none of
      that is provided here.
   3. Grant's simplification (grant_time is instantaneously shared state,
      not a delivered message) may itself remove the interesting behavior
      -- if the real bug requires a stale client belief formed from a
      grant_time the client observed before some intervening event, this
      model cannot represent that, since there is only one grant_time
      variable rather than a per-client recorded value.
   4. SANY parse errors (stray commas, module name vs. file name mismatch,
      EXTENDS list) have not been checked at all; treat the syntax itself
      as unverified.
   This list is deliberately not exhaustive; running SANY/TLC will likely
   surface more. *)
