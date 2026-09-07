#!/usr/bin/env python3
"""dist_skeptic_probe.py -- independent review probes over dist_checker.py.

Written by the skeptical reviewer. Imports the shipped checker unchanged and
runs three probes:
  1. exhaustive sweep of the lease_mutex safety boundary (tests the shipped
     claim protocol_is_safe_iff = "guard >= max_lag");
  2. Paxos with a CUMULATIVE definition of "chosen" (a value is chosen if a
     quorum ever accepted (b,v)), replacing the simultaneous-snapshot ghost;
  3. Paxos with a non-durable promise but a recovery that repairs
     promised := accepted_ballot (restoring maxBal >= maxVBal for free).

Run:  python dist_skeptic_probe.py     (deterministic, no randomness, ~20 s)
"""
import sys, itertools
sys.path.insert(0, '/private/tmp/claude-501/-Users-benjamincox-dsm-haidaa/da69bd54-2977-4882-88ec-9babfab4e89d/scratchpad/commons/dist/artifacts')
from dist_checker import bfs_check, LeaseMutex, PaxosNondurablePromise, VAL_NAME

print("=== 1. lease boundary sweep: is 'safe iff guard >= max_lag' right? ===")
bad = []
for L in (2,3,4,5):
    for guard in range(0, L+1):
        for ml in range(0, 5):
            m = LeaseMutex(L=L, guard=guard, max_lag=ml, horizon=3*L)
            r = bfs_check(m, max_depth=200, max_states=500000)
            pred = guard < ml            # predicted unsafe
            obs = r["violation_found"]
            flag = "" if pred == obs else "  <<< MISMATCH"
            if pred != obs:
                bad.append((L,guard,ml,obs,r["state_space_exhausted"]))
            if flag:
                print("L=%d guard=%d max_lag=%d pred_unsafe=%s obs=%s exh=%s%s" % (L,guard,ml,pred,obs,r["state_space_exhausted"],flag))
print("mismatches:", len(bad))

print()
print("=== 2. Paxos with CUMULATIVE chosen (value chosen if a quorum EVER accepted (b,v)) ===")
class PaxosCumulative(PaxosNondurablePromise):
    def initial(self):
        acc = ((0,0,0),(0,0,0),(0,0,0))
        props = tuple((1,0,0,0) for _ in range(self.nprop))
        return (acc, props, ())          # third slot: sorted tuple of (b,v,acceptor) votes
    def _chosen_set(self, votes):
        cnt = {}
        for (b,v,a) in votes:
            cnt.setdefault((b,v), set()).add(a)
        return tuple(sorted({v for (b,v),s in cnt.items() if len(s) >= self.QUORUM}))
    def invariant(self, s):
        return len(self._chosen_set(s[2])) <= 1
    def _chosen_after(self, acc, chosen):
        return chosen                    # unused
    def actions(self, s):
        acc, props, votes = s
        out = []
        for label, s2 in PaxosNondurablePromise.actions(self, (acc, props, votes)):
            acc2, props2, _ = s2
            nv = votes
            for a in range(3):
                if acc2[a][1] != acc[a][1] or acc2[a][2] != acc[a][2]:
                    nv = tuple(sorted(set(nv) | {(acc2[a][1], acc2[a][2], a)}))
            out.append((label, (acc2, props2, nv)))
        return out
    def render(self, s):
        d = PaxosNondurablePromise.render(self, (s[0], s[1], ()))
        d["chosen_values"] = [VAL_NAME[v] for v in self._chosen_set(s[2])]
        d["votes"] = ["b%d/%s@%s" % (b, VAL_NAME[v], "ABC"[a]) for (b,v,a) in s[2]]
        return d

for dur in (False, True):
    r = bfs_check(PaxosCumulative(durable_promise=dur), max_depth=24)
    print("durable=%-5s states=%-7d depth=%-3d violation=%s at depth %s exhausted=%s"
          % (dur, r["states_explored"], r["max_depth_reached"], r["violation_found"],
             r["violation_depth"], r["state_space_exhausted"]))
    if r["violation_found"]:
        for st in r["violation_trace"]:
            print("   %2d %s" % (st["step"], st["action"]))

print()
print("=== 3. Paxos, non-durable promise but recovery repairs promised := accepted_ballot ===")
class PaxosRepair(PaxosNondurablePromise):
    def actions(self, s):
        acc, props, chosen = s
        out = [(l, s2) for (l, s2) in PaxosNondurablePromise.actions(self, s)
               if "CRASHES" not in l]
        for a in range(3):
            pr, ab, av = acc[a]
            if pr == ab:
                continue
            nacc = list(acc); nacc[a] = (ab, ab, av)
            out.append(("acceptor%s_CRASHES_and_RECOVERS (promised reset to accepted ballot %d)" % ("ABC"[a], ab),
                        (tuple(nacc), props, chosen)))
        return out
r = bfs_check(PaxosRepair(durable_promise=False), max_depth=24)
print("states=%d depth=%d violation=%s at depth %s exhausted=%s"
      % (r["states_explored"], r["max_depth_reached"], r["violation_found"],
         r["violation_depth"], r["state_space_exhausted"]))
if r["violation_found"]:
    for st in r["violation_trace"]:
        print("   %2d %s" % (st["step"], st["action"]))
