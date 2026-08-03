# Decision: `cumulative` defaults to `'auto'`

**Written:** 2026-08-03, executing the flip that `DECISION_OSCPROB_CUMULATIVE.md` §8 called
"the obvious next question".
**Machine:** the same 12-core box; timings are min-of-2.
**Status: SHIPPED**, as `osc_prob_energy_baseline(..., cumulative='auto')`.

---

## 1. The call

**Flip it — but to `'auto'`, not to `True`.**

`cumulative=True` does not fall back; it **raises** when the request does not fit (differing
energies, `t_slab_edges`, a baseline behind `L0`). Flipping the default to `True` would
therefore break every multi-energy call in the package. That is not a detail to work around: it
is the right design for an explicit request — a caller who names the cumulative scan should hear
that it did not happen — so the default becomes a third value that resolves per call.

## 2. Why flip at all: it is more accurate, not merely faster

Per-point versus cumulative, scored against `solve_ivp`/DOP853, sampled across each scan:

| N | solar 2ν, 5 MeV → R_sun | | 3ν exponential (nb03/43) | |
|---|---|---|---|---|
| | time | error, per-point → cumulative | time | error, per-point → cumulative |
| 2 | 0.75× | 2.3e-05 → 9.4e-07 | 0.84× | 1.7e-04 → 6.7e-07 |
| 10 | 0.87× | **9.7e-03** → 3.3e-06 | 2.22× | 5.2e-04 → 7.4e-06 |
| 25 | 2.65× | **5.6e-03** → 7.2e-06 | 4.63× | 4.5e-04 → 1.1e-05 |
| 100 | 7.01× | **2.6e-03** → 8.3e-06 | 14.5× | 5.3e-04 → 4.1e-06 |
| 500 | 34.3× | 1.7e-04 → 7.0e-06 | 35.5× | 4.9e-04 → 4.6e-06 |
| 1000 | 84.3× | 4.1e-04 → 3.3e-06 | 44.1× | 5.9e-04 → 4.5e-06 |
| 3000 | 158× | 6.7e-04 → 9.5e-06 | 52.7× | 5.0e-04 → 4.5e-06 |

The cumulative scan is more accurate **at every N measured**, usually by two or three orders.
And on the solar scan the per-point path returns answers **outside the requested 1e-3** at
N = 10, 25 and 100 — the same silent-tolerance failure recorded in
`NOTES_ADAPTIVE_REFINEMENT.md` §4b, here on the baseline axis. That makes this a correctness
change with a speed consequence, not an optimisation.

## 3. Where `'auto'` engages, and the two exclusions that are not obvious

`'auto'` resolves to the cumulative scan when **all** of:

1. the Hamiltonian varies with position;
2. every requested energy is equal;
3. no `t_slab_edges` were given;
4. every baseline is at or beyond `L0`;
5. there are at least `CUMULATIVE_AUTO_MIN_POINTS` (= 2) of them.

Conditions 2–4 are what `cumulative=True` already checks. The other two were added here.

**At least two baselines.** A single point has no prefix to reuse and would pay the
inherited-grid probe for nothing. This matters more than it sounds: every *single-point* call
through the wrapper layer — `osc_prob_2nu_sun(E, L, ...)` and its siblings — is served by
`osc_prob_energy_baseline`, so without this the flip would roughly double the cost of the
package's most ordinary call.

**A position-dependent Hamiltonian.** For vacuum or constant density, `osc_prob` detects the
constant and integrates it exactly on one slab. There is no traversal to share, and the
cumulative scan would replace an exact one-slab answer with an adaptive probe plus a walk —
strictly worse on both axes.

**This one was found by the docs build, not by the test suite.** With the rule as first written,
a three-baseline vacuum example in the documentation
(`osc_prob_3nu_vacuum(energy, baselines, ...)`) began emitting `MagnusConvergenceWarning`, which
`sphinx -W` escalates to an error. The full 665-test suite passed at that point. Two lessons
worth carrying: `-W` on the docs is a real test of user-facing behaviour and caught a case the
suite did not; and the fix was found by instrumenting the engagement path and rebuilding, rather
than by guessing which example was responsible — three guesses in a row were wrong.

## 4. The threshold is deliberately below the speed crossover

Cumulative becomes *faster* at N ≈ 25 (solar) and N ≈ 5 (3ν exponential), not at N = 2. Between
2 and the crossover it is at most ~1.3× slower — a few milliseconds — while being one to three
orders more accurate, and while removing the tolerance misses in §2. Trading milliseconds for
that is the right way round, so the floor sits at "any genuine scan" rather than at the
crossover. `CUMULATIVE_AUTO_MIN_POINTS` names it, and the reasoning is recorded there.

## 4b. Its reach is narrower than it looks, because the reorder gets there first

Measured after the fact, alternating between `main` and the branch to cancel machine drift
(a whole-suite comparison cannot: `test_adiabatic.py` moved 124 s → 69 s between two runs of
code neither state modifies):

| case | main | branch |
|---|---|---|
| single-point solar wrapper | 27–208 ms | 28–33 ms |
| **solar baseline scan through the wrapper, N = 400** | **8.9–11.5 s** | **9.7–14.8 s** |
| vacuum baseline scan, N = 300 | 43–46 ms | 34–97 ms |
| constant-density scan, N = 300 | 43–47 ms | 34–131 ms |

The scan the flip exists for does not move. Spying on the dispatch chain shows why:

```
through osc_prob_2nu_sun:                hybrid=1, cumulative=0, osc_prob_energy_baseline never reached
through osc_prob_energy_baseline direct: cumulative=1
```

Since `DECISION_DISPATCH_ORDER.md` put the hybrid strategy first, it answers a solar **baseline
scan** too — one independent `hybrid_propagator` call per point, at its ~28 ms floor — and
returns before `osc_prob_energy_baseline` is ever called. So on the wrapper families
(`osc_prob_*_sun` and siblings) this default changes nothing at all.

Where it does apply: direct `osc_prob_energy_baseline` callers, and any wrapper request the
dispatchers decline. The §2 table was measured through the direct entry point, so those numbers
stand for that path — they are simply not what a `osc_prob_2nu_sun` user gets.

**This is a missed win, not just a scoping note.** At N = 400 the wrapper spends ~11 s where the
cumulative scan answers in ~0.3 s (§2, N = 500: 299 ms) and more accurately. Teaching the
dispatch chain to prefer the cumulative scan for a single-energy baseline scan — before handing
the points to hybrid one at a time — is the obvious follow-up, and it is a change to the
dispatchers rather than to this default, so it is deliberately not made here.

## 5. Cost, stated plainly

**Results move.** The two paths build different grids, so any applicable baseline scan returns
slightly different numbers — within the requested tolerance, and on all evidence toward the
truth. `cumulative=False` reproduces the previous behaviour exactly, and is tested to.

**Silently.** No warning when `'auto'` engages. This matches how the package already treats its
other dispatch choices — the hybrid strategy, the interaction-picture fast path and the batched
scan engine all engage silently — and a warning on this path would fire on very ordinary calls.
The change belongs in the changelog, which is where a user comparing against stored results will
look.

**No test broke.** The full suite passed unchanged across the flip, including the tests that
compare batched scans against per-point ones, which is evidence that the movement really is
inside tolerance.

## 6. Honest summary of the runtime effect

| | before → after |
|---|---|
| single-point wrapper call | unchanged (excluded by `CUMULATIVE_AUTO_MIN_POINTS`) |
| vacuum / constant-density scan | unchanged (excluded: position-independent `H`) |
| solar baseline scan via the wrappers | **unchanged** — hybrid answers first (§4b) |
| baseline scan via `osc_prob_energy_baseline` | 2.65× at N = 25, 84× at N = 1000, and 1–3 orders more accurate at every N (§2) |
| notebooks 02 and 03 | unchanged — their converted cells pass `cumulative=True` explicitly |

So the flip is a correctness-and-speed win **for direct callers of the primordial baseline-scan
entry point**, and a no-op everywhere else until the dispatch chain is taught about it.

## 7. What would change this call

- **A caller who needs bit-reproducibility against pre-1.0.0 results** and does not know about
  `cumulative=False`. The changelog is the mitigation; if that proves insufficient, the
  fallback is a one-session warning on first engagement.
- **A profile where the inherited `n_acc` is badly wrong.** The grid is sized from one adaptive
  `osc_prob` call at the longest baseline (times `CUMULATIVE_N_ACC_SAFETY`), so it inherits that
  path's failure modes — including the coincidental-agreement stop that `strict_convergence`
  exists to reject. Nothing measured here shows it, but the coupling is real and worth
  remembering.
