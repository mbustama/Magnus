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

**This was a missed win, and it is now taken** — see §4c.

## 4c. Reaching the wrapper layer, and the trap on the way

`_osc_prob_hybrid_dispatch` now returns `NotImplemented` for a single-energy baseline scan when
`strategy == 'auto'`, so the caller falls through to `osc_prob_energy_baseline` and
`cumulative='auto'` engages. Declining there is sufficient: `ip_exp` requires every baseline
equal and the separable engine a single shared baseline, so both decline a scan as well.
`strategy='hybrid'` is an explicit request and still gets hybrid; a single point still gets
hybrid, since a scan of one has nothing to reuse. Both are tested.

The justification, on solar profiles against `solve_ivp`, hybrid versus cumulative:

| E, N | hybrid | cumulative | |
|---|---|---|---|
| 5 MeV, 50 | 999 ms, 1.06e-05 | 170 ms, 5.36e-06 | 5.9× |
| 5 MeV, 800 | 17 635 ms, 1.01e-05 | 206 ms, 4.03e-06 | 85.6× |
| 20 MeV, 200 | 8 352 ms, 3.01e-05 | 220 ms, 8.37e-07 | 38.0× |

**The trap.** Making hybrid stand aside on that evidence alone would have shipped a regression.
At **10 MeV** — the configuration where the adaptive ladder stops on a coincidental agreement —
the cumulative scan was **outside the requested tolerance at every scan size**:

| 10 MeV | N = 2 | N = 5 | N = 50 | N = 400 |
|---|---|---|---|---|
| hybrid | 2.21e-07 | 1.42e-05 | 1.21e-05 | 1.66e-05 |
| cumulative, loose probe | **3.10e-03** | **2.87e-03** | **6.68e-03** | **5.20e-03** |

The cause is the coupling recorded in §7: the grid is sized by one ordinary adaptive `osc_prob`
probe at the longest baseline, and at 10 MeV that probe returns `n_slabs = 3298` on a
coincidence — 6596 after the safety factor, against the ~20 000 the case needs. Handing the
whole scan to a grid built on one early-stopping call turns a single bad point into N of them.

**The fix is `strict_convergence` on the probe**, which is now always applied there whatever the
caller asked for their own points. It is the one call whose convergence decides the entire
scan's grid, and its cost — one extra refinement level on a single call — amortises over every
baseline:

| E, N | loose probe | strict probe |
|---|---|---|
| 5 MeV, 400 | 1.5e-05 | 2.35e-06 |
| **10 MeV, 400** | **5.20e-03** | **1.01e-06** |
| 20 MeV, 400 | 3.3e-06 | 8.12e-07 |

That is the feature this project shelved as too costly to enable by default, earning its place
in the one spot where the economics invert.

### Measured through the wrapper, alternating between `main` and the branch

| case | main | branch |
|---|---|---|
| single-point solar wrapper | 20.5 / 40.8 / 25.7 ms | 26.5 / 25.5 / 27.9 ms |
| **solar baseline scan, N = 400** | **7632 / 10178 / 10086 ms** | **328 / 325 / 328 ms** |
| vacuum baseline scan, N = 300 | 37.9 / 38.2 / 37.9 ms | 52.4 / 37.7 / 40.3 ms |
| constant-density scan, N = 300 | 39.4 / 37.9 / 38.3 ms | 54.6 / 37.7 / 39.2 ms |

**23–31× on the baseline scan**, reproducible across all three rounds, with accuracy improving
from ~1e-05 to ~1e-06 and the 10 MeV case from 5.2e-03 to 1.0e-06. Everything else is unchanged
within scatter.

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
| solar baseline scan via the wrappers | **23–31× faster**, 1e-05 → 1e-06, and 5.2e-03 → 1.0e-06 at 10 MeV (§4c) |
| baseline scan via `osc_prob_energy_baseline` | 2.65× at N = 25, 84× at N = 1000, and 1–3 orders more accurate at every N (§2) |
| notebooks 02 and 03 | unchanged — their converted cells pass `cumulative=True` explicitly |

So the flip is a correctness-and-speed win on both the primordial entry point and the wrapper
families, and a deliberate no-op for single points and position-independent Hamiltonians.

## 7. What would change this call

- **A caller who needs bit-reproducibility against pre-1.0.0 results** and does not know about
  `cumulative=False`. The changelog is the mitigation; if that proves insufficient, the
  fallback is a one-session warning on first engagement.
- **A profile where the inherited `n_acc` is badly wrong.** The grid is sized from one adaptive
  `osc_prob` call at the longest baseline (times `CUMULATIVE_N_ACC_SAFETY`), so it inherits that
  path's failure modes — including the coincidental-agreement stop that `strict_convergence`
  exists to reject. Nothing measured here shows it, but the coupling is real and worth
  remembering.
