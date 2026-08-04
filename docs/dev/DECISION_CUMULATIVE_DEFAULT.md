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

### The threshold at the dispatcher is not the threshold at the entry point

Yielding from N = 2 would have been a regression. Below the crossover the cumulative scan's
near-constant cost -- its strict probe -- is not amortised, while hybrid is both accurate and
~20 ms per point, so a short scan would have paid several times over for precision it did not
need. Measured through `osc_prob_2nu_sun`, main against a branch that yielded from N = 2:

| N | 2 | 3 | 5 | 10 | 25 | 400 |
|---|---|---|---|---|---|---|
| main | 37 ms | 58 ms | 93 ms | 184 ms | 464 ms | 7459 ms |
| yield-from-2 | 283 ms | 259 ms | 269 ms | 282 ms | 263 ms | 292 ms |
| | **7.6x slower** | 4.5x | 2.9x | 1.5x | 1.8x faster | 25.6x faster |

Hence `HYBRID_YIELDS_TO_CUMULATIVE_MIN_POINTS = 25`, larger than `CUMULATIVE_AUTO_MIN_POINTS`.
The two guard different trades: at the entry point the alternative is the general per-point
path, which is *wrong* on solar profiles at small N (9.7e-3 at N = 10), so taking the cumulative
scan from N = 2 is right there; through the wrapper the alternative is hybrid, which is not.

Being the larger of the two also keeps the fall-through safe: when the dispatcher declines on
this count, `'auto'` is guaranteed to accept, so a scan can never be refused by both and land on
the general per-point method.

### `strategy='magnus'` opts out

That strategy is documented to reproduce the behaviour Magνs had before the adiabatic strategy
existed, *unconditionally*. The cumulative scan is Magnus machinery but postdates the promise
and builds a different grid, so the three wrappers now pass `cumulative=False` under it.
Without that, the escape hatch quietly stopped being one for exactly the case -- a baseline scan
-- where someone reproducing older numbers would reach for it. The test suite did not catch
this: the one test that pins `strategy='magnus'` uses a single point, where the cumulative scan
never engaged anyway.

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

Confirmed across the families the gate actually serves, not just 2ν standard potential -- all
three wrappers reach `_osc_prob_hybrid_dispatch` at every flavor count. Below the threshold the
errors are **bit-identical** to main, which is the evidence that nothing changed there:

| family | N = 5 | N = 25 | N = 200 |
|---|---|---|---|
| 2ν std | 90 → 92 ms, err identical | 448 → 70 ms, 2.2e-05 → 6.2e-07 | 3599 → 82 ms, 2.3e-05 → 2.2e-06 |
| 3ν std | 127 → 156 ms, err identical | 889 → 331 ms, 5.2e-05 → 4.2e-07 | 5141 → 332 ms, 9.3e-05 → 2.1e-06 |
| 2ν NSI | 106 → 129 ms, err identical | 529 → 91 ms, 2.4e-05 → 8.9e-07 | 4196 → 82 ms, 2.7e-05 → 1.6e-06 |
| 2ν LIV | 126 → 125 ms, err identical | 602 → 114 ms, 1.7e-05 → 3.1e-07 | 5355 → 168 ms, 2.4e-05 → 2.7e-07 |

## 4d. The accuracy sweep that found a regression, and the fix

A wide sweep through the wrapper -- 48 configurations, 1-100 MeV, N of 8/40/150, over 0.4 and
1.0 R_sun, each scored against `solve_ivp` -- looking for regressions rather than wins. At
`CUMULATIVE_N_ACC_SAFETY = 2` it found three, all at high energy over the shorter baseline:

| configuration | hybrid (main) | safety 2 | safety 4 |
|---|---|---|---|
| 60 MeV, N = 150, 0.4 R_sun | 1.57e-05 | **5.03e-05** | 8.56e-07 |
| 100 MeV, N = 150, 0.4 R_sun | 2.51e-05 | **3.77e-05** | 6.11e-07 |
| 100 MeV, N = 40, 0.4 R_sun | 9.13e-06 | **1.10e-05** | 5.58e-07 |

The diagnosis contradicted the obvious guess. The error is *not* at the short baselines where a
uniform grid would be expected to under-resolve: on the 60 MeV case it sits at the **longest**
baselines (5.03e-05 there against 2.18e-06 over the shortest third), and a probe at the short
end asks for a density within 1% of what the grid already provides. The scan was simply as
accurate as a tolerance-1e-3 Magnus grid, while hybrid happened to be better there.

So the lever is total resolution, and `CUMULATIVE_N_ACC_SAFETY` moved from 2 to 4. Two had been
chosen when the alternative was the general per-point path, against which the scan was already
124x faster and 11x more accurate; routing wrapper scans here makes the alternative *hybrid*,
which is a higher bar. Four removes all three regressions and beats hybrid on each, improves the
untouched configurations about twenty-fold as well, and costs ~1.4x on a path still tens of
times faster than what it replaces.

Re-run at safety 4, the same 48 configurations:

| | safety 2 | safety 4 |
|---|---|---|
| branch more accurate | 28 | **32** |
| effectively equal | 17 | 16 |
| **branch less accurate** | **3** | **0** |
| worst error (main 3.21e-05) | 5.03e-05 | **2.93e-05** |
| outside the requested 1e-3 | 0 | 0 |

And the threshold still holds at the higher cost -- across N from 2 to 400 at 5 and 10 MeV,
**no configuration is more than 10% slower**, errors below the threshold are bit-identical to
main, and at N = 25 the scan is already 1.3x faster with 190x-2800x better accuracy.

**Unitarity loosens, from 8.9e-16 to 8.2e-12.** The hybrid strategy returns an exactly unitary
operator by construction; the cumulative scan accumulates a product over tens of thousands of
slabs, so it is unitary to roughly 1e-11 rather than to machine precision. Far inside the 1e-9
the suite asserts, and inherent to a long product rather than fixable by tuning -- but it is a
change in kind, not only in magnitude, and worth knowing before anyone tightens a unitarity
assertion.

## 4e. The probe must not warn about answers it throws away

An ordinary wrapper baseline scan began emitting `MagnusConvergenceWarning` where `main`
(hybrid) was silent. Tracing it: the warning comes from the **probe**, not the traversal --

```
probe alone      : ['MagnusConvergenceWarning']
traversal alone  : (none)
```

-- and the probe keeps only a slab count; its probabilities are discarded. So the warning
described a result nobody receives, while the grid it sizes produces no such warning when
actually walked. That is misleading rather than merely noisy, so it is suppressed for that one
call, for that one category. Verified narrow in both directions: a scan that genuinely cannot
meet its tolerance still raises `ToleranceNotAchievedWarning`, and a deliberately coarse call
elsewhere still raises `MagnusConvergenceWarning`.

Re-running the notebooks afterwards: **no figure changed** (a warning cannot move a number), and
warnings disappeared from exactly the cumulative cells -- notebook 02 cells 64 and 92, notebook
03 cells 57, 71 and 99 -- with none gained. This also retires the cost recorded earlier on this
branch, that notebook 03 cell 57 had started warning where it previously did not.

**A related fact, documented rather than fixed.** Over a full solar radius at 5 and 10 MeV the
strict probe reaches `max_n_slabs` (20 000) without two successive levels agreeing, so `n_acc`
is derived from a ceiling rather than from a converged requirement. The resulting scans measure
~5e-08 against `solve_ivp`, so the safety factor is carrying that; but a caller who lowers
`max_n_slabs` lowers the scan's resolution in proportion. Recorded in `CUMULATIVE_N_ACC_SAFETY`.

**Amended 2026-08-03** (`FINDINGS_ADVERSARIAL_VALIDATION.md` §4). This paragraph originally ended
"...in proportion, **without a separate warning**", which overstates the exposure. Measured on a
60-point solar scan at 10 MeV over a full solar radius:

| `max_n_slabs` | probe `n_slabs` | scan error | warned? |
|---|---|---|---|
| 500 | 500 | 6.194e-02 | `ToleranceNotAchievedWarning` |
| 2000 | 2000 | 4.657e-03 | `ToleranceNotAchievedWarning` |
| 5000 | 5000 | 3.053e-06 | `ToleranceNotAchievedWarning` |
| 20000 (default) | 20000 | 9.546e-09 | — |

The degradation is real and proportional, as stated. But it is **not silent**: there is no
*cumulative-specific* warning, and the probe's `MagnusConvergenceWarning` is suppressed here
(§4e above), yet `ToleranceNotAchievedWarning` survives that suppression and reaches the caller
at every capped level. That is the signal that matters, and it is exactly what §4e's narrowing
was designed to preserve.

## 4f. What two further adversarial passes established

Aimed at surfaces the earlier passes could not reach: parameter combinations the new path had
never run on, and numerical oracles other than a single `solve_ivp` comparison.

**Dimensions verified, all previously untested on this path.** Channel selection (`nu_i`/`nu_f`
matches the full matrix exactly for all four channels), antineutrinos (6.84e-09 against
`solve_ivp`), 4ν and 5ν scans (unitary to 1.8e-11 and 1.1e-11), Magnus orders 2/3/6
(3.75e-07 to 2.13e-09), Python lists rather than arrays (bit-identical), and a **non-zero
`L0`** -- a scan starting at 0.2 R_sun matches `solve_ivp` to 1.24e-07, so the traversal's
origin handling is right.

**Quadrature methods.** `trapezoid` and `simpson` cannot reach the default tolerance on a full
solar radius by *any* route: a single adaptive call at the far end gives 1.0e-01, having
exhausted both `max_n_slabs` (2000) and `max_n_tpts_per_slab` (500). Against that, the
cumulative scan is a large improvement rather than a regression -- 5.19e-03 where the per-point
path it replaces gives 1.00e-01. Recorded because the first reading of this measurement looked
like a regression and was not.

**Tighter tolerances are delivered.** Requested 1e-3, 1e-5 and 1e-7 all return ~1e-08, at a flat
~368 ms: the probe is already capped, so a tighter request costs nothing and is already
exceeded. Beyond what the cap can deliver the caller is still told -- at 1e-9 and 1e-12 the scan
raises `ToleranceNotAchievedWarning`, confirming that suppressing the probe's convergence
warning (§4e) did not swallow the signal that matters. On the same request the new path is 500x
more accurate than the old one (1.13e-08 against 5.77e-06).

**Independent oracles.** Against `expm` for a constant Hamiltonian -- the time-ordering and
indexing check, not an accuracy one -- the cumulative scan agrees to **2.31e-14** over 60
baselines. At extreme low energy, where the accumulated phase is largest, 0.5 MeV gives 2.74e-06
and 0.2 MeV 5.49e-06.

**One change was reverted for lack of evidence.** The traversal takes `n_tpts_per_slab` from a
fixed default rather than from the probe, which looked like an oversight for the quadrature
methods, where accuracy depends on it jointly with the slab count. Inheriting it was implemented
and then measured across eight configurations: the probe reports the 500 cap every time, and the
answer is identical either way, because the error is dominated by the slab count. The change was
removed rather than kept on principle -- unmeasured complexity in a dispatch path is what this
document exists to argue against.

## 4g. Multiple resonances, the hardest case for this dispatch choice

The two candidate paths fail in opposite ways here, so a profile with many non-adiabatic
crossings is the sharpest available test. The cumulative scan has **no resonance detection at
all** -- one uniform grid, sized by a single probe. The hybrid strategy locates each window with
a Hellmann-Feynman diagnostic and patches it. The obvious expectation is that hybrid wins.

Profile: 3ν with NSI `eps_et = 3.0`, on an exponential solar decay modulated by a strong sine so
the resonance density is crossed repeatedly. `hybrid_propagator` reports **ten windows** across
0.5–7 l_scale. Baseline scan of N = 60, everything scored against `solve_ivp`/DOP853:

| path | time | max error |
|---|---|---|
| **cumulative (the new default)** | **468 ms** | **1.02e-05** |
| per-point general | 7.8 s | 4.00e-02 |
| hybrid, one call per baseline | 101 s | **2.86e-01** |

Through the user-facing wrapper, `main` against this branch:

| N | main | branch |
|---|---|---|
| 8 (below the threshold) | 7.48e-02 | 7.48e-02 — unchanged, hybrid keeps it |
| 60 (above it) | **2.86e-01** | **1.02e-05** |

**A 28 000× improvement**, and the expectation was exactly backwards: the path with no resonance
detection is the one that gets it right, because a grid fine enough for the whole trajectory
resolves the crossings without needing to find them.

### The reason hybrid loses is worse than being inaccurate

It **certifies while wrong**. Per baseline, against `solve_ivp`:

| L (l_scale) | windows found | certified | error |
|---|---|---|---|
| 1.54 | **0** | yes | 4.32e-02 |
| 3.36 | 4 | yes | 6.48e-02 |
| 5.18 | 6 | yes | 7.48e-02 |
| 7.00 | 10 | yes | 7.40e-02 |

Eight of eight sampled baselines certified, every one wrong by far more than the requested 1e-3.
That it reports **zero** windows at 1.54 l_scale while being wrong by 4.3e-02 points at the
detector missing crossings rather than mis-patching them. The certification is self-referential
-- successive refinements of its own patches, never against truth -- so it can converge to its
own wrong answer, which is the failure class `NOTES_ADAPTIVE_REFINEMENT.md` documents on the
Magnus ladder, here on the adiabatic path.

This is **pre-existing and not introduced by this branch**: on `main` every baseline scan went to
hybrid, so `main` returns 2.86e-01 at N = 60. The routing here rescues scans of 25 points or
more; single points and shorter scans still reach hybrid and are still exposed. Logged as
separate work rather than fixed here, since the defect is in `magnus.adiabatic`, not in the
dispatch.

**Resolved 2026-08-03**, in `magnus.adiabatic` as this section anticipated. The mechanism was
that a result with *no* window is certified on the strength of two adiabatic answers agreeing
with each other — which they always do, since successive iterations differ only in the transport
grid. `hybrid_propagator` now additionally requires the adiabaticity parameter itself to be small
enough for the requested tolerance before an empty window list may be certified (see
`GAMMA_TO_ERROR`), and lowers the threshold until a window opens otherwise. The claim in
`docs/source/adiabatic_strategy.rst` has been corrected rather than merely flagged: the per-pair
half of it is true and now verified at 2, 3, 4 and 5 flavors, while the two genuine limits — a
profile that is not smooth at the probe scale, and a feature narrower than the probe spacing —
are stated there explicitly. See `FINDINGS_ADVERSARIAL_VALIDATION.md` §3.2 and §6.

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
| solar baseline scan via the wrappers, N >= 25 | **23–31× faster**, 1e-05 → 1e-06, and 5.2e-03 → 1.0e-06 at 10 MeV (§4c) |
| solar baseline scan via the wrappers, N < 25 | unchanged — hybrid keeps it (§4c) |
| accuracy across 48 wrapper configurations | **32 better, 16 equal, 0 worse**; unitarity 9e-16 → 8e-12 (§4d) |
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
