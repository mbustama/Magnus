# Handover: the supernova shock silent band  — CLOSED

> **STATUS, 2026-08-05: every item this document listed as open has landed.**
> `dev-robustness` (PR #28), `routing-regression-test` (PR #29) and `example-notebooks`
> (PR #30) are all **merged into `main`**. This file is kept for the reasoning and the
> measured numbers, which remain valid; its task list does not. See "What is left" below
> for the current state.
>
> The document was originally written about a **single-point dispatch exposure**, and
> §13.17-13.20 of the findings retracted most of that: applying the *observable* (the
> phase-averaged probability) collapses the solar error 53x and leaves no silent miss
> anywhere on the physical population.

**Written:** 2026-08-05, at the close of the session that executed
`HANDOVER_PHYSICAL_PROFILES.md`; reconciled the same day after PRs #28-#30 merged.

---

## 0. Verify the base before starting

```bash
git -C ~/Research/magnus log --oneline -1        # expect the PR #30 merge or later
python -c "import sys; sys.path.insert(0,'src'); import magnus.oscprob as o; \
  print(o.HYBRID_YIELDS_TO_CUMULATIVE_MIN_POINTS, o.CUMULATIVE_AUTO_MIN_POINTS)"
# must print: 8 2
```

If that prints `1 1`, someone has already applied the change this document argues against
shipping as-is. Read §3 before going further.

**Read first:** `FINDINGS_ROBUSTNESS_PROGRAMME.md` §13, then this.

---

## 1. The job, in one paragraph

**Read §13.17-13.19 of the findings first - they retract most of what the rest of this document
was written to support.** Applying the *averaged* observable, which is what solar and supernova
physics actually measures, collapses the solar error by 53x and the turbulence error by 23x,
leaving both inside tolerance and not user-facing. Only a **supernova shock** survives averaging,
and on sharp shocks the package already warns every time. What is left is one narrow, genuinely
open question was **the silent band around a 70 km shock front** - and §2 below now records that
it was mapped and **does not exist**: 18 configurations, 2 outside tolerance, both warned, none
silent. The dispatch fix this document originally argued for lost its mandate along with the
solar case; it is kept in §3 for the record, not as a recommendation.

**What this document listed as left — all four are DONE:**

1. ~~`UnmarkedDiscontinuityWarning`'s advice is unestablished for single points.~~ **Done, PR
   #28.** The warning now separates the two cases and carries the numbers for each.
2. ~~Nothing covers the *routing* the two shipped bugs exposed.~~ **Done, PR #29** — and writing
   the test found a live defect: the hybrid path stood aside for the cumulative scan even when
   the caller had disabled it with `cumulative=False`, costing 256x at the seam
   (1.157e-05 -> 2.966e-03). Fixed, with the `cumulative=False` semantics clarified:
   `strategy='magnus'` is the exact route to pre-1.0.0 numbers at every N, and
   `cumulative=False` guarantees only that the cumulative scan is not used.
3. ~~`implementation_details.rst` says the seam is 25.~~ **Done, PR #28** — four stale claims
   corrected; the one surviving `N = 25` is marked as historical.
4. ~~Example notebooks.~~ **Done, PR #30** — `13_magnus_tabulated_solar_model.ipynb` and
   `14_magnus_supernova_shock.ipynb`, a matched pair where the contrast is the lesson.

**What is genuinely left**, none of it urgent:

* **Nothing detects broadband sub-grid roughness** (§13.7, §13.14). The cheap statistic that
  would see it is described there; it needs its own false-positive measurement before it could
  ship, and the turbulence errors it would catch are already caught by the convergence machinery.
* **The single-point dispatch question** (§13.15), much less urgent since §13.17 withdrew the
  solar case. The surgical form is: have the hybrid path yield *to the cumulative scan
  specifically* rather than decline into the chain.
* **The detector chain works by accident, not design** (§2 below): on the rows that are outside
  tolerance, the 6400-point probe grid says *resolved* and the 200-point grid is what catches
  them. Fragile to any change in the probe ladder.
* **GitHub Actions cannot run** — every job reports a billing block, so nothing in CI has been
  verified since before PR #28. All three PRs were verified locally only (758 tests, `ruff`, docs
  under `-W`). The long-standing "GitHub Pages is disabled" item is **indistinguishable from the
  billing block** until Actions run again.

---

## 2. THE SILENT BAND - MAPPED, AND IT DOES NOT EXIST

**Answered before this handover was needed.** `shock_silent_band.py`, 9 widths (2100 km down to
7 km) x 2 energies, on the **averaged** observable: 18 configurations, **2 outside tolerance,
both warned, 0 silent**. The w = 1e-3 spike is real - 9.840e-04, a 6x excursion above its
neighbours, reaching 98 % of tolerance while quiet - but it never crosses. The package is honest
across the whole sweep.

Two things that came out of it and are still live:

* **The detector chain works by accident, not design.** On both outside-tolerance rows the
  **6400-point** probe grid says *resolved*; what catches them is the **200-point** grid saying
  *unresolved*. The coarse probe rescues the fine one. The fine grid's verdict genuinely is
  uninformative about whether the transport can integrate the front - so the masking is real but
  fragile, and anything that changes the probe ladder could unmask it.
* **`t_breakpoints` is not the cure it is documented to be at single points.** On the averaged
  observable across those 18 configurations it improved 7, worsened 11, and pushed 2 answers that
  were *inside* tolerance to *outside* it (w = 5e-4 at 15 MeV: 6.048e-04 -> 1.275e-03;
  w = 1e-4 at 45 MeV: 2.858e-04 -> 1.682e-03). Same mechanism as §13.6 - declaring breakpoints on
  a single point moves dispatch to the general Magnus ladder, which is often worse there.
  **`UnmarkedDiscontinuityWarning` tells the user to pass `t_breakpoints`, and for a single-point
  call that advice is unestablished.** That is the smallest well-defined piece of work left.

Full data: `FINDINGS_ROBUSTNESS_PROGRAMME.md` §13.20, rows in `shock_band_rows.npy`.

---

## 2c. The original open question, for the record

| shock width | what it represents | error (averaged) | resolution test | warns |
|---|---|---|---|---|
| 1e-6 (0.07 km) | a real hydrodynamic front | 0.2135 | **144/144** | yes |
| **1e-3 (70 km)** | **a front smeared over simulation cells** | **9.773e-04** | **0/144** | **no** |

A real shock is mean-free-path thin, so w = 1e-6 is the physical object and it is loud. But
nobody hands this package a real shock - they hand it a **simulation snapshot**, where the front
is smeared across a few grid cells, i.e. tens of km. **The silent band is plausibly the most
likely form a user's shock actually arrives in**, and 9.773e-04 against a 1e-3 tolerance is two
per cent of margin: luck, not headroom. Shift the width, the energy or the tolerance and it goes
outside, silently.

**Hypothesised mechanism** (two pieces of evidence, not a measurement): at 70 km the front spans
about six cells of the 6400-point **probe** grid, so ``_profile_is_resolved`` calls it resolved,
while the far coarser **transport** grid still straddles it. If so, the detector is checking a
grid the answer does not depend on.

**The job:** sweep width x energy x tolerance across w = 1e-2 ... 1e-4, find where "unflagged but
outside tolerance" lives and how far outside it gets, on the **averaged** observable
(``avg_check2.py`` is the harness). If the band is real, the fix is to key the resolution test to
the transport grid rather than the probe grid - much smaller and better targeted than anything
else in this document.

---

## 2b. The two exposures as originally written (solar now WITHDRAWN)

Both measured against `solve_ivp`/DOP853 at `rtol = 1e-12`, oracle verified to 5e-10 at every
energy (`physical_battery.py oracle_check`).

| | A — real solar model | B — supernova shock |
|---|---|---|
| profile | BS2005-AGS,OP table, log-interpolated | Fogli et al. envelope, 70 km fronts |
| configuration | d = 2, **5 MeV**, single point | d = 3, 15 MeV, single point |
| error | **1.380e-03** (linear: 1.178e-03) | **1.095e-03** |
| `strategy_info['certified']` | `True` | `True` |
| warnings | **none** | **none** |
| what the cumulative scan gives | **1.707e-05** | **2.586e-06** |

**A is WITHDRAWN.** It read as the most ordinary input in the population — a real published solar
model at an ⁸B energy in the standard two-flavour treatment — but the error is phase, not
envelope: averaged over six oscillation lengths it is **2.603e-05**, inside tolerance by 38x.
Solar physics measures the averaged survival probability, so no user sees it. B's *sharp* end
survives averaging (0.2135) and is the real result; B's 70 km end is the open question in §2.

Reproduce either with `bs05_energy_band.py`, or directly:

```python
import physical_profiles as pp, harness as H, magnus.oscprob as op
f = next(x for x in pp.families() if x['label'] == 'BS05(AGS,OP) cubic')
info = {}
P = op.osc_prob_matter_std_potential(2, f['ne'], 5.0e6, f['l1'], H.params_for(2),
        L0=f['l0'], density_is_of_number_of_electrons=True, strategy_info=info)
# info['engine'] == 'hybrid', info['certified'] is True, and P is 1.38e-3 from the truth
```

---

## 3. The fix that works, and the reason it is not shipped

**Both exposures disappear** if `HYBRID_YIELDS_TO_CUMULATIVE_MIN_POINTS` goes 8 → 1 **and**
`CUMULATIVE_AUTO_MIN_POINTS` goes 2 → 1. Verified: A 1.380e-03 → 1.707e-05, B 1.095e-03 →
2.586e-06. Both constants are needed — 8 is a threshold on *point count*, and a single point
cannot cross any such threshold however far it is lowered.

**And it is not a cost trade.** `physical_battery.py seam_cost`, three alternating rounds with a
control that returned 1.00x / 1.01x / 1.06x, reading minima rather than medians:

| N = 1 | cumulative / hybrid cost |
|---|---|
| BS05 solar model | **0.02x** (681 ms → 15 ms) |
| turbulence | **0.01x** |
| Earth crust | 0.15x |
| supernova shock | 0.27x |
| tabulated | 0.62x |

Cheaper on every profile, median 0.15x. One row in the whole 25-row sweep costs anything
(`tabulated` at N = 2, 1.61x).

**So why is it reverted?** Three reasons, in order of weight.

1. **Standing aside does not mean the cumulative scan answers.** The dispatch order is
   hybrid → interaction picture → separable → cumulative. When the hybrid path declines, the next
   engine that *applies* takes the request, and on several workloads that is `ip_exp`. Two tests
   caught it as `assert 'ip_exp' == 'hybrid'`. **The cost measurement above compares hybrid
   against cumulative. It does not measure hybrid against whatever comes next, which is what
   actually happens.** This is the sharpest objection and it invalidates the simple change, not
   the goal.
2. **The cumulative branch has never carried default traffic.** Making it the default surfaced
   two latent defects within minutes — both fixed in this session, both **pre-existing and
   reachable on the shipped tree** through an explicit `cumulative=True` (§4). Two found that
   quickly is evidence about the branch's test coverage, not about those two bugs.
3. **It is a change of default, not a dominant engine.** The cumulative scan's worst error over
   the 76 physical workloads it serves is 5.10e-03, and on one of them the hybrid path is 15x
   better.

**The surgical fix**: have `_osc_prob_hybrid_dispatch` yield **to the cumulative scan
specifically** — decline only when the cumulative scan will actually take the request — rather
than decline into the dispatch chain and hope. That preserves the accuracy gain without handing
requests to `ip_exp` or `separable` as a side effect. `_cumulative_scan_would_serve` already
exists and is the natural place to build on.

---

## 4. Two shipped bugs found on the way, both fixed here

Both were unreachable by default only because `CUMULATIVE_AUTO_MIN_POINTS = 2`; both are
reachable today with an explicit `cumulative=True`, and no test covered the combination.

* **Missing scalar squeeze.** `osc_prob_3nu_earth(E, ..., cumulative=True)` with scalar energy
  and baseline returned `(1, d, d)` instead of `(d, d)`, so `P[nu_i][nu_f]` silently selected a
  row of the wrong array. The cumulative branch of `osc_prob_energy_baseline` was the only return
  site in that function that did not apply `__getitem__(0 if return_float else slice(None))`.
* **`convergence_info` forwarded to an engine that rejects it**, raising
  `TypeError: magnus_expansion_multislab() got an unexpected keyword argument` instead of
  returning a probability. Dropped now, for the same reason `strict_convergence` already was:
  the traversal walks a fixed grid and has no refinement ladder to report on.

**Neither has a regression test yet.** That is the smallest useful piece of work in this document
and should probably be done first.

---

## 5. What is still open

* **The two exposures** (§2), unfixed.
* **Broadband sub-grid roughness is invisible to both structural detectors.** A statistic that
  would see it is nearly free — `find_hidden_features` already computes total variation on the
  fine and reference grids and discards their ratio, which is 1.000–1.003 on all 17 families with
  nothing to hide against 1.59 and 3.19 on the two turbulent ones. Not built, because the
  turbulence errors are already caught by the convergence machinery, so a new detector needs its
  own false-positive measurement first (§5c discipline).
* **`docs/source/implementation_details.rst` says the seam is 25** in three places (lines 110,
  116, 255). Stale since the previous tranche — the constant has been 8 throughout this session.
* **Example notebooks for A and B**, requested by the user: one each, showing the failure, the
  diagnosis and the fix. Only writable once there is a fix to show.
* **GitHub Pages is still disabled**, so the docs workflow fails on every commit. Nothing in the
  codebase can fix it.

---

## 6. Traps that cost real time in this session

* **A `SIGSTOP` corrupts a `perf_counter` stopwatch.** Pausing the timing battery for 75 minutes
  made the in-flight call report the whole pause as its own execution time. Kill and re-run a
  timing measurement rather than resuming it; accuracy measurements resume fine.
* **`min`, not median, is the statistic for timing on a shared machine.** Interference only ever
  adds time, so the fastest observation is closest to the uncontended cost.
* **A monotone step is not a "hidden feature".** `find_hidden_features` computes
  `per_interval - endpoints`, and the two reference nodes bracketing a step see its full height,
  so the statistic is zero at *any* width. Measured: fine/coarse variation ratio 1.000 from
  w = 1e-2 to 1e-6. Labelling the shock family as hiding sub-grid structure was a construction
  error, caught only by re-deriving the ground truth by measurement.
* **The resolution test flags the shock widths that are already loud and misses the silent one**
  (0/144 at w = 1e-3; 144/144 at w = 1e-6). Any cure gated on that flag cannot fire on the case
  that needs it.
* **Ask "does *any* row satisfy this", not "does the *worst* row satisfy it".**
  `bs05_energy_band.py`'s first verdict reported the largest silent miss and then tested that one
  row for band membership, concluding "does not reach the solar band" when two of the three do.
* **`pgrep -f foo` matches the shell issuing it.** Kill by PID from a targeted `ps`.
* **Never leave two waiters polling the same file.** One `until … sleep … done` per condition.
