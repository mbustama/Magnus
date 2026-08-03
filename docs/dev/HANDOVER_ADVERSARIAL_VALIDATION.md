# Handover: adversarial validation of the dispatch, cumulative and adiabatic changes

**Written:** 2026-08-03, at the close of the session that produced PR #23 and the eleven commits
on `notebook-breakpoints-and-cumulative`.
**Starting point:** branch `notebook-breakpoints-and-cumulative`, eleven commits ahead of `main`
(which already contains PR #23). **Nothing is pushed.**
**Machine for every number quoted:** 12-core box, min-of-2-or-3 timings.
**Read first:** `DECISION_DISPATCH_ORDER.md`, `DECISION_CUMULATIVE_DEFAULT.md`,
`NOTES_ADAPTIVE_REFINEMENT.md` §4b–§4c. This brief assumes all three.

---

## The job

Six audit passes have been run on this branch. They found and fixed nine defects. The last two
passes found no numerical errors, only a wrong test expectation of mine and one change worth
reverting. **That convergence is not proof of correctness, and this brief exists because it
should not be treated as such.**

Build and run **at least six batteries** of hard tests: a full regression over cases already
seen, cases that are hard and have never been run, and cases *designed to break this code*. Be
aggressive. A battery that passes everything on the first attempt is more likely to be a weak
battery than a correct codebase — six passes running at roughly two defects per pass says the
prior on "nothing left" is low.

**Report what breaks. Do not fix and move on quietly** — a defect found here is worth more as a
documented, reproduced failure than as a silent patch.

---

## What changed, and therefore what is under test

Eleven commits, three distinct pieces of work.

| # | change | where |
|---|---|---|
| 1 | `strict_convergence=True` — two consecutive agreements, opt-in, default off | `osc_prob` |
| 2 | Notebook breakpoints: castle wall gains `l_ini`/`l_fin`; four PREM cells gain `t_breakpoints` | notebooks 02/03 |
| 3 | `cumulative='auto'` becomes the default | `osc_prob_energy_baseline` |
| 4 | Hybrid stands aside for single-energy baseline scans of ≥ 25 points under `strategy='auto'` | `_osc_prob_hybrid_dispatch` |
| 5 | The cumulative grid's probe is **always strict**; its `MagnusConvergenceWarning` suppressed | `osc_prob_energy_baseline` |
| 6 | `CUMULATIVE_N_ACC_SAFETY` 2 → 4 | module constant |
| 7 | `strategy='magnus'` opts out of the cumulative scan | three wrappers |
| 8 | γ is swept **along the probe grid**, not only at gap extrema; one window per contiguous run | `find_nonadiabatic_windows` |

Current constants, which several tests below key off:

```
CUMULATIVE_AUTO_MIN_POINTS             = 2
HYBRID_YIELDS_TO_CUMULATIVE_MIN_POINTS = 25
CUMULATIVE_N_ACC_SAFETY                = 4
cumulative default                     = 'auto'
strict_convergence default             = False
max_n_slabs: 20000 ('gl'), 2000 ('trapezoid'/'simpson')
```

Baseline to beat: **672 passed**, ruff clean, docs clean under `-W`.

---

## The oracles, and which is which

This is the single most important section. The project has been burned three times by comparing
one of its own paths against another.

- **`solve_ivp`/DOP853 is the only accuracy oracle.** Use `rtol=1e-12, atol=1e-14`; verify
  convergence by tightening to `rtol=1e-13` and confirming the answer moves by ≪ the error you
  are quoting. Below ~1 MeV over a solar radius it costs ~2 min and ~9M evaluations — affordable,
  but budget for it.
- **`expm` for constant `H`** tests time-ordering and indexing, not accuracy. The cumulative scan
  currently matches it to 2.31e-14 over 60 baselines.
- **Identical-grid comparison with `rtol=atol=None`** is the *bookkeeping* oracle: it tests
  ordering and indexing only.
- **Unitarity** is necessary, never sufficient. Every wrong answer this session was exactly
  unitary.
- **Never** score one Magnus path against another. "Cumulative agrees with per-point" tells you
  nothing about which is right.

---

## Traps, every one of which cost real time in this session

Instrument bugs, not code bugs. Each one produced a confident wrong measurement.

**Physics / API**

- `matter.vcc_func_from_rho_func(rho, L0, ratio, e_fraction, nubar, in_g_per_cm3, is_number_density)`
  — the **7th positional argument** must be `True` for `gd.NUM_DENSITY_E_SUN_CENTRAL`. Passing
  `nubar` there yields a potential ~10⁹ too small, i.e. a silently *vacuum* reference that looks
  perfectly converged. This invalidated a whole first measurement round.
- Splatting `gd.OSC_PARAMS_PREDEFINED[...]` forwards `name`/`description`, which the probability
  functions reject. Filter to the six mixing parameters.
- **`H` closures must take exactly one parameter.** `lambda l, E=E: ...` is rejected with
  "the provided H_func is a function of more than one parameter", or worse, mis-dispatched. Use a
  factory returning a genuine one-argument closure.
- `t_slab_edges` takes **pairs** `[[t0,t1],[t1,t2],…]`, not a flat edge array — and supplying it
  **silently discards `t_breakpoints`**. A flat array produced `|ΔP| ≈ 0.5` and looked like a
  physics result.
- `solve_ivp` without `t_eval` stores `y` at every accepted step; harmless at high energy,
  ruinous at low.

**Measurement discipline**

- **This machine's run-to-run variance is enormous.** `test_adiabatic.py` moved 124 s → 69 s
  between two runs of *identical* code. Never compare a `main` run against a branch run taken at
  a different time — **alternate**, several rounds, and treat any effect smaller than the scatter
  on untouched work as noise. One whole-suite comparison (392.9 s vs 331 s) was pure artefact.
- Run nothing else while timing. Short concurrent commands contaminated two early rows.
- `pgrep -f <script>` matches the waiting shell's own command line, so
  `until ! pgrep -f foo; do …` never exits.
- Foreground commands are killed at 10 minutes. Use `run_in_background` plus an
  `until <condition>; do sleep N; done` waiter, and check repo state afterwards if a chained
  `git stash`/`checkout` run was interrupted mid-flight.
- Always run measurements under `(ulimit -v 10000000; …)`.
- `gh` is snap-confined: it cannot read `/etc/gitconfig` — prefix `GIT_CONFIG_SYSTEM=/dev/null`
  — and cannot read a `--body-file` under `/tmp`; pipe through `--body-file -`.

**Tooling**

- `nbclient` dispatches through `async_execute_cell`; overriding the sync `execute_cell` records
  **nothing** and the run still reports success. Cost a full 14-minute notebook pass.
- When attributing an `osc_prob` call to a scan axis by AST, use the **nearest** enclosing loop.
  Taking every enclosing loop mixes verdicts on the PREM cells (a `for i in range(len(costhz))`
  around a comprehension) and shifted the cost split by 3 percentage points with 12% unclassified.

---

## Where I expect it to break, ranked

Written down in advance so the battery can be aimed, and so that "we found nothing" can be
weighed against what was actually probed. **These are my honest suspicions, not known failures.**

1. **The detector's fixed probe grid.** `find_resonance_candidates` /
   `find_nonadiabatic_windows` sample `n_probe = 200` linear points (doubling to 6400 under
   refinement). Structure narrower than `(l1-l0)/200` can be stepped straight over. My new γ
   sweep inherits exactly this grid, so it has the same blind spot — it fixed the case where γ is
   *broad* and missed at the extrema, and does nothing for a genuinely **narrow** feature. This
   is the same class as the `suggest_n_slabs` blind spot in `NOTES_ADAPTIVE_REFINEMENT.md` §2.
   **This is where I would look first.**
2. **The 2 ≤ N < 25 band.** Scans in that range still go to hybrid. On the multi-resonance
   profile, N = 8 measured 7.48e-02 both before and after the routing change — the routing
   rescues only N ≥ 25. Commit 8 should have improved this band too; confirm it did, on profiles
   other than the one it was developed against.
3. **`n_acc` from a capped probe.** Over a full solar radius at 5 and 10 MeV the strict probe
   reaches `max_n_slabs` without converging, so `n_acc` is ceiling-derived. Lowering
   `max_n_slabs` therefore lowers a whole scan's resolution silently, in proportion. Untested
   below the default.
4. **Cost of the γ sweep on hostile profiles.** One window per contiguous run, grown from its
   peak — a profile where γ exceeds threshold over most of the path gives one enormous window
   patched by exact Magnus. Could be extremely slow. Untested.
5. **Unitarity loosened** from 8.9e-16 (hybrid, exact by construction) to ~8e-12 (cumulative, a
   product over tens of thousands of slabs). Fine against the suite's 1e-9, but untested at
   N ≫ 10⁴ where the product is longer.
6. **Very large N.** Memory measured flat at ~56–59 MB up to N = 4000. Untested at 10⁵–10⁶.
7. **Quadrature methods.** `trapezoid`/`simpson` cannot reach 1e-3 on a full solar radius by any
   path (a single adaptive call gives 1.00e-01, exhausting both caps). Everything about them is
   "better than the alternative", never "good".

---

## The batteries

Six are required; seven are sketched. Sub-tests are the unit of work. Each battery should state
its pass criterion **before** running, and report every configuration, not only failures.

### Battery 1 — Full regression: every past case, against the new code

Re-run everything this branch was ever measured on, and check nothing drifted.

1. **The bit-identity set.** These must be *bit-identical* to `main`; anything else is a
   regression. Dump raw probabilities on both and diff element-wise:
   `strategy='magnus'` scan; a sub-threshold scan (N < 25); a single point; a vacuum scan; a
   constant-density scan; an energy scan at fixed baseline; `osc_prob_earth` PREM;
   `average=True`; an explicit `cumulative=False` scan. All nine were bit-identical when last
   checked — with the caveat that commit 8 (the γ sweep) landed *after* that check and touches
   the hybrid path, so **the single-point and sub-threshold rows may legitimately have moved
   now**. Establish which, and justify each.
2. **The 48-configuration accuracy sweep**: 1–100 MeV × N ∈ {8, 40, 150} × L ∈ {0.4, 1.0} R_sun,
   scored against `solve_ivp`. Last result: 32 better, 16 equal, 0 worse; worst 2.93e-05 against
   main's 3.21e-05.
3. **The 50-configuration dispatch grid**: standard/NSI/LIV × 0.5–100 MeV × {0.9, 1.0} R_sun.
   Confirm hybrid still certifies where it did and that nothing newly fails to certify.
4. **Notebooks 02 and 03**: re-run with `nbclient`, hash every embedded PNG, and check every
   changed figure against `solve_ivp` rather than assuming a finer grid is better — error is
   **not** monotone in `n_slabs`.
5. **Cross-entry-point consistency**: `osc_prob` vs `osc_prob_energy_baseline` vs the wrappers on
   the same physics, within tolerance.
6. **Every item above at d = 2, 3, 4 and 5.** The sweeps behind these commits were 2ν and 3ν
   only; re-running them at 4ν and 5ν is regression testing, not new work. See Battery 6.

### Battery 2 — Aimed at the detector's fixed probe grid (most likely to break)

The premise: `n_probe = 200` linear samples cannot see structure narrower than the spacing.

1. **A single narrow resonance.** Place one crossing of width `w` at a random position, sweep
   `w / (l1-l0)` from 1e-1 down to 1e-5. Find the width at which the detector stops finding it,
   and check whether `hybrid_propagator` still reports `certified=True` there. **A certified
   answer that misses a resonance is the headline result if it occurs.**
2. **Aliasing.** A sinusoidally modulated profile whose period is close to `(l1-l0)/200`,
   `(l1-l0)/100`, and rational multiples thereof — the probe should alias. Scan the period
   through the resonant values.
3. **Edge crossings.** A resonance within one probe spacing of `l0`, and of `l1`. Also one
   exactly at a probe point, and one exactly midway between two.
4. **Many crossings.** 50, 100, 200 crossings across the range — more than the probe has points.
5. **Clustered crossings.** All crossings packed into 1% of the range, the rest quiet.
6. For each: does refinement (`n_probe` doubling to 6400) rescue it, and does `certified` tell
   the truth? Score against `solve_ivp` throughout.
7. Repeat the narrow-resonance sweep at d = 4 and 5. More level pairs means more chances for
   one pair's window to mask another's missed crossing, which cannot happen at d = 2.

### Battery 3 — The routing seams

Every threshold is a discontinuity, and discontinuities are where behaviour hides.

1. **N = 1, 2, 3, 24, 25, 26, 100** on each of: solar exponential, multi-resonance, castle wall
   with breakpoints, PREM, noisy. Accuracy should not jump discontinuously across N = 25.
   **A large accuracy discontinuity at the seam is a defect even if both sides are inside
   tolerance.**
2. The 2 ≤ N < 25 band on multi-resonance profiles specifically — see suspicion 2 above.
3. Escape hatches at every N: `strategy='magnus'`, `strategy='hybrid'`, `cumulative=False`,
   `cumulative=True` (must raise when inapplicable, never fall back silently).
4. Position-independent Hamiltonians must never reach the cumulative scan; confirm at every N
   and every d — the exclusion tests `isinstance(H_first, Callable)`, which is d-agnostic, so
   a failure here would be a surprise worth having.
5. `t_slab_edges` given (cumulative must decline); `t_breakpoints` given (cumulative must
   accept and honour them — verify they are actually in the grid, not silently dropped).

### Battery 4 — Extreme numerics

1. **Energy**: 0.05, 0.1, 0.5 MeV and 1, 10, 100 GeV. At the low end confirm the `solve_ivp`
   oracle is itself converged before trusting any error.
2. **N**: 10⁴, 10⁵, 10⁶ — wall time, and peak memory via `tracemalloc` (expect O(output), and
   confirm the chunking really is flat in N). Watch for the `_check_output_fits` refusal.
3. **Tolerances**: 1e-1, 1e-3, 1e-6, 1e-9, 1e-12. Confirm each is either met or
   `ToleranceNotAchievedWarning` is raised — the suppression added in commit 5 must not have
   swallowed it. Last checked: 1e-9 and 1e-12 both warn correctly.
4. **`max_n_slabs` lowered** to 500, 2000, 5000 on a scan whose probe would otherwise cap. Does
   the answer degrade silently? Suspicion 3.
5. **Geometry**: `L0` ≠ 0 (including `L0` mid-profile); baselines spanning eight orders of
   magnitude; unsorted, duplicated and degenerate baseline arrays; `L == L0` exactly.
6. **Flavors**: see Battery 6, which makes this a first-class axis rather than a sub-test.
7. Antineutrinos throughout, not as an afterthought.

### Battery 5 — Designed to break, not to confirm

Adversarial construction. The goal is a reproducible wrong answer.

1. **Fabricate a profile that defeats the γ sweep.** γ just under threshold over a long stretch,
   so nothing opens but the accumulated non-adiabaticity is large. This is precisely the
   mechanism that made hybrid certify at 4.3e-02, and my fix raises the bar without removing the
   mechanism — a *sub-threshold* accumulation is still invisible.
2. **Fabricate a profile where the cumulative probe stops early.** The probe is strict, but
   `NOTES_ADAPTIVE_REFINEMENT.md` §4b's "frozen grid" mode (agreement between bit-identical
   answers) survives strictness. A user-supplied `t_breakpoints` set that dominates the grid
   should reproduce it — and now it would misplace a whole scan, not one point.
3. **Discontinuous and pathological profiles**: step functions with unmarked edges; a profile
   with a kink but no jump; `C⁰` but not `C¹`; a profile with a genuine singularity approached
   but not reached.
4. **Random-profile fuzzing.** Generate a few hundred random smooth profiles (random Fourier
   sums with controlled bandwidth), run each through every entry point, and score against
   `solve_ivp`. Report the empirical distribution of error, and every case outside its requested
   tolerance **that does not warn**. This is the highest-yield sub-test in the battery.
   Draw the flavour count randomly too, from {2, 3, 4, 5}, with random sterile mixings — a
   fuzzer restricted to 3ν would repeat the blind spot this brief is trying to close.
5. **Adversarial `n_slabs`/`n_tpts_per_slab`** supplied by the caller, deliberately absurd
   (1, 2, 10⁶).
6. **Concurrent/repeat determinism**: same call twice bit-identical; shuffled input order gives
   identical answers (last measured 0.00e+00); `n_jobs > 1` agrees with `n_jobs = 1`.

### Battery 6 — Flavor count as a first-class axis (required)

**Every measurement behind these eleven commits was 2ν or 3ν.** 4ν and 5ν were touched exactly
once, in Pass 5, and only structurally — shape, finiteness and unitarity (1.8e-11 and 1.1e-11).
**Neither was ever scored against `solve_ivp` on the new path.** There is therefore no accuracy
evidence whatsoever for the sterile cases, and this battery exists to supply it.

Why it is not merely "run the same thing with a bigger matrix":

- The detector loops over **every level pair**: d(d−1)/2 = 1, 3, 6, 10 for d = 2…5. The γ sweep
  added in 9c7945a runs per pair and appends windows from each, which are then merged across
  pairs. Over-merging across *positions* was already a real bug at 3ν, caught by
  `test_crossings_too_close_in_phase_are_reported`; over-merging across *pairs* is untested.
- Sterile mixings put resonances at densities the three-flavour case has none at, so 3+1 and 3+2
  are genuinely different physics, not bigger arithmetic.
- `ip_exp` is gated to 2ν precisely because a neglected term's coefficient jumps three orders
  from 2 to 3 flavours (see the gate comment in `oscprob.py`). Whether anything analogous grows
  with d on the adiabatic or cumulative paths has never been asked.

Sub-tests:

1. **Accuracy against `solve_ivp`, 4ν and 5ν, on the new path** — single points and scans either
   side of N = 25, standard matter, over a solar profile. This is the gap; do it first.
2. **4ν/5ν with NSI**, and **with LIV** (`n_liv` = 0 and 1). The wrapper families reach
   `_osc_prob_hybrid_dispatch` at every flavour count, so all three dispatch sites need it.
3. **Multi-resonance at 4ν and 5ν.** Repeat Battery 2's constructions with sterile mixings large
   enough to produce extra crossings; count windows per pair and check the merged set against a
   dense per-pair γ scan. Confirm no crossing is swallowed by a window belonging to a different
   pair.
4. **Degenerate and near-degenerate levels.** Set two mass splittings equal, or nearly so, so
   that a gap is ~0 over a stretch. `_point_adiabaticity` returns `inf` when the gap is exactly
   zero — check what the sweep then does, and that `find_resonance_candidates`' bisection on
   `f_jk` behaves.
5. **Cost scaling.** One `eigh` per probe point plus a per-pair sweep: measure how the γ sweep's
   1.31× at 2ν grows at d = 4 and 5, where there are 10× the pairs.
6. **Antineutrinos at 4ν and 5ν**, which compounds two axes neither of which is covered.
7. Confirm `ip_exp` still declines for d > 2 at every entry point — it is the reason 3ν/4ν/5ν
   were unaffected by PR #23, and that must remain true.

Pass criterion: same as the rest — met tolerance or a warning, and `certified=True` implying
accuracy. Flag separately any case where accuracy *degrades monotonically with d*, which would
suggest a d-dependent term is being neglected somewhere.

### Battery 7 (optional but recommended) — Cross-module and oracle diversity

1. `expm` for constant `H` across dimensions 2–5.
2. The analytic two-flavour vacuum formula, where it applies.
3. Composition: `U(0→L₂)` against `U(L₁→L₂) · U(0→L₁)` computed independently.
4. `avgprob` against a brute-force window average of the oscillating probability — the averaged
   solar result matched a 25-oscillation window average to 1.5e-04 when last checked.
5. The adiabatic windows against a dense, independent γ scan: every stretch with γ > threshold
   should lie inside a reported window. **This directly tests commit 8's premise.**

---

## Pass criteria

State them before running. Suggested:

- **Regression**: no configuration less accurate than `main` by more than 5%, and none newly
  outside its requested tolerance.
- **Correctness**: every answer either meets the requested tolerance or raises
  `ToleranceNotAchievedWarning`. A silent miss is a failure regardless of magnitude.
- **Certification honesty**: `certified=True` from `hybrid_propagator` must imply the answer is
  within the requested tolerance of `solve_ivp`. This is the property that was false before
  commit 8; prove it is true now, or find where it still is not.
- **Performance**: no configuration more than 10% slower than `main`, measured by alternating.
- **Unitarity**: better than 1e-9 everywhere.

---

## Deliberately not done, and why

- **Splitting the branch.** Commit 9c7945a (the adiabatic fix) is a different module from the
  other ten and arguably wants its own PR. Left as one branch so this validation sees the whole
  change at once; split afterwards if the batteries pass.
- **Pushing.** Nothing is pushed. PR #23 is already merged to `main`; these eleven commits are
  local only.
- **The 2 ≤ N < 25 band.** Known exposed on multi-resonance profiles. Not addressed because the
  right fix is in `magnus.adiabatic`, and commit 8 may already have fixed it — Battery 3 should
  determine that rather than assume it.
- **`strict_convergence` as a default.** Measured at 1.53× median; declined because the errors it
  fixes are largely phase errors invisible to an averaged solar observable
  (`NOTES_ADAPTIVE_REFINEMENT.md` §4c). Note it is now **load-bearing internal machinery**
  regardless: the cumulative probe always uses it, so weakening it would silently degrade every
  baseline scan.
- **GitHub Pages** is still disabled on the repository, so the "Documentation Deployment"
  workflow has failed on every commit for months. Settings → Pages → source "GitHub Actions".
  Nothing in the codebase can fix it.
