> **Archived 2026-09-29.**  Part A merged as PR #89.  Part B became issue #120, merged as PR #128
> and closed.  The folder that section 8 describes was split: #120's harness went to
> `tools/auto_tight/`, the measurement records to `docs/dev/measurements/`, and the rest (logs,
> process files, drafts, the sessions' transcripts) was archived outside the repository.  Paths
> below are as they were.

# Handover: `strategy='auto'` takes the Magnus ladder at tight tolerances (plus the patch-fix wrap-up)

Written 2026-09-24, about 23:50 local, by session `93bacad5` for a fresh session. **Read this whole
file first, then act.** Everything here was verified unless it is marked OPEN. Do not re-derive it.
The author is very low on credits. Keep turns few, run long jobs in the background (a background
run costs nothing), read outputs with `tail`/`grep`, and launch no agents.

---

## 0. Reloading cheaply

- **This folder**, `resources/handover_auto_tight/`, is untracked on purpose. Never commit it, and
  delete it once the work is merged. It sits outside every lint and test path (CI's ruff checks
  `src/magnus/ tests/ docs/ notebooks/make_notebooks.py`; the file-tree test only sees tracked
  files). §8 lists its contents.
- **Dialogue of the day that led here:** `session/tail_convo.txt` (132 KB, from 2026-09-24 16:34
  UTC to the move to the cloud; `[#####]` = line in the full transcript). Everything before that is
  in `session/last_compact.txt`. The full transcript is
  `~/.claude/projects/-home-mbustamante-Research-magnus/a0d2e3ee-05d1-4f5e-8936-3df18ac27060.jsonl`,
  403 MB. Never read it whole; pull single lines with `sed -n 'Np'` plus `json.loads`.
- The auto-loaded memory (`MEMORY.md`) holds the standing rules. The ones that bite here are
  restated in §1.

## 1. The author's standing rules (non-negotiable)

1. **No agents or workflows** unless the author answers yes to a question that names how many and
   the rough cost. Ultracode being on is for reasoning depth only; it is **not** consent.
2. **Never rewrite cached data or recorded numbers or timings** without their decision:
   `notebooks/paper_figure_cache.json`, `notebooks/external_*.json`, the timing outputs of
   notebooks 13/24/25/28, and measurements quoted in docstrings. Report them instead.
3. **Commit, push, open a PR or file an issue only when asked.**
4. **Be brief and lead with the answer.** Verify every claim before stating it. Put both sides of
   any comparison on the same basis (fresh process vs warm, same engine or engine named).
5. **Fix the default path; never work around it in an example.** The counterexample: I forced
   `strategy='magnus', magnus_exp_order=8` into Listing 1 instead of fixing `auto`, and the author
   rejected it ("Why are you forcing strategy 'magnus' … Stop wasting my time").
6. **Do not edit `resources/paper/main.tex`.** Give tex for the author to paste.
7. **Writing style:** avoid the ", and" clause-join tic and the trailing ", which". A series of
   nouns is fine.
8. **Run the gates (§6) before saying anything is done.**
9. On long sessions, restate the plan and your position in it every turn.

## 2. The goal (the author's words, 2026-09-24)

> "you said that, after the latest edits, the only way for Listing 1 to run fast was by passing
> 'magnus' at order 8 explicitly. But we had been working so that this behavior was automatically
> reached by passing 'auto' instead (and explicitly asking for order 8). … This work is
> unfinished. Finish it."

Confirmed end state: **Listing 1 passes no `strategy` (so `'auto'`) and adds
`magnus_exp_order=8`. `auto` hands it to the Magnus ladder at order 8 by itself.** `strategy_info`
shows the hybrid declining with `'auto prefers the ladder'`, and a test locks this in.

**The order stays the user's dial** (default 4). Letting `auto` choose the order was rejected: the
paper (≈ line 614) says the package defaults to `'gl'` at order 4, and that raising the order "pays
only once the tolerance is tight enough". Without the order keyword, `auto` may still pick the
ladder at order 4 (≈ 5 s on Listing 1), but only if the measurement says so.

It all started with one sentence for the paper (§7): how long Listing 1 takes as printed and at the
default 1e-3. At 1e-3 `auto` took 9 s, which became #70 (merged). At 1e-12, `auto`'s hybrid was
2e-11 off while certified, which became the patch fix (Part A).

## 3. State at handover (verified about 23:40 local)

- **Branch `fix-hybrid-patch-tolerance` at `baabe81`.** The working tree is clean apart from this
  untracked folder.
  - `105703d`, made by the local session: the patch fix plus its tests.
  - `baabe81`, made by the cloud session after the move and fast-forwarded locally at 23:09 from
    `origin/fix-hybrid-patch-tolerance-4s8qlf`: a CHANGELOG "Fixed" entry and one sentence in
    `docs/source/adiabatic_strategy.rst`. Its numbers match `hp/b_*.jsonl`: the 0.3 R☉ chord is
    9.7× outside at rtol 1e-9 on the hybrid's own criterion, and 2.96e-10 off and warned after the
    fix.
  - `origin/fix-hybrid-patch-tolerance` is still at `105703d`, so local is one commit ahead. Delete
    the cloud's branch `…-4s8qlf` only on the author's go.
- **`main` = `origin/main` = `d9c4a3d`:** PR #72 merged, #70 closed, CI green.
- **`stash@{0}`** = "author's main.tex: Listing 1 edit (strategy='magnus', magnus_exp_order=8), to
  decide later". This is my rejected workaround. Dropping it is the author's call. The listing in
  `HEAD`'s `main.tex` is the original (no strategy, no order).
- **Branch `fix-auto-tolerance`** is merged but not deleted, locally or on origin.
- **Issues.** Open: #52 and #53 (hold), #71 (the separable engine can land outside a loose
  tolerance). Not yet filed, waiting for the author's go:
  - (a) the patch-tolerance defect;
  - (b) `docs/source/index.rst` says `auto` handles the full Sun, which holds only at loose
    tolerances. At 1e-7/1e-9 the hybrid declines and the ladder answers with
    `ToleranceNotAchievedWarning`, which is warned, not silent.
- **The cloud session** (`claude.ai/code/session_01Tt2xUDcnHiNLTmVavvXSRr`) can't be read from
  this machine. The local session tools and the browser pane can't reach it (not signed in). Only
  `baabe81` came back. OPEN: whether it ran 4c, `turb_convergence` or the gates. Assume it didn't.
- **Environment** (identical locally; pin it in any cloud run): Python 3.12.7, numpy 1.26.4,
  scipy 1.15.3, numba 0.60.0, llvmlite 0.43.0, matplotlib 3.10.8, mpmath 1.3.0, pytest 9.0.3,
  pytest-xdist 3.8.0, nbclient 0.10.4, nbformat 5.10.4, sphinx 9.1.0, ruff 0.15.1. Linux x86_64,
  glibc 2.39, 12 cores, ≈ 6 GB RAM free. **CI runs NumPy 2.x.**

## 4. Background (verified facts)

### 4.1 The patch-tolerance fix (done, on the branch)

- **Defect.** `adiabatic._hybrid_propagator_once` called `_local_evolution_operator(...)` without
  the requested tolerance, so every Magnus patch converged to its default `patch_atol = 1e-7`. The
  hybrid certifies by comparing successive levels, and two levels that hold the same patch agree
  trivially. It therefore returned `certified=True` below about `atol + rtol = 1e-6` while wrong.
  Examples:
  - Listing 1 at 1e-12 (3ν, 2 MeV): 2.0e-11 off DOP853.
  - A 2ν solar chord (0.3 R☉, 10 MeV): 4.8e-9 off at every tolerance.
  - Proven by intervention: changing only `patch_atol` took the error from 2.0e-11 (1e-7) down to
    3.0e-14 (1e-13).
- **Fix.** `hybrid_propagator` passes `patch_atol = min(1e-7, (atol + rtol)/10)` down through
  `_hybrid_propagator_once(..., patch_atol=1e-7)`, a keyword whose default is the old value. There
  is no floor, on purpose: a floor would re-create the bug. A patch that can't converge within its
  cap leaves the result uncertified (`HybridCertificationWarning`), and `auto` then falls back.
  Tests: the two fakes in `tests/test_adiabatic.py` (≈ lines 356 and 380) take `patch_atol`, and
  two regression tests were added that fail on `main` and pass on the fix.
- **Verified:**
  - 4a: bit-identical at `atol + rtol ≥ 1e-6`, 29 hybrid-path results (`hp/verify_a.py`,
    `a_head.npz` vs `a_fix.npz`).
  - 4b: 140 cases (137 with DOP853 references) at rtol 1e-7/1e-9/1e-12, `atol = rtol/100`. Silent
    misses: **78 on `main`, 0 on the fix**, under both `'hybrid'` and `'auto'` (`hp/verify_b.py`,
    `b_head.jsonl`, `b_fix.jsonl`).
- **Reach.** Only `oscprob`'s hybrid scan (≈ `oscprob.py:6211`) reaches the hybrid below 1e-6.
  `avgprob`'s two calls (≈ `avgprob.py:634`, `_crossing_from_windows` via `level_crossing_matrix`)
  run at the default 1e-3, so they are unchanged. The full Sun can't regress: the hybrid was already
  uncertified there at 1e-7/1e-9, since its patches can't converge within the cap.

### 4.2 The #70 route (merged in PR #72): where Part B changes code

All in `src/magnus/oscprob.py`. Line numbers are approximate; grep for the names.

- **Constants (≈ 708–775)**, all `versionadded:: 1.1.1` (unreleased, so free to change) and all
  exported in `__all__` (≈ 23375). A constant missing from `__all__` breaks the strict docs build
  (ten broken references last time).
  - `AUTO_LADDER_MAX_PHASE = 1e4`: the phase is ∫(λmax − λmin) dl, from 17 probes at up to 5
    energies.
  - `AUTO_LADDER_MAX_FLOOR_FRACTION = 0.25`: the ladder's starting slab count must be at most a
    quarter of its cap. The Gauss-Legendre cap is 20 000, and every full solar path starts at
    8 300–21 000 slabs, which keeps the Sun on the hybrid.
  - `AUTO_LADDER_MIN_TOLERANCE = 1e-6`: **the cut-off to lift or replace.** Its docstring says "the
    hybrid strategy reaches about 1e-9 … at no extra cost", measured with the broken patches, so it
    is now false. Rewrite it only with the author's OK.
  - `AUTO_LADDER_TOLERANCE_MARGIN = 10`: the ladder runs at tol/10 (#71: 5.8e-3 at a requested
    1e-3).
- **`_PreferLadder` (≈ 5791)** carries `min_n_slabs`. Its `request()` (≈ 5815) divides rtol and atol
  by the margin and raises `min_n_slabs` to the convergent floor, capped at the resolved
  `max_n_slabs`. The ladder refines ×1.5 per rung.
- **`_estimated_phase` (≈ 5835)** returns `(phase, n_floor)`.
- **`_auto_prefers_ladder` (≈ 5877)** has **the line to change**:
  `tols = [t for t in (rtol, atol) if t > 0.0]; if not tols or min(tols) < AUTO_LADDER_MIN_TOLERANCE: return None`.
  It then checks phase and floor and runs the resolution test (`UnmarkedDiscontinuityWarning`).
  Finally it calls `_note_engine('hybrid', answered=False, reason='auto prefers the ladder',
  estimated_phase=..., min_n_slabs=..., tolerance_margin=...)`.
  **Note:** the tolerance tested is `min(rtol, atol)`, which for Listing 1 is atol = **1e-14**,
  not 1e-12.
- **Callers.** `_osc_prob_hybrid_dispatch` and `_osc_prob_hybrid_dispatch_generic` return the
  marker under `strategy == 'auto'`. The four call sites apply `request()`, skip the
  interaction-picture engine (`ip_exp`, 2ν only), then try separable and finally the general
  ladder:
  - `osc_prob_matter_std_potential`
  - `osc_prob_matter_nsi`
  - `osc_prob_liv`
  - `_osc_prob_with_potential`

  `magnus_exp_order` reaches the ladder untouched (checked in the old session, `[110663]`).
  `integration_method` is a local variable in both dispatchers (`[109260]`).
- **Docs that describe the route:**
  - the `strategy` docstring (≈ 7908);
  - `docs/source/engines.rst` (dispatch table and ladder-route section);
  - `docs/source/adiabatic_strategy.rst`;
  - `CHANGELOG.md` `[Unreleased]` (the #70 "Changed" entry says "1e-6 or looser");
  - `docs/dev/DECISION_DISPATCH_ORDER.md` (#70 addendum).
- **Tests:**
  - `tests/test_auto_ladder.py`: 14 tests, some asserting that a tight tolerance stays on the
    hybrid. Replace those.
  - `tests/test_engines.py`: a fixture keeps `auto` on the hybrid for six engine tests.
  - The old runtime inventory logged at least 32 hybrid calls below 1e-6 in the suite
    (`hp/hp_tests.log`). Those may now reroute.

### 4.3 Listing 1 numbers so far (hybrid before the patch fix; warm; four curves × 140 energies)

| at rtol 1e-12, atol 1e-14 | time | worst error vs DOP853 at 1e-13 |
|---|---|---|
| `auto` → hybrid, order 4 (the listing as it was) | 7.7 s | 2.0e-11, certified (now fixed) |
| hybrid at order 8 | 9.3 s | – |
| ladder, order 4 | 4.8–5.0 s | at DOP853's level |
| ladder, order 6 | 0.84 s | at DOP853's level |
| ladder, order 8 | **0.43 s** | ≤ 9.6e-14 at all 26 reference energies per curve (median 1.2e-14) |

- **At 1e-3,** `auto` takes the route: 0.040 s warm (3–22 ms per curve), 0.17 s in a fresh
  process. The 0.12 s one-off start-up lands on whichever call comes first (checked by reversing
  the curve order). The import (0.24 s) sits outside both. The hybrid at 1e-3 took 8.0 s, with an
  answer bit-identical to its 1e-12 one.
- **Fresh-process, as printed** (hybrid era): median 8.2 s (7.96–8.46), per curve
  1.4/1.9/2.0/2.9 s.
- **Why the ladder must win on Listing 1.** Every Listing 1 point has one non-adiabatic window
  covering the whole path. The hybrid there costs its window search (≈ 87% of a call at loose
  tolerance) plus at least two full-path patches, so the ladder wins at any tolerance. In general
  the ladder can lose only where windows cover under half the path and its full-path cost exceeds
  the window search: solar chords and long adiabatic baselines.
- **Fig. 1 needs no rerun.** Its data are nb28's cached `oracle` entry, computed with
  `oscprob.osc_prob` (ladder, order 4, 1e-12), which none of this touches. The order-8 listing
  reproduces the plotted curves to ≤ 3.7e-12, at the level of the figure's own residuals
  (2.7e-13/2.0e-12/1.3e-12/2.7e-12). An order-8 bottom panel would sit under the oracle floor.

### 4.4 Smoke run of `tight_measure.py` (3 cases, rtol 1e-12 only, 3 workers; NOT a result)

- All three were eligible (the route is taken at 1e-3): L1-3nu 2 MeV, sun2-exp-0.1R 5 MeV, and the
  L1-3nu scan. **Partial solar chords are eligible**, so any rule change reaches them. Measure them
  with care.
- Ladder time as a fraction of the hybrid's at the same order (median / max):
  - L4/10: 0.28 / 0.89
  - L4/1: 0.09 / 0.23
  - L8/10: 0.16 / 0.19
  - L8/1: 0.16 / 0.16
- No silent miss in any configuration.
- Honest warnings:
  - L4/10 warned on all three. OPEN which warning; probably `ToleranceNotAchievedWarning`, since
    ÷10 at 1e-12 asks the ladder for 1e-13/1e-15, near round-off.
  - L4/1 warned on one.
  - L8 warned on none.
  - H4 warned on one (the solar case, presumably uncertified, then fell back).
- Listing 1's 3ν scan (26 energies) at 1e-12:
  - H4: 0.66 s (hybrid)
  - H8: 0.40 s (hybrid)
  - L4/10: 0.18 s (separable)
  - L4/1: 0.054 s (separable)
  - L8/1: **0.008 s** (separable)

  All within their limits.

---

## 5. THE PLAN (in order; **STOP** marks where to report and wait)

### Part A: finish the patch fix (branch `fix-hybrid-patch-tolerance`)

- **A1 (old step 4c): hybrid cost before and after.** Report only.
  - "After" is the H4/H8 columns of B1.
  - "Before" means the same calls on `main`:

    ```bash
    git archive d9c4a3d | tar -x -C $S/headtree
    PYTHONPATH=$S/headtree/src:$S/headtree/tests python ...
    ```

    Export both `src` and `tests`, because `conftest.py` defeats a bare `PYTHONPATH` (memory
    `magnus-conftest-defeats-pythonpath`).
  - Minimum: H4 and H8 on Listing 1's four scans at 1e-7/1e-9/1e-12 on both trees. Run
    `tight_measure.py` against the export with a subset of jobs, or adapt `i70/warm_both.py`.
- **A2 (old step 7): `turb_convergence`.** This is nb28's cached entry at rtol 1e-8 that can reach
  the hybrid. Recompute it in scratch only, main against the fix, and report the difference.
  **Never write `notebooks/paper_figure_cache.json`.** Find its compute function in
  `notebooks/make_notebooks.py` (`grep -n turb_convergence`) and call that function directly from a
  scratch script. Don't go through `cached()`, which computes and stores on a miss.
- **A3: report, don't edit,** the recorded numbers that are now suspect:
  - the `hybrid_propagator` docstring in `src/magnus/adiabatic.py`: ≈ 1744, "converges in 7
    iterations" at 1e-12, and ≈ 1753, a γ table at 1e-9;
  - the `AUTO_LADDER_MIN_TOLERANCE` docstring: "about 1e-9 at no extra cost".
- **A4: the paper's `prob_vs` passage,** which cites rtol 1e-8, "0.8 s vs 15.5 s" and "better than
  $10^{-10}$". It ran on the hybrid below 1e-6. Its declining example's phase is ≈ 8e4 rad, so it
  stays on the hybrid even after Part B. Check whether the patch fix moves its numbers. Report
  only.
- **A5: gates (§6). STOP.** Ask the author:
  - to push (this fast-forwards `origin/fix-hybrid-patch-tolerance` to `baabe81`);
  - to open a PR;
  - to file issues (a) and (b).

  Suggest landing Part A as its own PR, since it is a correctness fix. Part B would then go on a new
  branch (say `auto-tight-tolerance`) from this one, or from `main` after the merge. The author
  decides.

### Part B: `auto` takes the ladder at tight tolerances

- **B0: the pre-registered question.** On requests that pass #70's phase and floor conditions, at
  tolerances tighter than 1e-6: is the ladder, run as the route runs it at the user's order, inside
  the tolerance and faster than the fixed hybrid at the same order? And with which margin?
- **B1: measurement** (background; nothing written to the repo):

  ```bash
  S=<this session's scratchpad>
  python resources/handover_auto_tight/tight_measure.py run $S/tight.jsonl 6 > $S/tight.log 2>&1   # run_in_background
  python resources/handover_auto_tight/tight_measure.py score $S/tight.jsonl
  ```

  - **Jobs:** 144 in all, the 140 step-4b points (§8) plus the 4 Listing 1 scans over their 26
    reference energies.
  - **Configurations:**
    - H4 and H8: `auto` as on the branch, which means the hybrid;
    - L4/10, L4/1, L8/10 and L8/1: the route forced open by setting
      `AUTO_LADDER_MIN_TOLERANCE = 0` and the margin to 10 or 1. Only on eligible cases.
  - **Tolerances:** rtol 1e-7/1e-9/1e-12, with atol = rtol/100.
  - **Run time:** expect 30–60 min on 6 workers. Step 4b's 840 calls took about 5 min.
  - **Waiting:** wait on the background job's completion notification. A `pgrep` loop matches
    itself and never ends (memory `magnus-wait-loop-self-match`).

  **Pass conditions, fixed now, before the data:**
  - **Accuracy:** zero silent misses. A silent miss is an error above
    `max(atol + rtol·P_ref, DOP853 spread)` with none of `HybridCertificationWarning`,
    `ToleranceNotAchievedWarning` or `UnmarkedDiscontinuityWarning` raised.
  - **Speed:** on every case the rule sends to the ladder, ladder time ≤ 1.2 × the hybrid's at the
    same order.
  - **Warnings:** a ladder configuration that warns `ToleranceNotAchievedWarning` where the hybrid
    certified is a regression for that case. Count those; among passing configurations, take the
    one that warns least.
  - The scorer prints all of this per tolerance and configuration, with the worst cases listed.
    Re-time borderline ratios (1.0–1.5) alone on an idle machine before trusting them.
- **B2: choose the rule,** the simplest that passes, in this order of preference:
  1. **Lift the cut-off entirely.** The route applies at any tolerance under the existing phase and
     floor conditions, plus a margin rule. The smoke run suggests ÷10 hurts at 1e-12, so for
     example the margin applies only while `min(rtol, atol) ≥ 1e-6` (the #70 regime, unchanged bit
     for bit) and is 1 below that. The argument for dropping the margin: in the asymptotic regime,
     the difference between rungs overestimates the error by about 1.5^p − 1 (≈ 4 at p = 4,
     ≈ 25 at p = 8).
  2. **If some eligible workloads lose by more than 1.2× at tight tolerance,** add a phase limit
     scaled by the tolerance: `phase ≤ AUTO_LADDER_MAX_PHASE × (tol/1e-6)^(1/p)`, with p the
     user's `magnus_exp_order`. The ladder's slab count grows as tol^(−1/p) at order p, while the
     hybrid's window search stays fixed. Check that it separates the measured cases.
  3. **If neither separates them,** compare an estimate of the ladder's needed slab count with its
     cap.

  **Scope:** Gauss-Legendre only (the default `integration_method`). Trapezoid, Simpson and
  cumulative keep today's behaviour. The order is never chosen by `auto`.

  One more fact belongs in the rationale. Where the fixed hybrid is uncertified at a tight
  tolerance (partial solar chords), `auto` already ends on the ladder after paying for a discarded
  hybrid attempt. Going straight there is strictly cheaper unless the ladder then fails.

  **STOP.** Show the author the score table and the rule with its numbers, briefly, and get an OK.
  The approach itself is already approved; this stop is for the numbers. Docstrings that record
  the new measurements need that OK too.
- **B3: implement** (`oscprob.py`, dispatch layer only):
  - `_auto_prefers_ladder`: replace the `min(tols) < AUTO_LADDER_MIN_TOLERANCE` return with the
    chosen rule. Pass in what it needs: `magnus_exp_order` and `integration_method` from both
    dispatchers.
  - `_PreferLadder.request()`: the margin rule.
  - Constants: repurpose or rename `AUTO_LADDER_MIN_TOLERANCE`, for example as the tolerance below
    which the margin is 1. Update the docstrings with the B1 numbers (after the OK) and keep
    `__all__` in step.
  - Record the order and margin used in the `strategy_info` detail. Update the `strategy`
    docstring (≈ 7908).
- **B4: tests** (in `tests/test_auto_ladder.py`):
  - **Listing 1 as it will be printed** (no `strategy`, `magnus_exp_order=8`, rtol 1e-12,
    atol 1e-14; a few energies per curve keeps it fast):
    - the ladder answers (`'separable'`);
    - the hybrid declined with `'auto prefers the ladder'`;
    - error ≤ 1e-12 against DOP853 at the worst point (3ν, 2 MeV; DOP853 at 1e-13 takes 0.15 s
      there).
  - The same call without the order keyword: whatever the rule decides, within tolerance.
  - `strategy='hybrid'` still runs the hybrid; `strategy='magnus'` is unchanged.
  - A non-Gauss-Legendre `integration_method` keeps today's behaviour.
  - The full Sun stays on the hybrid at 1e-9, e.g. 2ν at 10 MeV over 0.9 R☉.
  - Replace the old "tight tolerance stays on the hybrid" tests, and check that the
    `test_engines.py` fixture still pins the hybrid at tight tolerances.
  - Run the full suite and justify every changed test on its own.
- **B5: docs.** `engines.rst` (table row and ladder-route section), `adiabatic_strategy.rst`,
  `CHANGELOG.md` `[Unreleased]` (amend the #70 entry; don't duplicate it), the
  `docs/dev/DECISION_DISPATCH_ORDER.md` addendum, and `index.rst` if it states `auto`'s tolerance
  behaviour. The author's style: brief.
- **B6: notebook route audit** (scratch, report only):
  - **Mirror:** `rsync -a --exclude .git ./ $S/mirror/` (≈ 560 MB; check free RAM first).
  - **Run:** execute the notebooks with a kernel hook that logs every call taking the route with
    `min(rtol, atol) < 1e-6`. Reference implementations:
    - `i70/audit_w.py`: a worker taking an explicit notebook list;
    - `i70/ipy/`: the IPython profile whose startup loads the hook;
    - `i70/audit.py`: the hook and runner;
    - `hp/hp_hook.py`, `hp/hp_plugin.py`, `hp/run_nb.py`: the hybrid-call variant, also usable as
      a pytest plugin.

    Their paths point into `/tmp/...a0d2e3ee.../scratchpad`, so fix those. **Load-test the hook in
    one kernel first:** a relative path once made a whole audit silently log nothing. Use 5
    workers under `setsid`, report each rerouted cell once, don't flood with monitors, and wait on
    completion.
  - **Known candidates:**
    - nb24: 4 hybrid calls at rtol 1e-7…1e-10, all uncertified, so another engine answered.
      **nb24 is timed; never rebuild it without the author's decision.**
    - nb28: cache-only, so it can't move except through A2.
  - **Report** which notebooks change, and rebuild nothing without an OK. If a rebuild is approved,
    restore the committed timing outputs (as was done for nb13 cells 27/29/31/35). If only
    markdown changes, use the prose-only patch (memory `magnus-timed-runs-stay-untouched`).
- **B7: gates (§6). STOP.** Report, then ask about commit, push and PR.
- **B8: the Listing 1 numbers for the paper** (§7). Time the listing exactly as it will be
  printed, on an idle machine, with the engine read from `strategy_info` every time:
  - **Two bases:** fresh processes (median of 5, the protocol of `i62/time_lst_validation.py`)
    **and** warm (like `i70/run_listing_warm.py`).
  - **Two tolerances:** 1e-12/1e-14 and the default 1e-3. At the default, run it both with and
    without `magnus_exp_order=8`, and settle with the author what "at the default" means.
  - `i70/run_listing.py` and `run_listing_warm.py` time **my rejected `strategy='magnus'`
    listing**. Edit their `common` dict before reuse.

  Give the author the filled sentence and state the basis.

## 6. Gates (all of them before "done")

- **`git add` new files first.** `test_file_tree` compares against `git ls-files`, so an untracked
  file is invisible to it.
- **`pytest tests -n auto`** (≈ 8.5–10 min, in the background). Read the summary line: the exit
  code of a piped `tail` is always 0.
- **Strict docs, into a fresh directory** (an incremental build can't re-find a cached warning):

  ```bash
  sphinx-build -n -W --keep-going -b html docs/source $S/docs_strict_N
  ```

- **`python tools/lint_notebook_cells.py`**
- **`ruff check src/magnus/ tests/ docs/ notebooks/make_notebooks.py`**, exactly CI's command.
- **If the generator changed:** `pytest $(grep -l make_notebooks tests/*.py)`.
- **Re-run only what a change can move** (memory `magnus-rerun-only-what-can-move`).
- **After a push:** `gh run list --commit <FULL sha>` (the short sha shows nothing), with
  `GIT_CONFIG_NOSYSTEM=1` and the issue body on stdin (memory `magnus-gh-in-sandbox`). CI runs
  NumPy 2.

## 7. The paper (the author edits `main.tex`; hand over tex)

- **Listing 1, `common` dict:** add only `magnus_exp_order=8` after `rtol=1e-12, atol=1e-14`, with
  no `strategy`. The line must stay within 56 characters, the listing width.
- **The timing sentence** for the main-text paragraph describing Listing 1. Slots are filled from
  B8.

  ```latex
  Run as printed, Listing~\ref{lst:validation} takes about [T_TIGHT]~s on a laptop, [T_CURVE_MIN] to [T_CURVE_MAX]~s per curve of 140 energies. It asks for ${\tt rtol} = 10^{-12}$ and ${\tt atol} = 10^{-14}$, far beyond what any oscillation analysis requires, only to match the comparison in \figu{validation}; at the default of $10^{-3}$, the four curves take about [T_DEFAULT]~s.
  ```

  Both numbers go on one basis. If the engines differ between the two tolerances, name that. An
  earlier suggestion for the 1e-3 clause: "\magnus\ hands the curves to its slab ladder instead
  (Sec.~\ref{sec:engines})". With order 8 in the listing, a pointer to where the order is
  introduced also helps (§ "Choosing a quadrature", `sec:quadrature`).
- **Caption parenthetical.** The author's draft is "(The defaults are $10^{-3}$; here they are
  cranked up for illustration.)". Suggested: "(The defaults are $10^{-3}$; here they are tightened
  to match the comparison in \figu{validation}.)", or drop it, since the sentence covers it.
- **≈ line 266:** "The deviation … is $10^{-12}$ or smaller" becomes "a few times $10^{-12}$ or
  smaller". The plotted residual reaches 2.7e-12.
- **≈ line 1034** (engines section, after "…otherwise the request passes to the next engine."):
  the #70 clause proposed earlier needs rewriting after Part B, to cover tight tolerances and the
  margin. The earlier version:

  > {\tt 'auto'} also hands the Magnus ladder the requests it answers faster---loose tolerances over
  > moderate accumulated phases, such as Listing~\ref{lst:validation} at the default
  > tolerance---and runs it at a tenth of the requested tolerance.

  The sentence after it ("The accuracy of the default is therefore the accuracy of the certified
  answers") then needs "and of the ladder's".
- **The `prob_vs` passage,** if A4 moves its numbers.

## 8. What is in this folder

- **`tight_measure.py`:** the B1 runner and scorer. Smoke-tested on 2026-09-24 (§4.4). It doesn't
  depend on where it's run from, and it imports `hp/verify_b.py`.
- **`hp/`**, the patch-fix work:
  - `verify_a.py`, `a_head.npz`, `a_fix.npz`: 4a, bit-identity.
  - `verify_b.py`: the 140-case list (`cases()`), `call_for()` and the DOP853 reference code.
  - **`b_fix.jsonl`: the reference set.** Per case: `ref12`, `ref13` and `phase`, plus results
    on the fix. **Never overwrite it.**
  - `b_head.jsonl`: results on `main`.
  - `hp_hook.py`, `hp_plugin.py`, `run_nb.py`, `ipy/`: the runtime inventory hook.
  - `hp_notebooks.log`, `hp_tests.log`, `tests_run.log`, `w*.out`: the inventory's results.
- **`i70/`**, the #70 work (paths inside point into the old `/tmp` scratch):
  - Campaign and route checks:
    - `campaign.py`, `campaign.jsonl`: 20 workloads with DOP853 references at 1e-12, 87 records;
    - `score.py`, `crossover.py`;
    - `validate.py`, `validate_sun.py`, `validate*.jsonl`: route validation;
    - `measures.py`, `measures.jsonl`: old and new phase measures;
    - `routes.py`: which engine answers under the rule.
  - Notebook route audit: `audit.py`, `audit_w.py`, `ipy/`, `route_audit.log` (the 4 notebooks
    rerouted under #70: 13, 20, 22, 27).
  - Listing 1 checks:
    - `check_l1.py`, `rungs.py`, `l1_ref.npz`;
    - `run_listing.py` (fresh process), `run_listing_warm.py`, `warm_both.py`;
    - `engines_both.py`: every engine at both tolerances;
    - `order_tight.py`: orders 4/6/8 at 1e-12;
    - `worst_dop.py`;
    - `ladder_check.py`, `ladder_check_*.json`: the order-8 listing against DOP853 at 26
      energies;
    - `fig1_vs_listing.py`: the listing against the cached Fig. 1 data.
  - Patch fix and solar: `patch_probe.py` (the root-cause probe: vary only `patch_atol`),
    `sun_patch.py`, `sun_auto.py` (solar behaviour before and after the fix), `sun10.py`,
    `sun_floor.py`.
  - Other: `nb24_shock.py`; `patch_*.py` and `blk_*.txt` (the #70 edits, already applied;
    history only); `suite*.log`, `docs_strict*.log`.
- **`i62/`:** `time_lst_validation.py`, `time_lst_default.py`, `time_lst_default_reversed.py` and
  their logs. This is the original fresh-process protocol for timing Listing 1 (median of five).
- **`session/`:** `tail_convo.txt` and `last_compact.txt` (see §0).
- **Not copied, because they can be regenerated:**
  - the `main` export: `git archive d9c4a3d | tar -x -C DIR`;
  - the repo mirror: `rsync`;
  - the NumPy-2 venv (`python -m venv`, then numpy 2.x), needed only for cache-key portability
    checks.

## 9. Traps that cost time before

- A notebook's code sits inside string literals in `make_notebooks.py`, so ruff never sees it.
  `tools/lint_notebook_cells.py` covers that gap.
- CI executes the committed `.ipynb` on NumPy 2, not the generator on numpy 1.26. Verify the
  artifact CI runs.
- `git checkout -- notebooks/` also reverts the generator. Restore only `'notebooks/*.ipynb'`.
- `make_notebooks.py --only 2` matches twelve notebooks. `--no-execute` destroys outputs.
- `nohup`/`setsid` record the launcher's PID, not Python's. Find the real one before waiting on it.
- A monitor that prints every call floods the conversation. Deduplicate per cell.
- Scratch execution renders figures at a different size than `make_notebooks.py`. Never splice
  scratch outputs into committed notebooks. Rebuild with the builder.
- "Agreement between two refinements" proves nothing if both hit the same cap. A flat plateau means
  a fixed-width slab, not convergence.
- Two numbers compared must use one basis: warm vs fresh, same engine, same tolerance meaning
  (`min(rtol, atol)`).
