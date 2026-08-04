# Handover: the robustness programme (five items) and the implementation-details page

**Written:** 2026-08-04, at the close of the session that produced PRs #24 and #25.

**Where to work: branch `dev-robustness` from `main`, and do everything in it.**

**Base is `main`.** PRs #24, #25 and #26 are all merged; `main` contains both tranches as of
2026-08-04 (`d2e83ed`). The stack did briefly go wrong -- #25 was based on
`cumulative-and-notebooks`, which `main` had already taken via #24, so #25's merge never reached
`main` until #26 -- and the guard below is kept because that failure was silent:

```bash
git fetch origin
git show origin/main:src/magnus/adiabatic.py | grep -c '^GAMMA_TO_ERROR'   # must be 1
```

If that prints 0 you are on a `main` without the adiabatic tranche -- the resolution test, the
measured constants, the bounded patch budget, `UnmarkedDiscontinuityWarning` and 11 regression
tests -- and every item below assumes it is present. **Do not start on a `main` that fails that
check**: Item 1's acceptance test compares the current code against `main`'s older

behaviour, so it needs both.
**Read first:** `FINDINGS_ADVERSARIAL_VALIDATION.md`, then this. The batteries and harness it
describes are committed under `docs/dev/adversarial_batteries/` and you should reuse them
rather than rebuild them.

**Baseline to hold:** 685 passed on `main`, `ruff check src/ tests/` clean, docs build under
`-W`. (`ruff check .` reports 63 findings, all in `notebooks/`, identical before this work --
pre-existing.) The suite takes ~11 minutes; see expectation 5 below for why, and for the one
test whose tolerance is arguably over-specified.

**Definition of done for the tranche**, beyond each item's own criterion:

- the baseline above still holds, with the test count grown rather than shrunk (a *falling*
  count means something stopped being collected -- that happened once here and was nearly
  missed);
- every new constant or threshold carries the population it was measured on, in its docstring,
  as `GAMMA_TO_ERROR` and `RESOLUTION_RATIO` now do;
- every new warning has had its false-positive rate measured before shipping;
- anything asserted in a docstring or a decision document has been *run*. Three claims were
  found false this way -- a comment claiming a reuse that did not happen, a documented cure that
  had never been tested, and a warning telling users to ignore it when it mattered.

**A practical note on the batteries.** `docs/dev/adversarial_batteries/README.md` explains how to
run each one and how to drive `bitident.py`/`attribute.py` across two worktrees. Their *outputs*
(`.txt`, `.npy`) are deliberately gitignored, so the raw distributions behind the numbers in
`FINDINGS_ADVERSARIAL_VALIDATION.md` are not in the repo -- re-run to regenerate them. Budget
for that: the smooth fuzzer is ~25 minutes at 150 cases, dominated by the `solve_ivp` oracle.

**Branches.** `main` is the base. `notebook-breakpoints-and-cumulative` is historical -- it holds
this brief's own history and is merged nowhere; do not build on it.

---

## The job

Six deliverables. Five are robustness work agreed after the adversarial validation; the sixth is
a documentation page that ties them together. Work them on `dev-robustness`; they touch the same
files as each other, so a single branch is right.

The five exist because of a specific finding, not because they are good practice in general.
Across this session I added four safeguards and **three of them were wrong on first attempt**,
each caught by a *different* instrument and none by the one I would have guessed. The
programme below is aimed at that: the codebase's weak spot is not its physics, it is
**heuristics with hand-set calibration constants, certified by methods that check themselves.**

---

## Running this alongside other work

**Yes -- keep the machine loaded.** None of items 1, 2, 4, 5 or 6 needs a quiet box: they are
accuracy and correctness work, and every number they produce is deterministic. Three caveats,
all learned the hard way here.

**Ratios survive load; absolute times do not.** The alternating harness
(`adversarial_batteries/timing.py`) carries two control workloads the change cannot touch. If the
controls come back at ~1.00x, the ratios are trustworthy however busy the box is -- that is not a
theory, it is how the reuse-fix timing was taken, with a fuzzer pegging a core throughout, and
the controls came back at exactly 1.00x. **Never quote an absolute millisecond figure from a
loaded run.**

**Where a quiet machine is still needed.** Only two places. Item 2 for any constant whose
trade-off is cost against accuracy -- `threshold0` is exactly that, since lowering it opens more
windows and each window is patched -- and Item 3, where the point is to size a job to a CI
budget. Both want the alternating harness anyway.

**Watch memory, not CPU.** The box has ~15 GB and the batteries run under `ulimit -v 10000000`
(10 GB each). Two heavy jobs contend for memory long before they contend for the 12 cores, and
the failure mode is a kill rather than a slowdown.

**Make stalls distinguishable from slowness.** Under load, "no output for ten minutes" is
ambiguous, and it was misread twice in the previous session -- once concluding a fuzzer had died
when it was running at 100% CPU, once the reverse. Give every long job a progress line with
`flush=True`, and check CPU time (`ps -o time=`) rather than elapsed before concluding anything.

**Do not run two pytest suites at once.** Both write, both look alike in `ps`, and killing "the"
one kills the wrong one.

---

## What changed already, so you are not surprised by it

| constant | was | now | why |
|---|---|---|---|
| `adiabatic.GAMMA_TO_ERROR` | 1.0 (+ `GAMMA_SLACK` 2.0) | **0.85**, slack removed | measured over 149 configurations |
| `adiabatic.RESOLUTION_RATIO` | (new) 0.75 | **0.70** | 192 smooth vs 30 piecewise configurations |
| `_local_evolution_operator` `max_n_slabs` | 500 000 | **32 768** | legitimate patches ≤ 12 800; a 102 400-slab one should decline |

New in `magnus.adiabatic`: `_profile_is_resolved`, `_H_on_grid`, `find_nonadiabatic_windows(...,
info=)`. New in `magnus.oscprob`: `_n_required_params`, `_resolve_cumulative_kwarg`,
`UnmarkedDiscontinuityWarning`.

---

## Item 1 — Cross-method agreement instead of self-agreement

**The finding.** Every silently-wrong result found this session came from a method certifying by
comparing itself with itself. `hybrid_propagator` refines its own knobs and checks the answers
agree; `osc_prob`'s ladder does the same. When the method has a blind spot, both sides of the
comparison share it and the agreement carries no information. That is not a bug in either
comparison — it is a limit of the *shape* of the check.

**The idea.** The package contains genuinely independent engines: the hybrid strategy, the
general Magnus ladder, the cumulative scan, the two-flavour interaction-picture fast path, and
`scipy.linalg.expm` (exact for a constant `H`). Running two and comparing needs **no oracle at
all**, and detects exactly the class self-certification cannot.

**Build.** A public diagnostic in `magnus.oscprob` — working name `cross_check_strategies` — that
takes the same arguments as a wrapper call, runs whichever engines apply, and returns each
answer plus the pairwise spread. It must:

- never be on by default (it multiplies cost by the number of engines);
- report *which* engines actually ran, since most decline on any given request — reuse the spy
  pattern in `adversarial_batteries/battery3.py`;
- treat a large spread as a finding to report, not an exception to raise.

**Pass criterion.** On the profiles in `FINDINGS §3` where a method was silently wrong *before*
the fixes, the spread must be large. Check this by running against `main`'s `magnus` (a worktree
at `978663a`) — if the diagnostic cannot see the defects we already know about, it does not work.
That is the acceptance test, and it is stronger than any assertion on today's code.

**Then wire it into CI** over a small profile matrix, asserting the engines agree within
tolerance. This is the highest-value item; do it first.

---

## Item 2 — Audit the remaining constants

**The finding.** `GAMMA_TO_ERROR` was set from five configurations that all happened to sit in
one corner of the parameter space; it was wrong by up to 1.6×, then over-corrected when read
from the *unrestricted* population, and only came right when restricted to the regime the rule
actually governs. Of the four constants measured properly this session, **three were wrong.**

**The unaudited ones**, with the same provenance:

```
adiabatic:  threshold0=0.1   min_threshold=1e-6   n_probe0=200   max_n_probe=6400
            n_points0=201    max_n_points=12864   fd_step_frac=1e-6
            patch_atol=1e-7  n_slabs0=400         max_iters=12
oscprob:    CUMULATIVE_N_ACC_SAFETY=4   CUMULATIVE_AUTO_MIN_POINTS=2
            growth_factor_n_slabs=1.5   max_num_loops=50   min_n_tpts_per_slab=2
```

**Start with `threshold0 = 0.1`.** It decides whether a window opens at all, which makes it the
constant every other adiabatic safeguard is downstream of. `fd_step_frac = 1e-6` is second: it
sets the finite-difference step for every Hellmann-Feynman diagnostic, and nothing has ever
checked it against the step-size trade-off (truncation vs cancellation).

**Method, and it matters.** For each: identify *the regime in which the constant is consulted*,
sample that regime, and report the distribution — not one number. The `GAMMA_TO_ERROR` mistake
was reading a maximum over a population that included rows the constant never decides. Then
write the population into the docstring, as `GAMMA_TO_ERROR` and `RESOLUTION_RATIO` now do.

**Pass criterion.** Every constant either measured, or explicitly marked as unmeasured with the
reason. "It has always been 0.1" is not provenance.

---

## Item 3 — Fuzzers in CI, on aggregate assertions

**The finding.** The batteries found things the unit tests would not have, and nothing runs them.
They print tables and take tens of minutes.

**Build.** A seeded fuzz test — `tests/test_fuzz_statistics.py`, or a marked subset of the
existing suite — asserting on **statistics, not per-case values**:

- silent-miss rate ≤ some threshold (post-fix measurement was 4.1 % on smooth profiles);
- median error ≤ 1e-6;
- zero exceptions.

Per-case assertions on random input are brittle; aggregate ones are stable and still catch a
regression that moves the distribution. Size it to run in CI — 40–60 cases, not 250. Reuse
`battery5.sub4` and `battery8_piecewise.fuzz`.

**Trap.** The `solve_ivp` oracle dominates the cost, badly, at low energy and high flavour
count — measured at **25 s for 2 of 80 baselines** at 6.5 MeV, 3ν, over 1.4 scale heights, while
the package answered the whole case in 2.1 s. Restrict the fuzzer's energy range, or a CI job
will look hung when it is merely honest. If a fuzz run appears to stall, **suspect the instrument
before the package** — that was true every time this session.

---

## Item 4 — Property tests for cross-entry-point agreement

Oracle-free invariants, swept over a profile matrix. These encode "the answer must not depend on
which door you came in", and the dispatch seams are exactly where the defects lived.

- `strategy='auto'` vs `'magnus'` vs `'hybrid'`, within tolerance
- `cumulative=True` vs `False` vs `'auto'`, within tolerance
- a scan vs the same points computed one at a time
- shuffled baseline order → identical after unshuffling (currently exact, 0.00e+00)
- `n_jobs>1` vs `n_jobs=1` → identical
- composition: `U(0→L₂)` vs `U(L₁→L₂)·U(0→L₁)`
- unitarity via probability row/column sums (currently ≤ 1.6e-11 at N = 10⁵)

Several already exist as one-off checks in the batteries; the work is turning them into a
parametrised sweep. Overlaps Item 1 — build them together.

---

## Item 5 — Warnings: fix what exists, and surface certification

**5a. The warnings should say when and how to change something -- and one of them currently says
the opposite.**

`MagnusConvergenceWarning` ends with:

> *If a target tolerance (rtol/atol) was requested, the adaptive refinement narrows the slabs
> automatically and this warning can be ignored.*

That instruction is **false in exactly the cases where the warning matters**. Measured on a
sawtooth density with `rtol=atol=1e-3` explicitly requested:

```
strategy=auto     err=7.484e-03  (7x the requested 1e-3)   MagnusConvergenceWarning
strategy=magnus   err=7.484e-03  (7x the requested 1e-3)   MagnusConvergenceWarning
```

A tolerance *was* requested, the refinement *did* run, and the answer is seven times outside it --
while the message tells the user to ignore the warning. `DECISION_DISPATCH_ORDER.md` section 5
separately records it firing on rows accurate to 1.6e-06. So it is both over- and under-trusted:
it fires on good answers, and it disclaims itself on bad ones. **Fix that clause first; it is the
single most misleading string in the package.**

**The standard to hold every warning to.** Four things, in this order:

1. **What was detected** -- the condition, in the user's vocabulary.
2. **What it means for the answer** -- affected, unaffected, or unknown; and *by how much*, if
   the code knows. Several do not say this and could: `ToleranceNotAchievedWarning` already
   computes the disagreement between the last two refinement levels, so it can report how far
   from converged it stopped rather than only that it stopped.
3. **What to change** -- parameter name and direction, ideally magnitude.
4. **When it is genuinely safe to ignore** -- and only when that is true.

Graded against it, the package is already better than most:

| warning | verdict |
|---|---|
| `ScalarHamiltonianWarning` | **the model.** What, why (4.6x, measured), and the exact code change with a before/after snippet |
| `DensityUnitWarning` | excellent -- quantifies the consequence (~18 orders), names the tell (`P_ee = 1`), gives two ways out |
| `UnmarkedDiscontinuityWarning` | good -- says why refinement cannot help, and what to pass |
| `PhaseAveragingWarning` | good -- reports the s.e.m. and names the alternative function |
| `MagnusHighOrderCostWarning` | good -- states the trade rather than forbidding anything |
| `ToleranceNotAchievedWarning` | incomplete -- says which cap to raise, never how far off the answer is |
| `HybridCertificationWarning` | incomplete -- says the accuracy is uncertified, suggests **no action at all** |
| `MagnusConvergenceWarning` | **harmful** -- see above |

The first tranche of work is therefore: fix the one false clause, add the missing magnitudes
where the code already knows them, and give `HybridCertificationWarning` an action
(`strategy='auto'` for the automatic fallback, `t_breakpoints` if the profile has known
structure, or a looser tolerance).

**5c. New warnings are in scope -- both new instances of existing kinds, and entirely new ones.**

Several conditions are currently detected and then handled *silently*, which is defensible but
not obviously right. Candidates, offered as things to **measure**, not as a specification:

- **`UnmarkedDiscontinuityWarning` only fires on the cumulative path.** The same detector
  (`adiabatic._profile_is_resolved`) already runs inside `hybrid_propagator`, where an unresolved
  profile causes a silent decline to the general path. Firing there too is a new instance of an
  existing kind, and close to free.
- **The hybrid strategy declining under `strategy='auto'` is entirely silent by design.** That is
  the right default -- it happens on ordinary calls -- but a user debugging why a call got slow,
  or why a result moved, has no way to see it. An opt-in (`verbose`, or the `strategy_info` of
  5b) is probably better here than a warning.
- **A ceiling-derived `n_acc`.** When the cumulative scan's probe hits `max_n_slabs` without
  converging, the grid is sized from a cap rather than from a converged requirement, and the
  whole scan inherits it. `ToleranceNotAchievedWarning` fires, but says nothing about *the scan's
  resolution* being capped, which is the consequence that matters.
- **Marginal certification.** `hybrid_propagator` certifies when
  `GAMMA_TO_ERROR * gamma_max <= atol + rtol`. A result passing at 0.95 of that bound is
  certified on the same footing as one at 0.01, and the constant is itself only good to ~2x.
- **A profile the resolution test rejects at *both* densities** is a strong statement about the
  user's input, and currently produces only a quiet `certified=False`.

**The bar for any new warning**, which is the part not to skip:

1. It must meet the four-part standard above -- in particular it must say what to *change*.
2. Its **false-positive rate must be measured before it ships**, on the profile families the
   package actually serves. `UnmarkedDiscontinuityWarning` was measured at 0 false positives over
   192 smooth configurations and 12/12 true positives before being wired in; that is the
   precedent.
3. Run the full suite **and** `sphinx -W` afterwards. A new warning inside a `jupyter-execute`
   example is a docs build *error*, and the warning added last session crashed on a constant
   Hamiltonian -- caught by three existing tests, not by review.

Warnings caught **nothing** in the previous session, and the one added introduced a crash. That
is an argument for measuring them, not for avoiding them.

**5b. Expose the hybrid `certified` flag.** `hybrid_propagator` returns it, but under
`strategy='auto'` an uncertified result silently falls back and the caller cannot ask what
happened. Add a `strategy_info` out-parameter to the three scenario wrappers, following the
`convergence_info` convention already in `osc_prob`, reporting which engine answered and — for
hybrid — whether it certified.

**This was started and deliberately not finished.** The plumbing is more invasive than it looks:
the wrappers pass `**kwargs` into `scan_kwargs['kwargs']`, `_osc_prob_hybrid_dispatch` returns
`NotImplemented` if any *unrecognised* key remains there, and `osc_prob_energy_baseline` forwards
its `**kwargs` to `osc_prob`, which rejects unknown keywords. So a new keyword must be recognised
and popped at **every** layer it passes through, or it will silently disable the hybrid strategy
rather than instrument it. Nothing was committed; the tree is clean.

---

## Item 6 — `docs/source/implementation_details.rst`

A new page, **docs only — not the README**, laying out the strategies in sections and showing how
speed, accuracy and robustness are achieved. Add it to the toctree in `docs/source/index.rst`.

Suggested shape:

1. **The engines** — general Magnus ladder, two-flavour interaction-picture fast path, energy-
   batched separable scan, cumulative baseline scan, adiabatic+Magnus hybrid: what each assumes,
   when it applies, when it declines.
2. **Dispatch** — the order, the thresholds (`HYBRID_YIELDS_TO_CUMULATIVE_MIN_POINTS`,
   `CUMULATIVE_AUTO_MIN_POINTS`), and that the accuracy *steps* at N = 25 rather than varying
   smoothly.
3. **Speed** — with the measured numbers: 79× on scans ≥ 25 points, 0.79× of `main` on hybrid
   single points, median 2 ms across 164 Earth/solar configurations.
4. **Accuracy** — the oracle discipline (`solve_ivp` at `rtol=1e-12`, `expm` where exact, never
   one Magnus path against another), and the measured distributions.
5. **Robustness** — the safeguards and *what each cannot do*: `_profile_is_resolved` (a jump must
   be ≥ 1.33× the local smooth variation), `GAMMA_TO_ERROR`, the patch budget, and the one
   irreducible limit — a feature narrower than the probe spacing, cured by `t_breakpoints`.
6. **Warnings: what each one means and what to do about it.** A table of every warning the
   package can raise -- the condition, whether the answer is affected, and the concrete change to
   make -- plus the reasoning that a warning here is an instruction rather than a disclaimer.
   Worth its own section because it is the part of the API a user meets when something has
   already gone wrong, and because the audit in Item 5 will have produced exactly this material.
   Include the honest cases: `MagnusConvergenceWarning` reports slab width rather than accuracy,
   and `HybridCertificationWarning` means the answer is unverified, not that it is wrong.
7. **How the constants were set** — the provenance table from Item 2. This is the section that
   makes the page worth writing: it is the part no reader can reconstruct from the code.

Cross-link `/adiabatic_strategy` rather than duplicating it. **`sphinx -W` escalates any new
warning inside a `jupyter-execute` example to a build error** — `DECISION_CUMULATIVE_DEFAULT.md`
§3 records that catching a case the 665-test suite missed, so build the docs after every change.

---

## Traps — every one of these cost real time in *this* session

**Process and git**

- **A failed `cd` does not stop the rest of a `&&`-chain.** `cd nonexistent && git reset --hard
  main` ran the reset **in the main repository**, detaching thirteen commits and discarding the
  whole uncommitted working tree. Recovered from the reflog, but everything had to be re-applied.
  **Use `git -C <dir>` for every git command. Never `cd X && <destructive git>`.**
- **Commit before doing anything structural.** The above cost an hour because the work was
  uncommitted.
- **`pkill -f <pattern>` matches the process you are about to start**, and your own shell.
  Killing "pytest tests/" killed the run issued in the same command. **Kill by PID.**
- **`pgrep -f foo` matches the waiting shell's own command line**, so `until ! pgrep -f foo`
  never exits. (In the previous handover too; hit again anyway.)
- **Background waiters accumulate.** Spawning one per poll left 23 zombie `until…sleep` loops.
  Poll directly and sparsely; one long-running job needs at most one waiter.
- **Foreground commands are killed at 10 minutes.** Use `run_in_background` for the suite (~11
  min) and the batteries.
- Always run measurements under `(ulimit -v 10000000; …)`.
- `gh` is snap-confined: prefix `GIT_CONFIG_SYSTEM=/dev/null`, and pipe bodies via
  `--body-file -`.
- **Backticks in a heredoc-free `-m` message are shell-substituted.** Write commit messages to a
  file and use `-F`.

**Physics / API**

- `earth.density_matter_func_prem(r)` takes a **radius from the Earth's centre in km**; its second
  argument is `tol`, not `costhz`. Route through
  `earth.earth_radial_distance_from_depth(costhz, l/gd.UNIT_KM)`, as the package's own wrappers
  do. Getting this wrong produced 52 `ValueError`s that looked like a package defect and were
  entirely the harness's fault.
- `matter.vcc_func_from_rho_func`'s **7th positional argument is `density_is_of_number_of_electrons`,
  not `nubar`** — `harness.vcc_of` passes both by keyword so the mistake is unrepresentable. Use it.
- A `rho_func` must return a **scalar for scalar input**; a profile written array-first returns a
  0-d array and is rejected by validation. `harness.scalarize` handles it.
- `H` closures may now carry bound defaults (`def H(energy, l, VCC, _h=h)`) — that was fixed this
  session. Older documentation still steers you to factories; either works now.

**Measurement discipline**

- **Alternate, always, and carry a control.** Two workloads the change cannot touch (a vacuum
  scan and a constant-density scan) go in every timing round; if they do not come back at ~1.00×,
  the comparison is worthless. This is what let a timing run survive being taken on a busy
  machine.
- **A stalled battery is usually the oracle.** See Item 3.
- `adiabatic_propagator` at `n_points=12864` is a Python loop over `eigh` plus parallel transport
  — expensive inside a sweep. Budget for it or lower `n_points`.
- **Verify the oracle before believing a defect.** For a narrow feature an adaptive `solve_ivp`
  can step over it exactly as the detector does. `adversarial_batteries/verify_b2.py` checks a
  finding against four independent oracles; run it on anything surprising.

---

## Where I expect this to be harder than it looks

Written down in advance so "we found nothing" can be weighed against what was probed.

1. **Item 1's engines are not as independent as they look.** The cumulative scan sizes its grid
   from an `osc_prob` probe, so it inherits the general ladder's failure modes; and both share
   `magnus.magnus`. Genuine independence exists only between {hybrid} and {everything else}, and
   between all of them and `expm`. Do not report agreement between two engines that share a
   blind spot as evidence.
2. **`threshold0 = 0.1` may not be a single number.** The right value plausibly depends on the
   requested tolerance, exactly as `GAMMA_TO_ERROR` turned out to depend on the γ regime. Be
   ready for the answer to be a rule rather than a constant.
3. **Item 4's invariants will find disagreements that are correct.** `auto` and `magnus` *should*
   differ — they are different methods, and `auto` is usually the better one. The assertion is
   "within tolerance", and picking that tolerance per invariant is the actual work.
4. **5a may conclude the warning cannot discriminate.** `MagnusConvergenceWarning` reports a
   property of the slab width, not of the answer, so a genuinely predictive version may not
   exist. If so, the honest outcome is to downgrade or rename it -- but the false "can be
   ignored" clause has to go regardless, because that clause is a claim about the *answer*, and
   it is wrong.
5. **The suite is now ~11 minutes**, up from ~8, almost entirely
   `test_adiabatic_result_matches_the_engine_averaged_numerically` (6.7 s → 209 s). It asks for
   `rtol=1e-4` over five solar scale heights, where the tightened γ rule correctly forces
   refinement. Its assertion compares against a statistical reference whose own s.e.m. is far
   above 1e-4, so its tolerance is over-specified for its purpose — **loosening it is a
   reasonable call, but it is the maintainer's, not yours.**

---

## Deliberately not done

- **Items 1–5 themselves.** Only 5b was started, and nothing was committed.
- **Exposing `n_probe` through the wrappers.** Considered and rejected: `t_breakpoints` is the
  better lever for the case it would address, and is verified to work (2.9e-02 → 8.8e-04).
- **Merging #24/#25.** Both await review.
- **GitHub Pages** is still disabled, so the "Documentation Deployment" workflow fails on every
  commit. Settings → Pages → source "GitHub Actions". Nothing in the codebase can fix it.
