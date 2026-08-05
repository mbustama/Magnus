# Findings: the robustness programme (six deliverables)

**Written:** 2026-08-04, executing `HANDOVER_ROBUSTNESS_PROGRAMME.md` against `main` at
`e0b2bc9`, on branch `dev-robustness`.

**Base verified before starting**, as the brief's guard requires:
`git show origin/main:src/magnus/adiabatic.py | grep -c '^GAMMA_TO_ERROR'` → 1.
**Baseline reproduced:** 685 passed, `ruff check src/ tests/` clean, docs build under
`-W --keep-going`.

**Definition of done, held:**

| criterion | before | after |
|---|---|---|
| tests collected | 685 | **743** (+58; grown, not shrunk) |
| `ruff check src/ tests/` | clean | clean |
| `make html SPHINXOPTS="-W --keep-going"` | succeeds | succeeds |
| every new constant carries its measured population | — | `LOCAL_JUMP_RATIO`, `N_LOCAL_CONFIRM`, `MAX_LOCAL_CONFIRMATIONS`, plus provenance added for `fd_step_frac` and `threshold0` |
| every new warning's false-positive rate measured before shipping | — | 0/1440 (§3.2); the one whose rate could not be established was **not shipped** (§5c) |
| everything asserted in a docstring has been run | — | three docstring claims corrected against measurement while writing this |

---

## 1. The verdict, in one paragraph

All six deliverables are built. The programme's own premise held up better than expected:
**the instruments found three defects that the 685-test suite did not**, and two of them were
found by the instruments built for a *different* item. The cross-check diagnostic (item 1)
passes its acceptance test on the pre-fix tree — 7 of 7 known silent misses detected, each at
least four times the requested tolerance — and correctly fails to detect the one construction
the previous session recorded as unfixable. The invariant sweep (item 4) found a **pre-existing
false positive in the resolution test**, firing on 6 of 30 baselines of an ordinary scan over a
smooth Gaussian bump; that is now repaired and measured at **0 false positives over 1440
configurations**. And building item 1 exposed that **passing `cumulative` to a wrapper — even
the documented default `'auto'` — silently disabled three engines**.

The pattern the brief predicted repeated three more times: *the instrument was wrong before the
code was.* A single-slab Magnus comparison reported a composition error of 1.5; a per-position
normalisation reported a derivative error of exactly 1.0 at every step size; and a criterion
phrased as "these seven must be detected" would have failed on the fixed tree for the best
possible reason. All three were caught by looking at the number and asking whether it could be
true.

---

## 2. Item 1 — cross-method agreement instead of self-agreement

**Built.** `magnus.oscprob.cross_check_strategies(entry_point, *args, **kwargs)`. Runs whichever
engines apply to a request, reports each answer, the pairwise spread, which engines declined and
why, and which pairs are genuinely independent. Never on by default; a large spread is reported,
never raised.

**The independence question is answered explicitly, not assumed.**
`magnus.oscprob.ENGINE_FAMILIES` records that the general ladder, the cumulative scan and the
separable scan are one family — they share `magnus.magnus`, and the cumulative scan additionally
*sizes* its grid from an ordinary `osc_prob` probe — so their agreement is not evidence. The
result carries `max_spread_independent`, restricted to cross-family pairs, which is the number to
read. This is the brief's caution #1, encoded rather than remembered.

**Acceptance, against the pre-fix tree (`978663a`).**
`docs/dev/adversarial_batteries/crosscheck_acceptance.py`, run under both trees via `PYTHONPATH`.
It deliberately does not import `cross_check_strategies` — that function does not exist on the
pre-fix tree, and an acceptance test that only runs where the fix already is tests nothing.

| construction (FINDINGS §3) | silent error there | max cross-family spread |
|---|---|---|
| step function, unmarked edge | 5.395e-01 | **5.399e-01** |
| ten crossings (worst found anywhere) | 3.907e-02 | **3.913e-02** |
| sinusoid at span/7 | 1.672e-02 | **1.687e-02** |
| kink, C⁰ but not C¹ | 1.448e-02 | **1.448e-02** |
| singularity approached | 8.625e-03 | **8.613e-03** |
| sub-threshold bump, w = 1e-2 span | 7.701e-03 | **7.768e-03** |
| sub-threshold bump, w = 3e-2 span | 4.388e-03 | **4.594e-03** |
| narrow bump, w = 3e-5 span (§8.3) | 2.907e-02 | 3.5e-14 — **not detected** |

**Seven of eight, each at least four times the requested 1e-3.** The eighth is the documented
irreducible limit and is carried as a row *expected* to go undetected: when every engine is wrong
together there is nothing to disagree about.

**The pass criterion was restated to be tree-independent**, which matters more than it sounds.
Phrased as "these seven must be detected" it passes on the pre-fix tree and fails on the fixed
one — because the defects are gone. Phrased as **"whenever at least one engine is outside
tolerance, the maximum cross-family spread must be at least that tolerance"** it is meaningful on
both, and it states the instrument's actual reach: *it sees a wrong engine exactly when some
other engine got it right.* Both trees now report 0 failures of that criterion.

**Wired into CI** as `tests/test_engines.py`, over a four-configuration profile matrix, asserting
cross-family agreement within 5e-3 — a bound taken from the measured sweep (worst 3.8e-03,
between the hybrid strategy and the energy-batched scan, both inside their own requested 1e-3 of
the truth) rather than chosen.

### 2.1 A defect found while building it

`cross_check_strategies` reported that **no engine had answered at all**. The cause was not in
the diagnostic:

```
osc_prob_matter_std_potential(..., )                    -> hybrid,  P_ee = 0.2560247069
osc_prob_matter_std_potential(..., cumulative='auto')   -> magnus,  P_ee = 0.2560154231
```

`'auto'` is the **documented default** for `cumulative`. Passing it explicitly changed which
engine answered and moved the answer by 9.3e-06. The mechanism is the exact trap the brief
warned about for item 5b, already live in the shipped code: `_resolve_cumulative_kwarg` pops
`cumulative` out of `**kwargs` at the *last* call site, so the hybrid, interaction-picture and
separable dispatchers saw it first — and every one of them declines outright on an unrecognised
key.

Fixed by resolving `cumulative` at the top of each wrapper, before `scan_kwargs` is built. This
required settling what the flag means, and the answer is not symmetric:

* `cumulative=False` disables **the cumulative scan only**. It never meant "and also the
  adiabatic strategy"; `strategy='magnus'` is the opt-out for that, and it already implies
  `cumulative=False`.
* `cumulative=True` names one engine and is documented to **raise** rather than fall back, so
  the three engines tried before it now stand aside for it. Otherwise it would be silently
  substituted — the one thing that flag exists to rule out.

`osc_prob_sun`/`osc_prob_earth` had the same defect through `_osc_prob_with_potential`; fixed the
same way. Pinned by two tests in `tests/test_engines.py`.

---

## 3. Item 4 — property tests for cross-entry-point agreement

**Built.** `docs/dev/adversarial_batteries/invariants.py` measures; `tests/test_invariants.py`
asserts. The split is deliberate: some of these disagreements are *correct* — `auto` and `magnus`
are different methods — so each bound comes from a measured distribution, recorded in the test's
own docstring, rather than from a guess.

Two invariants are exact and asserted as `== 0.0` (shuffled baseline order, `n_jobs>1` vs 1);
the rest carry measured bounds. Over 60 configurations (5 profile families × d ∈ {2,3} × 2
energies × N ∈ {1, 8, 30}), before and after the resolution-test repair of §3.2:

| invariant | median | p90 | max, before | max, after |
|---|---|---|---|---|
| `auto` vs `hybrid` | 0.00e+00 | 2.1e-04 | 7.36e-02 | **6.21e-04** |
| `auto` vs `magnus` | 2.0e-04 | 1.4e-03 | 6.48e-03 | 6.48e-03 |
| `cumulative` True vs False | 1.2e-04 | 1.2e-03 | 3.05e-03 | 3.05e-03 |
| `auto` vs `cumulative=True` | 1.7e-09 | 2.1e-04 | 4.80e-04 | 4.80e-04 |
| a scan vs the same points singly | 1.9e-09 | 2.6e-04 | 1.10e-02 | **6.21e-04** |
| shuffled baseline order | 0.00e+00 | 0.00e+00 | **0.00e+00** | **0.00e+00** |
| `n_jobs=2` vs 1 | 0.00e+00 | 0.00e+00 | **0.00e+00** | **0.00e+00** |
| composition, hybrid | 4.8e-08 | 4.3e-04 | 1.30e-03 | 1.30e-03 |
| composition, general ladder | 1.1e-08 | 1.3e-05 | *harness error* | **8.79e-04** |
| unitarity (row/column sums) | 5.0e-13 | 8.7e-12 | 2.72e-11 | 2.72e-11 |

The two rows that moved are the two the repair touched, and they moved by 118× and 18×. The
rest are bit-identical before and after, which is the statement that the repair changed nothing
it should not have.

### 3.1 The sweep found a pre-existing false positive

The worst disagreement in the first sweep was **7.4e-02**, on a smooth Gaussian bump of width
1e-2 of the trajectory at d = 3, 50 MeV, N = 30. Scored against `solve_ivp`:

| path | engine | error | warnings |
|---|---|---|---|
| default (`auto`) | cumulative | **3.5e-09** | — |
| `strategy='hybrid'` | hybrid | **7.4e-02** | HybridCertification, UnmarkedDiscontinuity |
| `strategy='magnus'` | magnus | 4.2e-05 | MagnusConvergence |

Not a silent miss — the forced hybrid path warns twice — but the *reason* it declined was wrong:
`adiabatic._profile_is_resolved` declared a C^∞ Gaussian discontinuous, on **6 of the 30
baselines** of one ordinary scan.

**Mechanism.** The test is a maximum over probe intervals of "what fraction of this interval's
variation falls in one half". At a smooth **turning point** the two halves are genuinely
asymmetric — one nearly cancels — so that single interval's ratio is a draw in [0.5, 1] depending
on where the extremum happens to fall inside it, and one interval decides for the whole profile.
Refining the grid does not remove it; it *re-draws* it, which is why the existing two-stage
coarse-then-fine protocol could not separate a turning point from a jump either.

Measured at the offending interval: max ratio 0.7518 at n_probe = 200, at position 0.6985 of the
span, with the bump centre at 0.7035 — the interval containing the peak.

**Why the original measurement missed it.** `RESOLUTION_RATIO`'s population swept profile
families × dimensions × energies, always on the **full** trajectory. A baseline scan calls the
test once per baseline — thirty different sub-intervals of the same profile — and the statistic
depends on where the extremum falls relative to the grid. The axis that mattered was the one the
population did not have. This is the same shape of mistake as the original `GAMMA_TO_ERROR`
error, on a different axis.

### 3.2 The repair, and its measurement

A flagged interval is now a *candidate*, not a verdict. Each of the worst few flagged intervals
is re-sampled alone on `N_LOCAL_CONFIRM = 33` points, and the interval is a genuine
discontinuity only if the largest single adjacent step still carries more than
`LOCAL_JUMP_RATIO = 0.5` of its variation. A jump is not diluted by refinement — one step still
carries all of it; a smooth feature spreads over every step.

Measured (`docs/dev/adversarial_batteries/resolution_fp.py`), swept over **sub-intervals**, the
axis the earlier work lacked — 10 smooth and 3 piecewise families × d = 2…5 × 3 energies × 12
sub-intervals:

| population | flagged intervals | local statistic | declared unresolved |
|---|---|---|---|
| smooth, 1440 configurations | 79 | **≤ 0.087** | **0 / 1440** |
| piecewise, 432 configurations | 348 | **≥ 1.000** | 348 / 432 |

The separation is a factor of 5.7 below the threshold and 2.0 above it. The 84 piecewise
configurations not flagged are sub-intervals that end before the first jump — they contain no
discontinuity, so resolving them is correct.

Costs nothing on an ordinary call: the confirmation runs only on intervals the cheap test already
flagged, and a profile that flags none never reaches it.

### 3.3 A harness error, recorded because it looked exactly like a defect

The first sweep reported a composition failure of **1.5** for the general Magnus path —
`U(0→L) ≠ U(L/2→L)·U(0→L/2)`, which would be a time-ordering bug of the first order. It was the
harness: both sides were computed on a **single slab** over their whole interval, so neither was
converged, and two unconverged integrations of course disagree. With a converged grid on both
sides the row behaves. Third instance this session of "verify the oracle before believing a
defect".

---

## 4. Item 5 — warnings

### 5a. The false clause is gone

`MagnusConvergenceWarning` ended with: *"If a target tolerance (rtol/atol) was requested, the
adaptive refinement narrows the slabs automatically and this warning can be ignored."*

That is false in exactly the cases where the warning matters — measured on a sawtooth density at
`rtol=atol=1e-3` explicitly requested, the refinement ran and the answer was still 7.484e-03. It
is replaced by the four-part standard: what was detected (with **how far past π**, in three
buckets, so the message stays static enough for Python's once-per-session filter), what it means
(**unknown** — this reports a slab width, not an error, and it fires on answers accurate to
1.6e-06 as well), what to change (narrower slabs; `t_breakpoints` at any jump; *not* a higher
order), and when it is safe to ignore (when the answer has been checked another way).

The brief anticipated that 5a might conclude the warning cannot discriminate. **It cannot**, and
that is now what it says. Downgrading or renaming it was considered and rejected: the condition
it reports is real and actionable, and the class name is in users' warning filters.

### 5b. `strategy_info`

Added as an **explicit named parameter** to the three scenario wrappers, which sidesteps the
plumbing trap the brief describes entirely: a named parameter never enters `**kwargs`, so there
is no layer at which it must be popped. It reports the engine that answered, its family, the
hybrid certification flag, every engine that declined and why, and the full dispatch trace.

```
solar single point                 engine=hybrid      certified=True
solar single point, magnus         engine=ip_exp      certified=None
solar 30-pt scan                   engine=cumulative  certified=None
unmarked step (hybrid declines)    engine=magnus      declined=[('hybrid', 'the profile is
                                                        not resolved at the probe scale')]
unmarked step, strategy=hybrid     engine=hybrid      certified=False
```

This is also the brief's preferred answer to 5c's second candidate ("the hybrid strategy
declining under `auto` is entirely silent"): an opt-in, not a warning.

The same trace drives `cross_check_strategies`, so one instrument serves both items.

### 5c. New warnings shipped, and one not shipped

**Shipped: `UnmarkedDiscontinuityWarning` on the hybrid path.** The detector already ran there;
failing it caused a silent decline, so the caller heard about slab widths from whichever engine
answered instead — true, and pointing at the wrong knob. False-positive rate measured **before**
shipping, over the sub-interval axis: 0 of 1440 smooth configurations (§3.2). Getting to that
number required the detector repair; shipping the warning on the old detector would have fired on
6 of 30 baselines of an ordinary smooth scan.

**Shipped: a ceiling-derived accuracy grid in the cumulative scan.** When the probe that sizes
the grid hits `max_n_slabs`, `ToleranceNotAchievedWarning` already fires — but it describes *the
probe*, one call whose probabilities are discarded. The consequence that matters is larger: every
baseline in the scan inherits a grid sized by a cap. Fires exactly when
`probe n_slabs >= max_n_slabs`, which is a fact rather than a heuristic; verified to fire at
`max_n_slabs=200` and stay silent at the default.

**Not shipped: marginal certification.** `hybrid_propagator` certifies when
`GAMMA_TO_ERROR * gamma_max <= atol + rtol`, and a result passing at 0.95 of that bound is
certified on the same footing as one at 0.01, while the constant is itself good only to ~2×. A
warning here would fire on a *correct* result whose margin is thin, and nothing measured
distinguishes those from the rest — the population that would trigger it has no measured error
excess. Shipping a warning whose false-positive rate is unknown and probably high would fail the
brief's own bar. The margin is instead **exposed**: `hybrid_propagator(info=...)` returns
`gamma_max`, so a caller who cares can compute the ratio.

### Also improved

* `ToleranceNotAchievedWarning` now reports **how far from converged** the ladder stopped, as a
  multiple of the requested tolerance, in three buckets. The disagreement between the last two
  levels is computed anyway by the comparison that decides convergence; it was simply discarded.
* `HybridCertificationWarning` now says **unverified is not wrong**, and names three concrete
  actions (`strategy='auto'`, `t_breakpoints`, a looser tolerance) where it previously suggested
  none.

---

## 5. Item 3 — fuzzers in CI, on aggregate assertions

`tests/test_fuzz_statistics.py`. Two populations, both seeded, asserting on the **distribution**
— silent-miss rate, median, worst case — never on individual cases.

| population | n | median | p90 | max | outside 1e-3 | silent |
|---|---|---|---|---|---|---|
| random smooth (Fourier sums), `solve_ivp` oracle | 24 | 7.8e-09 | 4.9e-04 | 1.15e-03 | 2 | **2 (8.3 %)** |
| random piecewise-constant, `expm` oracle (exact) | 40 | 4.3e-04 | 6.9e-03 | 1.28e-02 | 10 | **0** |

The two smooth silent misses sit at 1.05e-03 and 1.15e-03 against a requested 1e-3 — the residue
`FINDINGS` §9.4 describes, not a new failure. At n = 24 a 4.1 % rate has an expectation of one,
so two is noise.

**The oracle trap the brief flags is designed around rather than survived.** The smooth
population is held to 20–200 MeV and d ∈ {2,3} because `solve_ivp` dominates the cost at low
energy and high flavour count; and the piecewise population uses `expm` composed across segments,
which is the *exact* operator for those Hamiltonians, so it costs nothing and cannot itself step
over a jump.

---

## 6. Item 2 — audit of the remaining constants

### `fd_step_frac = 1e-6`

Scored against the **analytic** `dH/dl` — available because
`H = h_vac/E + C·n_e(l)·P_ee` and every profile used has a closed-form `n_e'`, so the reference
carries no error of its own. Relative to the largest `|dH/dl|` on the path:

| profile | optimum | error at 1e-6 | error at the optimum |
|---|---|---|---|
| solar exponential | 1e-5 | 6.8e-11 | 1.9e-11 |
| sinusoid, period span/7 | 1e-6 | 3.4e-10 | 3.4e-10 |
| Gaussian bump, w = 1e-2 span | 1e-7 – 1e-8 | 2.5e-09 | 2.6e-11 |

**The optimum is not a single number** — it tracks the profile's shortest length scale, as the
truncation/cancellation trade predicts. But the curve is flat: anywhere in 1e-8 … 1e-5 the error
stays below 3e-09 relative on every profile measured, six orders below anything that could move a
probability at these tolerances. Outside that band it degrades fast in both directions (1.6e-04
at 1e-12 from cancellation, 0.23 at 1e-2 from truncation). **The band, not the value, is what to
preserve**, and that is now in the docstring.

*The first version of this measurement reported a relative error of exactly 1.000 at every step
size on the Gaussian profile.* It normalised per position, and on a narrow bump almost every
sampled position has `dH/dl ≈ 0`, so it divided ~0 by ~0. The instrument saturated before the
code did — second instance this session.

### `threshold0 = 0.1`

Measured as a cost-against-accuracy trade at three requested tolerances; see §7 below for the
result and the reframing it forced.

### The rest

Listed explicitly as **unmeasured, with the reason**, in
`docs/source/implementation_details.rst` — which is the brief's pass criterion. Four of them
(`max_n_probe`, `max_n_points`, `max_num_loops`, `max_iters`) are cost ceilings rather than
calibrations: reaching one is *reported* by `ToleranceNotAchievedWarning` rather than absorbed,
so "unmeasured" means something milder for them than for a threshold that silently decides an
outcome.

---

## 7. Item 6 — `docs/source/implementation_details.rst`

Added and in the toctree. Seven sections as the brief specifies: the engines and what each
declines; dispatch, with the N = 25 accuracy step stated rather than left to be discovered;
speed, with the alternating-harness numbers and the note that absolute times from a loaded box
are worthless; accuracy, with the oracle discipline and the measured distributions; robustness,
with **what each safeguard cannot do**; a warning table giving the condition, whether the answer
is affected, and the concrete change; and the constant provenance table, including the
unmeasured ones.

No `jupyter-execute` blocks, deliberately: the brief records that a new warning inside one is a
`sphinx -W` build *error*, and this tranche adds warnings. The page builds clean under
`-W --keep-going`.

---

## 8. What this did not cover

Stated so that "we found nothing" can be weighed against what was probed.

* **Timing was not re-measured.** Every change here is either a diagnostic that is off by
  default or a warning; the one change on a hot path is the local confirmation in
  `_profile_is_resolved`, which runs only on intervals the cheap test already flagged and
  therefore not at all on any smooth profile.
* **The fuzz populations are smaller than the batteries'** (24 and 40 against 150 and 250), by
  design, to fit CI. The batteries remain the place to run the large populations.
* **`d = 4` and `d = 5` are absent from the CI fuzz population** and present only in the
  resolution-test sweep. The adversarial validation found the failure rate flat in dimension, so
  this is a cost decision rather than a coverage claim.
* **The generic entry points (`osc_prob_sun`, `osc_prob_earth`) got the `cumulative` fix but not
  `strategy_info`.** They take a user-supplied Hamiltonian and do not build `scan_kwargs`, so the
  out-parameter would need separate plumbing; the engine trace already covers them, so this is
  work remaining rather than a decision.
* **`strategy='magnus'` still does not imply `cumulative=False` on the generic entry points**,
  where the three scenario wrappers do imply it. Noticed while fixing the kwarg shadowing;
  changing it would move answers for `osc_prob_sun(strategy='magnus')` scans, which is a
  maintainer's call rather than a bug fix.

---

## 11. Second tranche: the seven exposures §8 left open

Written the same day, after the §8 list was read back as "how robust is it now" and every item on
it was commissioned. **Three of the seven were closed, two were closed as measured rejections
(built, tested, removed), and two were closed by measurement alone.** The rejections are the
more useful half of this section.

Bit-identity was checked after every change: **0 of 11 workloads moved**, except where noted.

### 11.1 The narrow-feature blind spot is now detected — CLOSED

The one exposure the adversarial validation could not close: a feature narrower than every grid
the package lays down, wrong by 2.9e-02 on all engines at once, silently.

`adiabatic.find_hidden_features` looks at **the profile** rather than at the answers, which is
what lets it reach a class no cross-check can. Within each interval of a reference grid (the
refinement ceiling, `max_n_probe = 6400`), it compares the total variation a denser grid sees
inside that interval with the change its endpoints show; the excess is variation hidden between
reference samples, and the statistic is the largest such excess as a **fraction of the total**.

**The first version was wrong, and the way it was wrong is the point.** `TV_dense/TV_reference`
sends the denominator to zero on a sinusoid at exactly the probe spacing and the ratio to
1e13 — on a profile the package answers to ~1e-11. What separates a hidden bump from an aliased
sinusoid is not *how much* is hidden but *where*: the sinusoid hides some in every interval
(share ~1/n_ref), a bump hides all of it in one (share ~1).

| population | concentration |
|---|---|
| 67 smooth/resolvable profiles — solar, multi-resonance, noisy, sinusoids at 1×/2×/½× the probe spacing, 400 crossings, a declared step, 30 random Fourier sums, 30 random-width bumps | max **0.060** |
| features in the unresolvable band | 0.91 – 1.00 |

**0 false positives at every threshold from 0.2 to 0.6.** Shipped at 0.3, five times the
measured ceiling. Detection over 60 random positions per width: **0.68 at 3e-5, 0.90 at 1e-5,
0.82 at 3e-6, 0.73 at 1e-6.**

That is most of the class, not all of it, and the shortfall is structural: 3e-5 sits where the
refinement ceiling can *partly* resolve the feature, so the statistic is diluted; 1e-6 is far
below the dense spacing, so whether a sample lands inside is luck. Against a prior state of
**zero** detection, that is the improvement on offer, and `HIDDEN_FEATURE_CONCENTRATION`'s
docstring says so rather than implying a guarantee.

**Cost: 0.37 ms**, about 3 % of a 13 ms single-point call, and inside the noise floor of an
alternating measurement (controls at 1.03–1.07×). Three things made it affordable: the statistic
is provably **identical on the scalar potential** (H is affine in V_CC — verified bit-for-bit at
d = 2, 3, 5) which is 18× cheaper than sampling H; it runs **once per call**, not per (energy, L)
point, since the profile does not depend on energy; and 8 sub-steps per reference interval
measured as good as 32 for the statistic while the cost past 8 goes superlinear.

The action is `HiddenFeatureWarning`, naming the position and the exact `t_breakpoints` to pass.
Auto-inserting them was measured and **not** shipped: it improves the answer 3–46× (3.9e-03 →
8.5e-05) and stops it being silent, but it is a partial cure that also changes dispatch, and
choosing a grid is the caller's call. The numbers are in the warning so the choice is informed.

### 11.2 Cross-checking the default path — BUILT, MEASURED, REMOVED

The headline request, and it does not work. Below the N = 25 seam, `strategy='auto'` would
verify a window-free hybrid result against the general Magnus ladder, since certification with
no window rests on γ alone and `GAMMA_TO_ERROR` is good only to ~2×.

| what was measured | result |
|---|---|
| 200 random smooth profiles, shipped constant | 25 window-free certified, ladder agreed with **all 25** |
| the same, `GAMMA_TO_ERROR` made optimistic by 2× | 41 certified, 3 genuinely wrong — check fired **0 times** |
| non-circular trigger: verify **every** window-free result | still fires **0 times**, still misses the same 3 |

The middle row is the instructive one. The first trigger was `GAMMA_TO_ERROR·γ/(atol+rtol) >
0.3` — computed **from the very constant it was insuring against**, so mis-calibrating the
constant shrank the trigger in step with it. A self-referential check: exactly the failure shape
this whole programme exists to find, reproduced by the person writing the fix, and caught only
because the test deliberately broke the constant.

Removing the circularity does not rescue it, and the last row says why. **What is left in the
weak band is not disagreement between engines — it is the engines being wrong together.** A
cross-check detects the former by construction and can never detect the latter. It is the same
structural limit that stops `cross_check_strategies` seeing a sub-probe feature, and the reason
§11.1's instrument looks at the profile instead of at the answers.

Cost, had it shipped: 9 % of random calls, and 100 % of ordinary solar single points, which are
window-free (margin 0.057 at 2ν, 0.123 at 3ν). Zero measured benefit against that is not a
trade; the code is gone and the finding is recorded beside the dispatch constants.

### 11.3 `threshold0` as a rule — BUILT, MEASURED, REVERTED

Section 6 concluded the right value was a rule rather than a constant, and this tranche built
it: start at `(atol+rtol)/GAMMA_TO_ERROR`, the γ at which certification actually flips.

| workload | `t0 = 0.1` | the rule | |
|---|---|---|---|
| single point, solar | 1.624e-06 | 1.184e-10 | **13711× better** |
| sub-threshold scan, N = 8 | 3.220e-05 | 3.814e-05 | 1.2× worse |
| **energy scan at fixed baseline** | 2.509e-05 | **4.954e-04** | **20× worse** |

All inside 1e-3, but 4.95e-04 spends half the budget where 2.5e-05 spent a fortieth. Mechanism:
starting low opens a window on the first iteration, and `windows_next or windows_prev`
short-circuits the γ check, so agreement is accepted at a **coarser** transport grid.

**The sweep that justified the rule ran at a fixed baseline; the row that refuted it was an
energy scan.** That is precisely the mistake §6 congratulated itself on avoiding, committed one
section later. Reverted; the measurement and the reason are in `adiabatic.THRESHOLD0_PROVENANCE`.
The bit-identity check is what caught it — the three rows it moved were all hybrid-path rows.

### 11.4 The seven unaudited constants — CLOSED

Swept across **18 workloads spanning single points, baseline scans and energy scans** × 3 profile
families × d = 2, 3, with the reference computed once per workload
(`constants_audit2.py`). Worst error over all workloads, default in bold:

| constant | sweep | worst error across the sweep |
|---|---|---|
| `n_probe0` | 50 … 800 | 4.98e-04, 5.22e-04, **4.49e-04**, 3.38e-04, 3.36e-04 |
| `n_points0` | 51 … 801 | **4.49e-04** at every value, identical to three digits |
| `min_threshold` | 1e-4 … 1e-8 | 4.49e-04 at every value — **never reached** |
| `patch_atol` | 1e-5 … 1e-9 | 4.49e-04, 4.49e-04, **4.49e-04**, 3.04e-04, **2.08e-02** |
| `n_slabs0` | 100 … 1600 | 4.49e-04 at every value |
| `growth_factor_n_slabs` | 1.2 … 3.0 | 4.49e-04 at every value |
| `min_n_tpts_per_slab` | 2, 4, 8 | 4.49e-04 at every value |

Six are **not load-bearing**: they set where a doubling ladder starts, and the ladder reaches the
same place regardless. `min_threshold` is different and is written up as such — the sweep shows
the ladder never reaches the floor on any measured workload, which is evidence about the
population, not about the constant.

`patch_atol` at 1e-9 is the one real finding: most rows improve sharply, but one energy scan goes
to **2.08e-02**. The cause is not the constant — at 1e-9 the patch cannot converge within
`max_n_slabs`, the hybrid declines, and the **energy-batched separable engine** answers and is
that much worse on that profile. It warns, so it is loud rather than silent. What that row
actually measures is **fallback quality**, and it is left as an open question rather than chased
here.

### 11.5 Warning false-positive rates — CLOSED

`warn_fp.py` was killed unfinished last time because it put `solve_ivp` on every case. Rebuilt
with a split oracle — `expm` (exact, free) for piecewise profiles, `solve_ivp` only for smooth
families at ≥ 30 MeV and N ≤ 8 — it completes **168 of 168**.

| warning | fired | true positives | false positives | FP rate |
|---|---|---|---|---|
| `MagnusConvergenceWarning` | 70 | 17 | 53 | **76 %** |
| `UnmarkedDiscontinuityWarning` | 56 | 23 | 33 | 59 % |
| `ToleranceNotAchievedWarning` | 37 | 16 | 21 | 57 % |

Silent misses across the whole population: **2 of 168 (1.2 %)**.

The 59 % on `UnmarkedDiscontinuityWarning` needs reading carefully and its docstring now says so:
it reports a *condition about the input*, and on every one of those 33 the condition was real —
there was an undeclared discontinuity — and the answer survived anyway. Declaring the edges would
still have improved it by orders of magnitude.

**The mechanism measurement settles what §5a could only assert.** Of 66 single-point calls, some
refinement level exceeded π in 46 — but the level whose answer was *returned* did so in only
**7**. So **39 of 46 firings (85 %) describe an intermediate grid nobody receives**. Keying the
warning to the returned level would cut false alarms from 31 to 5. That change is mechanical and
is **deliberately not made here** — it touches the refinement loop and the warning plumbing
several tests depend on — but it is now specified with its numbers.

### 11.6 Fuzz power and d = 4, 5 — CLOSED

The piecewise oracle is `expm` composed across segments: exact, and free. So that population
carries the power and the flavour coverage at no cost.

| population | before | after |
|---|---|---|
| smooth (solve_ivp), d ∈ {2,3} | n = 24, 2 silent (8.3 %) | **n = 40**, 1 silent (2.5 %) |
| piecewise (expm), d ∈ {2,3} → **{2,3,4,5}** | n = 40, 0 silent | **n = 120**, 0 silent |

At 0 of 120 the 95 % upper bound on the piecewise silent-miss rate is about **2.5 %**. Cost: the
file takes about 6 minutes, nearly all of it the 40 smooth cases, and the docstring states that
so trimming `N_SMOOTH` is an informed decision rather than a discovery.

### 11.7 What is still open after this tranche

* **Features narrower than ~2e-5 of the trajectory** are detected only 73–82 % of the time, and
  the detector reports rather than cures. The residual is structural, not a tuning matter.
* **Fallback quality** (§11.4): when the hybrid strategy declines on an energy scan, the engine
  that answers can be two orders worse. Found incidentally; not characterised.
* **`MagnusConvergenceWarning`'s 85 % intermediate-level noise** is specified but not fixed.
* **`min_threshold`** governs a regime no measured workload enters, so it remains unexercised
  rather than validated.

---

## 12. Third tranche: the six items §11.7 left open

Commissioned as "implement 1 to 6, act as your own reviewer". **Two closed with real fixes, two
closed as measured rejections, one closed by construction, one closed by deletion.** The
reviewer pass found two defects in my own new code, which are recorded here rather than quietly
patched.

### 12.1 Fallback quality — CLOSED, and it was the dispatch order

The open item read "when the hybrid path declines on an energy scan, the separable engine can be
two orders worse; found incidentally, not characterised". Characterising it
(`fallback_quality.py`, 42 workloads across 7 profile families, d = 2 and 3, single points,
baseline scans and energy scans, every applicable engine forced and scored) found something
larger and closer to home.

| | |
|---|---|
| `'auto'` more than 10× worse than the best engine that applied | **30 of 42** |
| worst factor | **900 000×** |
| `'auto'` outside the requested tolerance | 3, of which **2 silent** |
| worst error, adiabatic hybrid, 42 workloads | 1.68e-03 |
| worst error, **cumulative scan, the 28 it serves** | **1.13e-07** |

`'auto'` picked the hybrid path on all 42. The cumulative scan was available on 28 of them and
is three to six orders more accurate — and it was declined every time by
`HYBRID_YIELDS_TO_CUMULATIVE_MIN_POINTS = 25`.

**That constant's docstring justified 25 on two claims, and both were wrong outside solar.** It
said yielding earlier would cost "several times slower (7.6× at N = 2) to buy accuracy that was
already two orders inside what the caller asked for". Re-measured:

* the accuracy given up is not two orders inside tolerance — it is up to 900 000×, and
  **outside** the tolerance, silently, on two of 42 workloads;
* the cost is not several times slower. Alternating with a control that returned 0.99×, the
  median cumulative/hybrid cost ratio is **0.87× at N = 2, 0.48× at N = 4, 0.25× at N = 8** —
  the cumulative scan is *cheaper at every size measured*, and 30× cheaper on a multi-resonance
  profile at N = 8. The 7.6× figure was solar, at N = 2, which is the best case for the hybrid
  path and the worst for the cumulative scan's strict probe.

**The seam is now 8**, chosen where the worst case stops mattering: at N = 8 the only profile
that pays anything is solar at d = 3, at 1.44×, against three to six orders of accuracy; at
N = 4 that worst case is 2.84× and at N = 2 it is 5.75×, which is a real price on the cheapest
requests. Both silent misses were baseline scans at N = 8 and route to the cumulative scan now.

This also corrects §11.2. The weak-band cross-check failed partly because it verified the hybrid
path against the **general ladder**, which shares the failure; the **cumulative scan** disagrees
with it by 1.7e-03 on exactly those cases. `ENGINE_FAMILIES` groups them together on shared
machinery, which is right about their blind spots and wrong about their accuracy.

### 12.2 `MagnusConvergenceWarning` keyed to the returned level — BUILT, MEASURED, REVERTED

85 % of its firings describe an intermediate grid nobody receives, so keying it to the returned
level looks obviously right. Implemented, then re-measured over the same 168 configurations:

| | before | after |
|---|---|---|
| fired | 70 | 53 |
| **true positives** | **17** | **4** |
| false positives | 53 | 49 |
| silent misses in the population | 2 | 2 |

It removed 13 true positives to remove 4 false ones. **"The ladder started far from convergence"
predicts a bad answer better than "the final grid is coarse" does.** The 85 % statistic was
real; the inference from it was not. Reverted, and the rates were re-measured afterwards to
confirm they returned exactly to 70/17/53. The mechanism is kept, private and documented, as
`magnus._deferred_slab_norm`.

### 12.3 Detection at the narrow end — CLOSED by understanding the two mechanisms

Measuring the *distribution* of the misses, rather than only their rate, split them cleanly:

| feature width | miss rate | median concentration of the misses |
|---|---|---|
| 3e-5 | 26 % | 0.100 (max 0.287) — **near-threshold** |
| 1e-5 | 12 % | 0.000, p90 0.168 — mixed |
| 3e-6 | 14 % | **0.000** — total |
| 1e-6 | 35 % | **0.000** — total |

Near-threshold misses at 3e-5 are the reference grid partly resolving the feature, which is
honest — 3e-5 is where the refinement ceiling *can* half-see it. Total misses are simply
unsampled, and no threshold change reaches them: **you cannot detect what you never sample.**

So the only lever is sampling density, and it is now spent where the call can afford it: the
scan runs once per call regardless of point count, so `n_sub` scales with the request — 8
sub-steps (0.37 ms) for a single point, 32 (2.85 ms) for sixteen or more, holding the scan under
about 7 % at every size instead of 20 % of the cheapest one. A single point keeps the cheapest
scan **by design**: what finer sampling buys is widths of 3e-6 and below, narrower than anything
physically plausible in a density profile.

### 12.4 `min_threshold` — CLOSED by constructing its regime

It was recorded as "unexercised rather than validated". The regime was then derived rather than
searched for: the floor is reached only when γ_max is *below* it (so no window can open however
far the threshold falls) **and** the tolerance is tighter than `GAMMA_TO_ERROR × γ_max` (so the
γ rule cannot certify either). An almost-flat profile (γ_max = 3e-7) at `rtol = atol = 1e-9`
satisfies both:

| `min_threshold` | error | windows | iterations | time |
|---|---|---|---|---|
| 1e-4 | 8.49e-13 | 0 | 9 | 3.3 s |
| **1e-6** | 8.49e-13 | 0 | 13 | 7.4 s |
| 1e-8 | 3.24e-12 | 1 | 13 | 7.7 s |
| 1e-10 | 3.24e-12 | 1 | 13 | 7.9 s |

It does change behaviour there — below γ_max a window opens — but not usefully:
`certified=False` at every value, the error is three orders inside the requested tolerance
either way, and the window costs 2.4× the time and makes the answer very slightly *worse*.

### 12.5 Suite time — CLOSED, partly

Measured rather than estimated: `test_smooth_profile_fuzz_statistics` 437 s,
`test_piecewise_profile_fuzz_statistics` 80 s, `test_fuzzing_raises_nothing` **72 s**.

The third was pure duplication — it re-ran 16 of the same cases to assert that nothing raises,
which the other two establish anyway by running 160 cases that cannot pass if one raises. Folded
into `collect` and deleted. The smooth population is left at 40: it is the expensive one, and
cutting it is exactly the statistical power the previous tranche added. The cost is stated in
the file's docstring so trimming `N_SMOOTH` stays an informed decision.

### 12.6 `strategy_info` on the generic entry points — CLOSED

`osc_prob_sun` and `osc_prob_earth` now accept it, through `_osc_prob_with_potential`, with the
same keys as the three scenario wrappers. A user-supplied Hamiltonian gets the same answer to
"which engine answered, and what stood aside" as a built-in scenario does.

### 12.7 What the reviewer pass found in my own new code

* **`find_hidden_features` crashed on a profile that returns a bare scalar** for array input —
  `values.shape[0]` on a 0-d array. A diagnostic that breaks the call it was inspecting is worse
  than no diagnostic, and one of my own tests claimed this case was handled (it passed only
  because the profile it used raised rather than returning a scalar).
* **A non-finite profile produced `concentration = nan`**, which compared `False` against the
  threshold by luck rather than by design and would have leaked a nan into `strategy_info`.

Both now refuse quietly, and five parametrised cases pin it.

* **Dead code**: `hybrid_kw` in `constants_audit2.py`, written and never called.
* **Public API for a rejected design**: `deferred_slab_norm` was exported before it was measured
  and reverted; it is now private.

### 12.8 Still open after this tranche

* **Features below ~2e-5 of the trajectory**: detected 73–82 % on a single point, and reported
  rather than cured. Structural.
* **`ENGINE_FAMILIES` groups the cumulative scan with the general ladder.** Right about shared
  blind spots, wrong about accuracy — §12.1 shows they differ by four orders on the same
  request. A cross-check between them is more informative than the grouping implies.
* **The suite is still long**; the smooth fuzz population is the remaining cost and trimming it
  trades away statistical power.

---

## 13. Fourth tranche: does any of this reach a real user?

Commissioned as `HANDOVER_PHYSICAL_PROFILES.md`. Every catastrophic number in §§1–12 came from
either a deliberate adversarial construction or from the fuzzers' random Fourier sum, which is a
convenient way to make a smooth positive function and has no physics in it. So the robustness
case rested on profiles nobody would compute. This tranche built a population a referee would
accept as physically motivated and re-ran the instruments against it.

**The answer is not the "mostly no" the brief expected, and the two families that pay are not the
ones it predicted.** A real published solar model produces a silent miss; the interpolation-kink
family, which the brief expected to pay, produces nothing at all.

### 13.1 The population, and how physical each family actually is

`adversarial_batteries/physical_profiles.py`, 19 families. §7.3 of the brief asks each to be
explicit about where it sits on the "physically motivated" spectrum, because a real tabulated
model is not the same evidence as an analytic shock with plausible parameters:

| family | n | provenance |
|---|---|---|
| tabulated, linear and cubic, N = 20…5000 | 7 | the package's own exponential, resampled. The *shape* is the package's; what is physical is the user behaviour — loading a table and interpolating it |
| BS05(AGS,OP) solar model | 2 | **a real published model, used directly** — the only family that is a measurement rather than a parametrization |
| SN shock, width 1e-2…1e-6 | 5 | analytic parametrization **fitted to simulations in the literature**, with radii from a named simulation |
| SN Kolmogorov turbulence, C\* = 1 %, 10 % | 2 | the construction used in the literature, with its spectrum, damping scale and amplitudes |
| Earth, non-PREM crust, 3 zenith angles | 3 | PREM is real; the three extra crustal layers are **invented** |

Sources, so the shapes are defensible rather than invented: Fogli, Lisi, Mirizzi & Montanino
(PRD 68, 033005) for the progenitor law ρ₀ = 10¹⁴ (x/km)^−2.4, the forward-shock jump ξ = 10 and
the rarefaction shape; Kneller & Kabadi (PRD 92, 013009) for the turbulence construction, the
k^−5/3 spectrum, the 100 km damping scale and the three discontinuity radii (reverse shock
1734 km, contact 12 348 km, forward shock 30 323 km) read off a 10.8 M☉ simulation at t = 3 s;
Bahcall, Serenelli & Basu (ApJ 621, L85) for the solar table.

### 13.2 The profiles were validated before anything measured with them

`validate_physical.py`, 43 checks, all passing: jump factors against the literature formulas,
spectral index, kink placement, resonance on the trajectory, and the API contract on every
family. **Three of its own checks were wrong**, and each would have produced a finding about the
construction rather than about the package:

* the contact-discontinuity jump read 0.72 instead of 2.5 at the widest front, because the check
  divided out an analytic background whose shape varies across a 700 km sampling interval;
* the turbulence rms read half of C\*, correctly — under a k^−5/3 spectrum most of the variance
  sits in the longest mode, and over a window shorter than its wavelength `np.std` removes it
  with the mean. C\* is an rms over *realizations*; measuring it over *space* was the error;
* the resonance was reported as absent from the supernova ray. It is present — but it is the
  **H** resonance, and `battery2.ne_res_for` hard-codes the 0–1 gap, which is the L resonance.
  At 15 MeV the crossing sits **at the forward shock**, the configuration the shock-effect
  literature studies.

### 13.3 The oracle was checked before anything was concluded from it

Every error here is measured against DOP853 at `rtol = 1e-12`, an eighth-order method whose error
estimate assumes a Taylor expansion exists — and the Earth chord carries PREM's steps while the
shock at w = 1e-6 is a 0.07 km ramp on a 7×10⁴ km ray. Earlier tranches sidestepped this with
`expm` composed across segments wherever the profile was piecewise *constant*; neither of these
is, so the oracle had to be checked instead. `physical_battery.py oracle_check` tightens to
`rtol = 1e-13` on the ten hardest cases: worst movement **4.9e-10** against the 1e-3 being
quoted. Six orders of margin.

### 13.4 P1 — two silent misses in 195 configurations

`warn_fp.py --physical`, all 195 scored, none skipped.

| | |
|---|---|
| outside the requested 1e-3, **and warned** | 31 |
| outside the requested 1e-3, **SILENT** | **2** |
| worst error anywhere | 2.03e-01 (SN shock w = 1e-6) — warned |

| family | worst error | outside | silent |
|---|---|---|---|
| tabulated linear/cubic, all N | **1.15e-04** | **0** | 0 |
| BS05(AGS,OP) linear | 6.79e-04 | 0 | 0 |
| BS05(AGS,OP) cubic | 1.44e-03 | 1 | **1** |
| Earth, non-PREM crust | 1.97e-03 | 2 | 0 |
| SN shock w = 1e-2 … 1e-6 | 1.55e-03 → **2.03e-01** | 23 | **1** (at w = 1e-3) |
| SN turbulence, both amplitudes | 1.39e-02 | 7 | 0 |
| every shock width, `t_breakpoints` declared, N = 8 | **5.93e-04** | **0** | **0** |

**The interpolation-kink family, which the brief expected to pay, paid nothing** — worst error
1.15e-04 across every node count from 20 to 5000, linear and cubic alike, not one configuration
outside tolerance. Ordinary tabulate-and-interpolate is not a hazard for this package.

### 13.5 Attribution — both silent misses are §12.1 again, in the half it could not reach

| | BS05(AGS,OP) cubic, d=3, 100 MeV | SN shock w=1e-3, d=3, 15 MeV |
|---|---|---|
| engine `'auto'` picked | hybrid | hybrid |
| `strategy_info['certified']` | **True** | **True** |
| error | 1.444e-03 | 1.095e-03 |
| cumulative scan, same request | **1.308e-06** | **2.586e-06** |
| ratio | **1100×** | **420×** |
| why it was not reached | single point; seam yields at N ≥ 8 | same |

§12.1 found the hybrid path answering while a far more accurate cumulative scan stood aside, and
moved `HYBRID_YIELDS_TO_CUMULATIVE_MIN_POINTS` from 25 to 8. **Eight is a threshold on the number
of points, so it repaired scans and left single points exactly as they were** — and both
remaining silent misses are single points.

It is not a limit of the method: the BS05 call at `rtol = atol = 1e-6` returns **5.911e-10**,
nine orders better, and warns. It is not a knife-edge either — sweeping energy on that profile
gives 3.96e-05, 4.61e-04, **1.44e-03**, 3.69e-04, **6.39e-03** at 70/85/100/115/130 MeV, silent
at every one. The miss is a band, and worse than the population found.

**And it reaches a physical energy — which is the answer to the brief's question.** The
population uses 30 and 100 MeV on solar profiles because `solve_ivp` is cheapest there, but no
solar neutrino has those energies, so a finding confined to them would not be a finding about
reachability. `bs05_energy_band.py` sweeps 5–100 MeV over both interpolants and both dimensions,
28 configurations, verifying the oracle at every one:

| | silent miss | in the ⁸B band? |
|---|---|---|
| cubic, d = 3, 100 MeV | 1.444e-03 | no |
| **cubic, d = 2, 5 MeV** | **1.380e-03** | **yes** |
| **linear, d = 2, 5 MeV** | **1.178e-03** | **yes** |

**Two of the three silent misses are at 5 MeV**, in the middle of the ⁸B spectrum, in the
two-flavour treatment that is the standard approximation for the solar problem — and they occur
with *both* interpolants, so this is not an artefact of choosing a spline. All 28 rows report
`certified: True`, and **not one of the 28 raises any warning at all**.

So a user computing solar-neutrino oscillations at 5 MeV from the BS2005-AGS,OP table, by either
of the two obvious interpolation choices, receives an answer ~1.2–1.4× outside the tolerance they
asked for, certified, in silence. That is a real user, a real model and a real energy.

*(The script's own first verdict said the opposite — "does not reach the solar band" — because it
reported the single largest silent miss and then asked whether **that one** was in the band. The
largest is at 100 MeV; the question is whether **any** is. The logic is fixed and the comment in
`bs05_energy_band.py` records why.)*

### 13.6 Two claims corrected against fuller data

Both were stated from partial results during the run and are wrong.

* **"`t_breakpoints` cures the shock" holds only for scans.** At N = 8 it does
  (2.03e-01 → 8.75e-06, zero outside tolerance at every width). At **N = 1 it makes the answer
  9.5× worse**: 1.095e-03 → 1.043e-02, because declaring the fronts moves dispatch off the hybrid
  path onto the general Magnus ladder. It does become loud, which is the better failure mode, but
  it is not a cure.
* **"The BS05 defect is the cubic spline" is wrong.** The linear interpolation of the same model
  is silently outside tolerance too — 2.906e-03 at 130 MeV — and the cubic only moves the onset
  earlier. Both report `certified: True`. The defect is the hybrid path's self-certification on a
  real tabulated model, not the choice of interpolant.

### 13.7 P2/P3 — the hidden-feature scan, and a category error in the brief

Ground truth here is **measured, not declared**: total variation on a grid 32× finer than the
6400-point reference grid, over that on the reference grid. A flag set by whoever built the
profile is how a construction error becomes a finding, and the first version of this made exactly
that error.

| | |
|---|---|
| **P2** false positives on families with no sub-grid structure | **0 of 17** |
| **P3** detection on families that genuinely hide something | **0 of 2** |

**The brief expected the shock family to exercise this scan; it cannot.**
`find_hidden_features` computes `per_interval − endpoints` — variation the fine grid sees inside
a reference interval, minus what the two bracketing reference nodes already show. A shock front
is a **monotone step**, so those two nodes see its full height however narrow it is, and the
difference is zero by construction. Measured: fine/coarse variation ratio **1.000 at every width
from 1e-2 to 1e-6**. What a narrow shock defeats is *resolution*, not visibility.

The turbulence is the real case — ratio 1.59 at C\* = 1 % and **3.19 at 10 %**, so the finest grid
the package reaches sees under a third of the profile's variation — and the scan misses it at
both amplitudes, at concentrations of 0.0015 and 0.0026 against a 0.30 threshold. **This is a
limit of the statistic, not a mis-set threshold**: concentration asks whether hidden variation
piles up in one interval, and a power-law spectrum spreads it evenly over all 6400. No threshold
reaches broadband roughness, just as none reaches a feature that was never sampled (§12.3).

### 13.8 The resolution test on the physical population

`resolution_fp.py --physical`, 2736 configurations:

| family | declared unresolved |
|---|---|
| tabulated linear, N = 20…5000, and cubic | **0 of 144** at every N |
| BS05(AGS,OP), linear and cubic | **0 of 144** |
| Earth, non-PREM crust | **144 of 144** at every zenith angle |
| SN shock w = 1e-2, 1e-3 | 0 of 144 |
| SN shock w = 1e-4 / 1e-5 / 1e-6 | 48 / 60 / **144** of 144 |
| SN turbulence, both amplitudes | 0 of 144 |

Earth's PREM steps are flagged every time, correctly — that is what the `t_breakpoints`
`osc_prob_earth` supplies exist to declare. But note the inversion that constrains any future
fix: **the widths that get flagged are the ones already loud, and the width that fails silently
is never flagged.** w = 1e-3 is 0 of 144 and is the silent miss; w = 1e-6 is 144 of 144 and warns
every time. A 70 km front on a 7×10⁴ km ray genuinely *is* resolved by the 6400-point grid — what
fails there is the certification, not the sampling. A cure gated on that flag could not fire on
the case that needs it.

### 13.9 Cross-check acceptance on the physical population

`crosscheck_acceptance.py --physical`, 38 configurations: **0 X1 failures**. Whenever one engine
was outside tolerance and another was not, the cross-family spread saw it — 15 of 15. Two rows
are X2, where **every** engine is wrong and there is nothing to disagree about: the shock at
w = 1e-5 and 1e-6 at d = 3, where hybrid is wrong by 0.854 and even the cumulative scan by
6.9e-03. Every engine warns on both.

### 13.10 Two documentation defects, found by building a real solar profile

Neither changes a number; both would mislead a reader. Docstrings only —
`NUM_DENSITY_E_SUN_CENTRAL` and `L_SCALE_SUN` are bit-for-bit unchanged.

* **`L_SCALE_SUN` was documented as "Electron number density at the center of the Sun"**, a
  copy-paste of the constant above it. It is a length.
* **`NUM_DENSITY_E_SUN_CENTRAL = 245 N_A` was documented as the Sun's central electron density.**
  It is the r → 0 intercept of the standard exponential *fit*. BS2005-AGS,OP gives 102.7 N_A at
  its innermost point. Measured against that table, the fit is high by 2.4× inside 0.05 R☉,
  agrees to 2.5 % only around 0.2–0.3 R☉, and departs by up to 89 % beyond 0.7 R☉.

The consequence is for this document rather than for users: **the trajectory every battery in
this directory uses — one solar scale height, 0.095 R☉ — lies where that fit is high by about
30 %.** The "solar" population is a legitimate exponential profile; it is not a solar model, and
no earlier tranche said so.

### 13.11 Fallback quality on the physical population — the single-point exposure, independently

`fallback_quality.py --physical`, 100 workloads (19 families × d = 2, 3 × the request shapes that
route differently). This measures something P1 does not: not "is the answer wrong" but "was a
better engine standing right there".

| | synthetic, before §12.1 | synthetic, after | **physical** |
|---|---|---|---|
| workloads | 42 | 42 | **100** |
| `'auto'` more than 10× worse than the best engine that applied | 30 | 17 | **31** |
| worst factor | 900 000× | — | **595 090×** |
| `'auto'` outside tolerance | 3 | 1 | **15** |
| of those, silent | 2 | 0 | **1** |
| cumulative scan applied | 28 of 42 | — | **76 of 100** |

**Of the 31 workloads where `'auto'` is more than 10× worse than an engine that applied, 27 are
single points and 4 are scans.** That is the §13.5 attribution arriving independently, from a
measurement that knows nothing about warnings: the dispatch-order penalty on physical profiles is
concentrated almost entirely where the seam cannot reach.

The worst case is a supernova shock at w = 1e-2, d = 3, single point: `'auto'` returns 1.55e-03
via the **general Magnus ladder** — the hybrid path declined — while the cumulative scan on the
same request returns 2.61e-09.

Two honest qualifications:

* **The cumulative scan is not universally better.** It answered 76 of 100 workloads with a worst
  error of 5.10e-03, and on one workload (`SN shock w=1e-3 d=2 L-scan N=8`) `'auto'` picked it and
  the hybrid path would have been 15× better. Any change here is a change of *default*, not the
  discovery of a dominant engine.
* **Whether reaching it at N = 1 is affordable is exactly what P4 measures**, and P4 has not been
  run. Nothing in this section is an argument for moving the seam until that number exists.

### 13.12 Pass criteria, stated before running and answered after

| | criterion | result |
|---|---|---|
| **P1** | no silent miss on any physical family | **FAILED — 2 of 195**, both attributed to §13.5, one reaching 5 MeV on a real solar model |
| **P2** | `find_hidden_features` stays at 0 false positives on families with nothing to hide | **PASSED — 0 of 17** |
| **P3** | detection rate on families that genuinely hide something | **0 of 2** — structural (§13.7), not a threshold |
| **P4** | the 25 → 8 seam still holds on physical profiles | **NOT RUN** — see `adversarial_batteries/RUN_P4.md` |

P4 is the one deliverable of this tranche that has not been measured. Its first round returned a
control drift of 0.99×, so the method is sound and only the run is outstanding; the command and
the reasoning are in `RUN_P4.md`. **Nothing in §13.11 should be read as an argument for moving
the seam until that number exists** — the case for reaching the cumulative scan at N = 1 rests
entirely on what it costs there, and that is precisely what has not been measured.

### 13.13 What this tranche establishes, in one paragraph

**The exposures do reach a real user, and the brief's prediction was wrong in both directions.**
The interpolation-kink family it expected to pay produced nothing — worst error 1.15e-04 across
every node count from 20 to 5000. The families that pay are a **real published solar model**,
where `strategy='auto'` returns a certified answer 1.2–1.4× outside the requested tolerance in
silence at 5 MeV with either interpolant, and a **supernova envelope**, where the error reaches
0.203 but the package always says so. Everything catastrophic is loud; what is silent is
marginal. Both silent misses are single points where a 400–1100× more accurate engine was
available and structurally unreachable, which makes them the same defect §12.1 found and
half-fixed — eight is a threshold on point count, and a single point cannot cross it.

### 13.14 Still open after this tranche

* **P4** — the one unmeasured pass criterion (`RUN_P4.md`).
* **The single-point half of the dispatch-order defect.** 27 of the 31 workloads where `'auto'`
  is more than 10× worse than an available engine are single points. Whether the seam should
  reach them depends on P4's N = 1 row.
* **Nothing detects broadband sub-grid roughness.** The turbulence family defeats both
  instruments (§13.7, §13.8). A statistic that *would* see it is nearly free —
  `find_hidden_features` already computes total variation on both the fine and the reference
  grid, and their **ratio** is 1.000–1.003 on all 17 families with nothing to hide against 1.59
  and 3.19 on the two turbulent ones — and it discards that ratio. Not built here: the turbulence
  errors are all *caught* by the convergence machinery already, so a new detector would need its
  own false-positive measurement before it could earn its place (§5c).
* **Auto-inserting `t_breakpoints` remains rejected, now with a sharper reason.** Beyond §12.2's
  "it changes dispatch": at N = 1 declaring the shock fronts makes the answer **9.5× worse**
  (§13.6), and the resolution test flags the widths that are already loud while missing the one
  that is silent (§13.8) — so a cure gated on that flag could not fire on the case that needs it.
* **The reverse shock was not carried.** The ray starts at 1e4 km, outside it, because the oracle
  cost diverges inwards as r^−1.4. Two of the three literature discontinuities are on the path.
* **The Earth crust layers are invented.** PREM is real; a genuine 3-D tomographic crust model
  would be better evidence and was not available.

### 13.15 The fix for the two silent misses — BUILT, MEASURED, REVERTED

P4 said the cumulative scan is cheaper than the hybrid path at N = 1 on every physical profile
(median 0.15×), and §13.5 said it is 400–1100× more accurate on exactly the two configurations
that fail silently. That looks like a free win, and it was implemented:
`HYBRID_YIELDS_TO_CUMULATIVE_MIN_POINTS` 8 → 1 **and** `CUMULATIVE_AUTO_MIN_POINTS` 2 → 1. Both
are required — 8 is a threshold on *point count*, and a single point cannot cross one however far
it falls.

**It works.** A: 1.380e-03 → **1.707e-05**. B: 1.095e-03 → **2.586e-06**.

**It is reverted anyway**, for three reasons, the first of which invalidates the measurement that
motivated it:

1. **Standing aside does not hand the request to the cumulative scan.** Dispatch runs hybrid →
   interaction picture → separable → cumulative, so when the hybrid path declines, the next
   engine that *applies* answers — and on several workloads that is `ip_exp`. P4 compares hybrid
   against cumulative; it never measured hybrid against *whatever comes next*, which is the
   thing that actually happens. Two tests caught it as `assert 'ip_exp' == 'hybrid'`.
2. **The cumulative branch has never carried default traffic**, and making it do so surfaced two
   latent defects in minutes (§13.16). That is evidence about the branch's coverage, not about
   those two bugs.
3. **It is a change of default, not a dominant engine** — the cumulative scan's worst error over
   the 76 physical workloads it serves is 5.10e-03, and on one the hybrid path is 15× better.

The right fix is **surgical**: have the hybrid path yield *to the cumulative scan specifically*,
declining only when that engine will actually take the request, rather than declining into the
chain. `_cumulative_scan_would_serve` already exists and is the natural place to build on.
Recorded in `HANDOVER_SINGLE_POINT_EXPOSURE.md`; both constants keep the full measurement in
their docstrings, in the same built-measured-reverted form as §11.3 and §12.2.

### 13.16 Two shipped bugs, found because the change exposed them

Both were unreachable by default only because `CUMULATIVE_AUTO_MIN_POINTS = 2`, and both are
reachable **on the shipped tree** through an explicit `cumulative=True`. No test covered the
combination. Both fixed here, and kept after the revert.

* **Missing scalar squeeze.** `osc_prob_3nu_earth(E, ..., cumulative=True)` with scalar energy
  and baseline returned `(1, d, d)` instead of `(d, d)` — so `P[nu_i][nu_f]` silently selected a
  row of the wrong array. The cumulative branch was the one return site in
  `osc_prob_energy_baseline` that omitted `__getitem__(0 if return_float else slice(None))`.
* **`convergence_info` forwarded to an engine that rejects it**, raising `TypeError` instead of
  returning a probability. Now dropped, for the same reason `strict_convergence` already was.

Neither has a regression test yet; that is the smallest open item in this document.

### 13.17 CORRECTION — P1 measured the wrong observable for the solar family

Every number in §§13.4–13.5, including the headline, is `max|P_package - P_oracle|` **at a single
baseline**. For solar neutrinos that is not the observable. The trajectory is **398 oscillation
lengths** at 5 MeV; the ⁸B production region is extended, the Sun–Earth vacuum phase is ~1e10
cycles, and detector energy resolution washes the remainder. The quantity a user sees is the
**phase-averaged** survival probability — which is why the package ships `magnus.avgprob`, whose
`averaged_probabilities_adiabatic` is the standard solar MSW path and **does not route through
the instantaneous propagation at all**.

Measured (`adversarial_batteries/avg_check.py`), averaging over a six-oscillation-length window
at the nominal baseline:

| BS05, d = 2, 5 MeV | instantaneous | **averaged** | reduction |
|---|---|---|---|
| cubic | 1.380e-03 | **2.603e-05** | 53× |
| linear | 1.178e-03 | **1.595e-04** | 7× |

Averaged `P_ee`: package 0.593732, oracle 0.593758 — five significant figures.

**So the headline of §13.4 is overstated and is corrected here.** The claim that survives is
narrower:

* For the **averaged** observable — solar physics as normally done — the package is **inside
  tolerance** at 5 MeV, by 6× to 38×. There is no silent miss for that user.
* For the **instantaneous** probability at a single baseline, it is outside tolerance and silent.
  That is a real defect for anyone who needs the coherent probability at a point, and the worst
  instantaneous error over the window is **3.4e-03**, worse than the 1.38e-03 at the nominal
  baseline — but it is a narrower class of user than §13.4 implied.
* The averaging is **not** a complete rescue. 1.595e-04 for the linear interpolant is only 6×
  inside the default tolerance and is still silent, so a caller asking for `rtol = 1e-5` on the
  averaged observable would be outside it with no warning.

**The methodological lesson is the larger one.** The brief's P1 — "`strategy='auto'` at the
default tolerance is inside 1e-3 **or** warns" — was executed exactly as written, and it is the
right criterion for a general-purpose propagator. But for families whose physical observable is
an average, it measures a quantity the user never sees, and it will therefore report failures
that are not failures. Any future tranche that includes solar, turbulence, or any other
fast-oscillation family should state, per family, **which observable the criterion is applied
to** before measuring anything. That question was not asked at the start of this one.

Unaffected by this correction: the two shipped bugs and their fixes (§13.16), the P4 cost
measurement (§13.12), the structural findings about the two detectors (§13.7, §13.8), and the
supernova shock results — where the loud errors reach 0.203 and are not phase artefacts.

### 13.18 The same test on the supernova families — and it splits the finding in two

§13.17 withdrew the solar headline as a phase artefact. The supernova numbers needed the same
test (`adversarial_batteries/avg_check2.py`, six-oscillation-length window, d = 3 so the window
is set by Δm²₃₁; the ray is 4729 oscillation lengths long):

| case | instantaneous | **averaged** | reduction | warns |
|---|---|---|---|---|
| BS05 solar, d = 2, 5 MeV | 1.380e-03 | **2.603e-05** | **53×** | no |
| SN shock w = 1e-3, d = 3, 15 MeV | 1.095e-03 | **9.773e-04** | 2× | **no** |
| SN shock w = 1e-6, d = 3, 15 MeV | 2.033e-01 (worst 6.770e-01) | **2.135e-01** | 3× | yes |

**The shock error does not average away, and it should not.** A shock front changes the
*adiabaticity* of the MSW level crossing, so it changes the conversion probability itself rather
than the phase of an oscillation — which is exactly the effect the shock-effect literature exists
to study. The averaged error at w = 1e-6 (2.135e-01) is in fact **larger** than the instantaneous
error at the nominal baseline (2.033e-01).

So the two exposures this tranche reported as one defect are not equivalent:

* **Solar (A) — withdrawn.** The observable is accurate to 2.6e-05. Not user-facing.
* **Supernova (B) — stands, and is the real result.** A 0.21 error in the averaged, observable
  conversion probability, on a profile taken from the literature. The package **warns**, which is
  correct behaviour, and declaring `t_breakpoints` cures it by three to five orders on scans.
* **The w = 1e-3 silent miss is marginal on the observable**: 9.773e-04 averaged, inside 1e-3 by
  2 %. Still silent, still worth closing, but it is not the headline it was written up as.

**What this changes about the proposed dispatch fix.** Its principal justification was case A,
and case A is not a user-facing error. The remaining justification is the marginal w = 1e-3 case
and the general point that a better engine was reachable and cheaper. That is a much weaker
mandate than §13.15 assumed, and it lowers the priority of that work below the two shipped bug
fixes of §13.16.

**And the methodological point, now with a measurement behind it.** Whether an error is
observable depends on the *family*, and cannot be decided once for a population: the identical
test collapses a solar error by 53× and leaves a supernova error untouched. Any future population
must state, per family, which observable its pass criterion applies to. The check itself costs
about three minutes per case and should be built in from the start, not bolted on at the end.

### 13.19 All three fast-oscillation families, and a diagnostic that falls out of it

Completing §13.18 with the turbulence family (25 points, 1587 s — it is by far the most
expensive, and `'auto'` routes it to the general Magnus ladder):

| family | instantaneous | **averaged** | reduction | warns | verdict |
|---|---|---|---|---|---|
| BS05 solar, d = 2, 5 MeV | 1.380e-03 | **2.603e-05** | **53×** | no | phase — not user-facing |
| SN turbulence C\* = 0.1, d = 3 | 1.701e-03 | **1.565e-04** | **23×** | yes | phase — not user-facing |
| SN shock w = 1e-3, d = 3 | 1.095e-03 | **9.773e-04** | 2× | **no** | marginal, and silent |
| SN shock w = 1e-6, d = 3 | 2.033e-01 | **2.135e-01** | 3× | yes | **real, and reported** |

**The reduction factor is itself a diagnostic**, and it separates the population cleanly with no
overlap: errors that shrink by 20× or more under averaging are *phase*, and phase is not
observable; errors that shrink by 2–3× are *envelope*, and the envelope is exactly what the
detector measures. The physics agrees — smooth or oscillatory structure perturbs the phase of a
fast oscillation, while a shock front changes the **adiabaticity of the level crossing**, which
moves the conversion probability itself. That is why the shock is the only family whose error
survives, and it is the effect the shock literature exists to study.

**So the reachability question, finally answered.** Of everything this tranche measured, exactly
one class of error reaches a user of the physical observable: a **supernova shock front sharp
enough to be unresolved** — and on those the package **warns every time**, flags them 144/144 in
the resolution test, and is cured to 8.75e-06 by declaring `t_breakpoints`. The single silent
case sits at 9.773e-04 on the observable, inside the requested 1e-3 by two per cent.

**The residual risk, and it is a real one.** That silent case is w = 1e-3 — a 70 km front, which
is what a shock looks like **after a simulation has smeared it over a few grid cells**. A real
hydrodynamic shock is mean-free-path thin (w ≈ 1e-6, loud); nobody hands this package a real
shock, they hand it a snapshot. So the silent band is plausibly the *most likely* form a user's
shock actually arrives in, and 2 % of margin is luck rather than headroom.

The mechanism is a hypothesis with two pieces of evidence rather than a measured fact: at 70 km
the front spans about six cells of the 6400-point **probe** grid, so the resolution test calls it
resolved (0 of 144), while the far coarser **transport** grid still straddles it. If that is
right, the detector is checking a grid the answer does not depend on. **Mapping that band —
width × energy × tolerance — is the one open question about whether this package is silently
wrong on a plausible real input**, and is the first item of `HANDOVER_SINGLE_POINT_EXPOSURE.md`.

### 13.20 The silent band does not exist — and `t_breakpoints` is not the cure it was called

`adversarial_batteries/shock_silent_band.py`, 9 shock widths (2100 km down to 7 km) × 2 energies,
each on the **averaged** observable (41 points over six Δm²₃₁ oscillation lengths), with the
resolution test's verdict and the `t_breakpoints` result recorded alongside. 18 configurations,
about 70 minutes.

| w | front | E | averaged err | res200 | res6400 | verdict |
|---|---|---|---|---|---|---|
| 3e-2 … 5e-3 | 2100–350 km | both | 6.9e-05 … 1.8e-04 | True | True | inside, warned |
| 3e-3 | 210 km | both | 1.2–1.3e-04 | **False** | True | inside, warned |
| 2e-3 | 140 km | 15 | 1.611e-04 | True | True | **inside, quiet** |
| **1e-3** | **70 km** | **15** | **9.840e-04** | True | True | **inside, quiet — 98 % of tolerance** |
| 5e-4 | 35 km | 15 | 6.048e-04 | True | True | inside, warned |
| 3e-4 | 21 km | 15 | **3.279e-03** | **False** | True | outside, **warned** |
| 1e-4 | 7 km | 15 | **3.933e-03** | **False** | True | outside, **warned** |

**There is no silent band.** Two configurations are outside tolerance on the observable and
**both warn**; only two configurations are quiet at all, and both are inside. The w = 1e-3 spike
is real — a 6× excursion above its neighbours, reaching **98 % of the requested tolerance** — but
it never crosses. The package is honest across the entire sweep.

**How it stays honest is worth recording, because it is not by design.** The 6400-point probe
grid declares the profile *resolved* on the two rows that are outside tolerance; what catches
them is the **200-point** grid saying *unresolved*. The coarse probe rescues the fine one. That
inverts the naive expectation and confirms the §13.19 hypothesis in part — the fine grid's
verdict genuinely is uninformative about whether the transport can integrate the front — but the
failure it would cause is masked, not present, because the detector chain asks two grids and
takes the pessimistic answer.

**And a correction to a claim made repeatedly in this document.** §13.4 and §13.19 say
`t_breakpoints` cures the shock, on the strength of the N = 8 rows in `warn_fp.py --physical`
(2.03e-01 → 8.75e-06). On the **averaged observable at single points** it does not:

| | of 18 configurations |
|---|---|
| `t_breakpoints` improved the answer | **7** |
| it made the answer **worse** | **11** |
| it pushed an answer that was **inside** tolerance **outside** it | **2** |

`w = 5e-4` at 15 MeV goes 6.048e-04 → 1.275e-03, and `w = 1e-4` at 45 MeV goes 2.858e-04 →
1.682e-03. This is the same mechanism as §13.6: declaring breakpoints at a single point moves
dispatch off the hybrid path onto the general Magnus ladder, which is often worse there. So the
honest statement is **`t_breakpoints` cures the shock on baseline scans and is a coin flip on
single points** — and the advice `UnmarkedDiscontinuityWarning` gives is right for the first case
and not established for the second.

### 13.21 An aliasing warning — BUILT AS A MEASUREMENT, REJECTED, REPLACED BY A STATISTIC

Proposed after §13.17–13.19: warn when a scan samples an oscillation too coarsely to represent
it, and point the caller at `magnus.avgprob`. The criterion is Nyquist and is objectively
correct — above half a wavelength per step the returned array cannot represent the oscillation,
which is wrong regardless of intent, so it needs no guess about what the user wanted.

**Cost measured before writing any of it** (`alias_cost.py`, alternating with a control, minima
read because interference only adds time): an 8-point probe costs **1.8 % of the cheapest scan**
and under 0.1 % of a substantial one; a 64-point probe costs 11 % and is out.

**Then the false-positive rate killed it** (`alias_fp.py`):

| | |
|---|---|
| realistic scan sizes the criterion would fire on | **44 of 45 = 98 %** |
| 8-point probe disagreeing with a 4096-point reference | **0 of 45** |
| baselines needed for Nyquist: Earth / solar / SN ray | 861 / 4 390 / **73 392** |

The probe was never the problem — it is accurate and cheap. **The concept fails**: a warning
firing on 98 % of calls is noise however right each firing is, and it would train users to
silence the category that also carries `UnmarkedDiscontinuityWarning`.

**Shipped instead**: `adiabatic.oscillation_sampling`, surfaced as `strategy_info['sampling']`
with `oscillation_length`, `cycles_over_trajectory`, `spacing`, `cycles_per_step`,
`nyquist_points` and `aliased`. Computed **only when the caller passes `strategy_info`**, so the
default path costs nothing — pinned by a test that spies on the helper and asserts it is not
invoked otherwise. Opt-in cost measured at 5.5 % of the cheapest scan.

**One defect found in my own prototype, and it is the interesting part.** The first version took
`min(diff(sort(λ)))` — the smallest *adjacent* gap. That is the **slowest** oscillation. Aliasing
is set by the **fastest**, which comes from the largest eigenvalue spread. The wrong statistic
looked healthy: it agreed with a 4096-point reference to 1.000× on every family, which I nearly
reported as evidence that a cheap probe suffices. It was evidence only that I was measuring
something smooth. With the right statistic the agreement is 1.022× worst case — genuinely
adequate, and for a reason: the largest spread tracks the matter potential and varies smoothly,
while the smallest gap has a sharp minimum at an MSW resonance and would need dense sampling.
A test now pins the distinction specifically.

Documented in `docs/source/averaged_probability.rst`, together with the phase-versus-envelope
table of §13.19 — which is the part a user needs, because it says when a large instantaneous
error is harmless and when a small one is not.
