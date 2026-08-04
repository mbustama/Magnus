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
