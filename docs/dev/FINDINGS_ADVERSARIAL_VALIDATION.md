# Findings: adversarial validation of the dispatch, cumulative and adiabatic changes

**Written:** 2026-08-03, executing `HANDOVER_ADVERSARIAL_VALIDATION.md` against
`notebook-breakpoints-and-cumulative` (eleven commits ahead of `main`).
**Machine:** the same 12-core box. **Nothing was pushed, and nothing was fixed** — the brief
asks for documented, reproduced failures rather than silent patches, and that is what this is.
**Baseline reproduced:** 672 passed; `ruff check src/ tests/` clean (the 63 findings under
`ruff check .` are all in `notebooks/` and are byte-for-byte identical on `main`).

---

## 1. The verdict, in one paragraph

**The eleven commits are sound. Ship them.** Across seven batteries and roughly 600 scored
configurations I could not construct a case where the branch is worse than `main` — not one.
The nine bit-identity rows are bit-identical; the two rows that moved moved toward the truth by
175 000× and 1 138×; 4ν/5ν accuracy evidence, which did not previously exist, is now supplied
and is clean.

**But the batteries did find defects — all of them pre-existing in `magnus.adiabatic`, none
introduced here, and none fixed by commit 9c7945a either.** Every failing configuration I found
is **bit-identical between `main` and the branch**, which is the strongest statement available
in both directions: the branch introduces no regression, and the γ-sweep fix is narrower than
its commit message implies. The headline is that the default `strategy='auto'` can return a
probability wrong by **0.54** with no warning at all.

The one number that carries the whole report: on 243 random profiles, **51 (21 %) came back
outside the requested tolerance with no warning — and every single one of them was at N < 25,
on the hybrid path. At N = 30 and N = 80, where this branch routes scans to the cumulative
scan, there were none, and the median error was five orders of magnitude better.**

---

## 2. What was run

| battery | subject | configurations | failures |
|---|---|---|---|
| 1 | Full regression, bit-identity vs `main` | 11 | 0 |
| 2 | The detector's fixed probe grid | 42 | **12** |
| 3 | Routing seams | 134 | **3 classes** |
| 4 | Extreme numerics | 28 | 0 |
| 5 | Designed to break | 32 + 250 fuzz | **3 classes** |
| 6 | Flavor count as a first-class axis | 72 | 0 |
| 7 | Cross-module and oracle diversity | 17 | 0 |
| — | Package unitarity, d = 2…5 × N = 25…10⁵ | 16 | 0 |

**Oracle discipline.** `solve_ivp`/DOP853 at `rtol=1e-12, atol=1e-14` is the only accuracy
oracle, and its convergence was verified per-configuration by tightening to `rtol=1e-13` and
confirming the movement is ≪ the error quoted. Where a case was adversarial enough that an
*adaptive* oracle could itself have stepped over the feature, three further independent oracles
were used (§3.1). No Magnus path was ever scored against another.

---

## 3. Defects found, ranked

### 3.1 A step-function profile with an unmarked edge returns P wrong by 0.54, silently

**Severity: highest.** Default settings, ordinary entry point, no warning.

```python
# electron density: 0.02*Ne0 for l < mid, 0.30*Ne0 for l >= mid; edge NOT declared
P = oscprob.osc_prob_matter_std_potential(
        2, step_profile, 50e6, L1, osc_params, L0=0.0,
        density_is_of_number_of_electrons=True)      # strategy='auto' is the default
```

| | P_ee | error | warnings |
|---|---|---|---|
| **`strategy='auto'` (default)** | **0.589270087** | **5.395e-01** | **none** |
| `strategy='hybrid'` | 0.589270087 | 5.395e-01 | none |
| `strategy='magnus'` | 0.049392132 | 3.610e-04 | `MagnusConvergenceWarning` |
| truth (`expm`, **exact** here) | 0.049753176 | — | — |

The Hamiltonian is constant on each half, so `expm` composed across the two pieces is not an
approximation but the exact answer. A naive `solve_ivp`, a piecewise `solve_ivp` restarted at
the discontinuity, and a Magnus run with the edge declared as a breakpoint all agree with it to
**1.3e-11**. The hybrid answer is wrong by more than half a probability unit.

**Mechanism.** `magnus.adiabatic`'s module docstring already states the method is "restricted to
a smooth `VCC_func` … a piecewise-discontinuous profile such as PREM breaks the finite-difference
diagnostics this method relies on". The intent is right; the *implementation* detects
discontinuity only through the proxy `did the caller pass t_breakpoints or t_slab_edges`
([oscprob.py:3587](src/magnus/oscprob.py:3587), [oscprob.py:3591](src/magnus/oscprob.py:3591)).
That proxy **fails open**: it declines when the user *tells* it about the discontinuity, and
accepts when they do not. The guard is exactly backwards with respect to the risk.

Two milder members of the same family, also silent under the default:

| profile | `auto` | `magnus` |
|---|---|---|
| kink, no jump (C⁰ but not C¹) | **1.448e-02, silent** | 6.188e-06, warns |
| singularity approached, not reached | **8.625e-03, silent** | 2.144e-04, warns |
| sawtooth (jump every 0.07 span) | 7.484e-03, warns | 7.484e-03, warns |

Note the contrast with `NOTES_ADAPTIVE_REFINEMENT.md` §4b, which classes an unmarked
discontinuity as a caller error that "no refinement strategy can help". That is true of the
*Magnus ladder*, and it was the right diagnosis there. It is not true here: the library's own
general path gets this profile right to 3.6e-04 and says so. Only the hybrid path is wrong, and
only the hybrid path is silent.

### 3.2 Sub-threshold γ accumulation: hybrid certifies while wrong, up to 3.9e-02

**This is the mechanism commit 9c7945a addressed, and the fix does not reach this case.**

The γ sweep added in 9c7945a fixed the case where γ is large *along the path* but small *at the
gap extrema*. It does nothing when γ never exceeds the **initial** threshold `threshold0 = 0.1`
anywhere: no window opens, successive refinements differ only in the adiabatic-transport grid,
they converge to the same wrong adiabatic limit, agree with each other, and `hybrid_propagator`
certifies. The threshold ladder (`threshold /= 3` down to `min_threshold = 1e-6`) never gets a
chance to help, because the agreement test fires first.

A deliberate construction sweeping a Gaussian resonance width, so γ_max crosses 0.1 (Battery 5.1):

| width / span | γ_max | hybrid err | windows | `magnus` err |
|---|---|---|---|---|
| 3e-01 | 4.312e-03 | **2.222e-03** | 0 | 3.434e-04 |
| 1e-01 | 1.294e-02 | **8.010e-03** | 0 | 7.919e-05 |
| 6e-02 | 2.156e-02 | **1.096e-02** | 0 | 1.062e-04 |
| 4e-02 | 3.234e-02 | **1.767e-02** | 0 | 1.738e-04 |
| 3e-02 | 4.312e-02 | 8.798e-12 | 1 | 1.247e-04 |
| 2e-02 | 6.468e-02 | 1.994e-12 | 1 | 1.082e-04 |
| 1e-02 | 8.625e-02 | 5.737e-12 | 1 | 2.454e-04 |
| 7e-03 | 1.848e-01 | 3.291e-12 | 1 | 2.123e-04 |

The transition is exactly at the threshold, and it is a cliff: **γ_max ≥ 0.043 → 1e-12;
γ_max ≤ 0.032 → up to 1.8e-02, certified, silent.** Every one of these is answered correctly to
~1e-04 by `strategy='magnus'`, so this is not a hard integration problem — it is a threshold
that is consulted once, too early.

Worst case found anywhere: **3.907e-02** (Battery 2.4, ten crossings), i.e. 39× the requested
1e-3, certified. A plain sinusoidal density at a **well-resolved** period (span/7, ~28 probe
points per period) gives 1.672e-02 — the failure has nothing to do with under-sampling.

**Sobering calibration:** the multi-resonance profile this branch was developed against sits at
γ_max = 7.7e-02 (2ν) and 7.9e-02 (3ν) — within 25 % of the threshold, on the correct side.

### 3.3 The narrow-feature blind spot the brief ranked first — confirmed, but second

`n_probe = 200` linear samples cannot see structure narrower than the spacing. Sweeping one
Gaussian resonance from 1e-1 down to 1e-5 of the domain, with γ scanned on a grid 500× denser
than the detector's:

| width / span | hybrid err | γ_max | dense pts with γ>0.1 | covered by a window? | mechanism |
|---|---|---|---|---|---|
| 3e-02 | 4.388e-03 | 8.6e-03 | 0 | n/a | sub-threshold (§3.2) |
| 1e-02 | 7.701e-03 | 2.6e-02 | 0 | n/a | sub-threshold (§3.2) |
| 3e-03 | 1.969e-11 | 8.6e-02 | 0 | n/a (a window opened anyway, from the extremum) | ok |
| 1e-03 … 1e-04 | ~1e-11 | 0.26 … 2.6 | 6054 … 10357 | YES | ok |
| **3e-05** | **2.907e-02** | 8.62 | 11742 | **NO (11742/11742)** | **detection miss** |
| **1e-05** | **9.318e-03** | 25.8 | 12846 | **NO (12846/12846)** | **detection miss** |

γ exceeds the threshold by **two to three orders of magnitude** over ~12 000 dense-grid points,
and **zero** windows are reported. Reproduced identically at d = 4 and d = 5, so it is a
property of the probe grid, not of the flavour structure.

This one is worse than §3.2 in one respect: it is wrong on **every** path, including the
cumulative scan (2.952e-02 over a 60-point scan) and `strategy='magnus'` (2.907e-02, which at
least warns). The `strict_convergence` probe cannot see the feature either — it is the
`NOTES_ADAPTIVE_REFINEMENT.md` §2 blind spot ("an integral is blind to structure that averages
out"), and as the brief's Battery 5.2 predicted, it now misplaces a whole scan rather than one
point. The documented cure — caller-supplied `t_breakpoints` — works, but the caller has to know.

**Counter-intuitive result worth recording: aliasing is not where the damage is.** Sinusoids at
exactly the probe spacing, at ½ and 2× it, and at 1.01/0.99× it (beat frequencies) all come back
accurate to ~1e-11. Clustering 50 crossings into 1 % of the range: 1.6e-07. Two hundred and four
hundred crossings: 5e-09 and 9.7e-09, and the machinery honestly reports `certified=False`.
The probe is far more robust to aliasing than expected; it is broad-and-smooth, and
extremely-narrow, that defeat it.

### 3.3b The fuzzer quantifies all of the above — and draws the line exactly at the seam

250 random smooth profiles (random Fourier sums, controlled bandwidth), with the **flavour count
drawn randomly from {2,3,4,5}** with random sterile mixings, random energy (5–200 MeV), random
span, random ν/ν̄, and random scan size N ∈ {1, 3, 12, 30, 80}; each scored against `solve_ivp`.
243 scored (7 raised `ValueError: rho_func must be non-negative`, which is correct validation of
a profile my generator drove negative, not a defect).

| | |
|---|---|
| error distribution | median **6.27e-05**, p90 5.19e-03, p99 3.17e-02, max **4.00e-02** |
| outside 1e-3, warned | 10 |
| **outside 1e-3, SILENT** | **51 (21 %)** |

Broken down by scan size, this is the cleanest result in the whole exercise:

| N | cases | silent misses | median error | path |
|---|---|---|---|---|
| 1 | 43 | **16 (37.2 %)** | 4.64e-04 | hybrid |
| 3 | 52 | **14 (26.9 %)** | 2.60e-04 | hybrid |
| 12 | 59 | **21 (35.6 %)** | 6.96e-04 | hybrid |
| 30 | 43 | **0** | **4.28e-09** | cumulative |
| 80 | 46 | **0** | **7.75e-09** | cumulative |

**Every one of the 51 silent misses is below the N = 25 seam. Above it there are none, and the
median error improves by five orders of magnitude.** On random profiles, the hybrid path is
silently outside its requested tolerance about a third of the time; the path this branch routes
scans to is not, on any of the 89 cases that reached it.

Stronger than "nothing warned-and-missed above the seam": across those 89 cases the **worst**
error is **2.57e-06**, two orders inside the requested 1e-3. So they are not a population sitting
just under the threshold — the cumulative scan is not close to failing on any of them.

The rate is flat in dimension — 16 / 11 / 11 / 13 silent misses at d = 2 / 3 / 4 / 5 — which is
independent confirmation that none of §3.1–§3.3 is a flavour-structure effect.

This is simultaneously the strongest evidence *for* these eleven commits and the sharpest
statement of the pre-existing defect they do not reach.

### 3.4 The 2 ≤ N < 25 band was **not** improved by commit 9c7945a

The brief asks directly: "Commit 8 should have improved this band too; confirm it did." **It did
not.** Multi-resonance profiles through the wrapper, scored against `solve_ivp` (Battery 3.2):

| n_cycles | N=2 | N=5 | N=8 | N=16 | N=24 | **N=25** | N=40 |
|---|---|---|---|---|---|---|---|
| 4 | 4.83e-04 | 4.83e-04 | 9.53e-04 | 9.88e-04 | 8.17e-04 | **7.80e-09** | 5.71e-09 |
| 6 | 4.27e-07 | 1.82e-03 | 1.25e-03 | 1.82e-03 | 1.58e-03 | **2.41e-09** | 3.56e-09 |
| 10 | 3.08e-06 | 1.69e-05 | 1.59e-04 | 4.39e-03 | 4.88e-03 | **2.16e-09** | 2.09e-09 |

Below the seam the band sits at 5e-04 – 4.9e-03, i.e. **outside the requested 1e-3** in six of
21 cells. These do emit `MagnusConvergenceWarning`, so they are not silent — but as
`DECISION_DISPATCH_ORDER.md` §5 already notes, that warning "also fires on rows accurate to
1.6e-06, so it does not discriminate". The branch's routing change rescues N ≥ 25 completely;
everything below it is exactly as exposed as on `main`.

### 3.5 Large accuracy discontinuities at the N = 25 seam

Criterion: a large accuracy discontinuity at the seam is a defect even if both sides are inside
tolerance. Measured (Battery 3.1):

| profile | err(N=24) | err(N=26) | ratio |
|---|---|---|---|
| solar exponential | 3.30e-05 | 2.13e-08 | **1 546×** |
| multi-resonance | 1.58e-03 | 2.86e-09 | **552 945×** |
| noisy | 6.27e-04 | 1.04e-08 | **60 418×** |
| castle wall + breakpoints | 2.80e-11 | 2.80e-11 | 1.0× (takes the cumulative path from N=2) |

Adding one baseline to a 24-point scan changes the answer by up to 1.6e-03 — across the
requested tolerance. This is the *intended* consequence of
`HYBRID_YIELDS_TO_CUMULATIVE_MIN_POINTS = 25`, and the jump is always toward the truth, so it is
a documentation matter rather than a numerical one. It deserves a sentence in the constant's
docstring, because a user scanning N and watching their answer step by 10⁵ in accuracy will
otherwise assume something is broken.

### 3.6 `cumulative=False` is unreachable through the wrapper layer

```
>>> osc_prob_matter_std_potential(..., cumulative=False)
TypeError: magnus.oscprob.osc_prob_energy_baseline() got multiple values for
           keyword argument 'cumulative'
```

The three wrappers now pass `cumulative=` positionally-by-name themselves
([oscprob.py:4923](src/magnus/oscprob.py:4923), and the two sibling sites), so a caller's
`cumulative` in `**kwargs` collides. `DECISION_CUMULATIVE_DEFAULT.md` §5 says "`cumulative=False`
reproduces the previous behaviour exactly", and §7 names exactly one mitigation for the
bit-reproducibility risk — a caller "who does not know about `cumulative=False`". That mitigation
is unavailable at the layer where the change is visible. `strategy='magnus'` is an equivalent
escape hatch and is tested, but it also opts out of hybrid entirely, so it is not a like-for-like
substitute. One-line fix (pop it from `kwargs` before forwarding); not applied here.

---

## 4. What passed, and what that buys

**Battery 1 — regression.** All nine bit-identity rows are **bit-identical** to `main`:
`strategy='magnus'` scan, sub-threshold scan (N=8), single point, vacuum scan, constant-density
scan, energy scan at fixed baseline, `osc_prob_earth` PREM, `average=True`, explicit
`cumulative=False`. The brief flagged that commit 9c7945a might legitimately have moved the
single-point and sub-threshold rows; **it did not** — on smooth solar profiles the γ sweep finds
nothing the extrema did not. The two rows that moved are exactly the two that should:

| row | `main` | branch | |
|---|---|---|---|
| 3ν solar scan, N=60 | 6.130e-04 | **3.504e-09** | 174 936× better |
| ν̄ solar scan, N=60 | 1.145e-05 | **1.006e-08** | 1 138× better |

**Battery 6 — flavor count (the required gap).** 72 configurations scored against `solve_ivp`
across d = 2,3,4,5 × {std, NSI, LIV(n_liv=0,1)} × {ν, ν̄} × N either side of the seam. Worst
error **6.049e-04**, zero silent misses, and **no monotone degradation with d** (median error by
dimension: 8.3e-07, 8.5e-05, 9.2e-05, 9.4e-05). Specifically:

- **No over-merging across level pairs** (6.3): at d = 2…5 the candidate count scales with the
  pair count (12, 54, 116, 182) and a dense per-pair γ scan finds **zero** exceedances outside a
  reported window. The failure mode the brief called untested is not present.
- **Degenerate levels** (6.4): `D41 == D31` exactly, and at 1+1e-9 / 1+1e-14, and `D21 == 0`, all
  return 2e-09 – 1.2e-08 with no `inf`/`nan` leaking out of `_point_adiabaticity`.
- **γ-sweep cost does not grow with dimension** (6.5): 1.16× at d=2, 0.98× at d=3, 1.09× at d=4,
  1.11× at d=5 against the extrema-only detector. The brief expected the 1.31× at 2ν to grow with
  10× the pairs at d=5; it does not, because the sweep is one vectorised `einsum` pass, not one
  pass per pair. **Suspicion #4 is refuted.**
- **`ip_exp` still declines for d > 2** (6.7) at all three dispatch sites, and for LIV at every d.

**Battery 4 — extreme numerics.**

- **The commit-5 suppression did not swallow the signal that matters**: `ToleranceNotAchievedWarning`
  fires correctly at requested 1e-9 and 1e-12, and never at 1e-3/1e-6 where the tolerance is met.
- **`max_n_slabs` lowered (suspicion #3)** — the coupling is real, but it is **not silent**:

  | `max_n_slabs` | probe `n_slabs` | scan error | warned? |
  |---|---|---|---|
  | 500 | 500 | 6.194e-02 | `ToleranceNotAchieved` |
  | 2000 | 2000 | 4.657e-03 | `ToleranceNotAchieved` |
  | 5000 | 5000 | 3.053e-06 | `ToleranceNotAchieved` |
  | 20000 (default) | 20000 | 9.546e-09 | — |

  `DECISION_CUMULATIVE_DEFAULT.md` §4e says the degradation happens "without a separate warning".
  It is worth amending: there is no *cumulative-specific* warning, but the probe's
  `ToleranceNotAchievedWarning` does reach the caller at every capped level.
- **Geometry** (4.5): `L0` mid-profile, spans of eight orders of magnitude, unsorted, duplicated,
  degenerate and `L == L0` baselines all correct to ≤ 2.5e-08.
- **Energy** 0.05 MeV – 100 GeV clean (worst 1.261e-04, at 1 GeV over 0.5 R_sun). At 0.05 MeV the
  oracle-convergence check honestly reports **ORACLE SUSPECT** — it moves by 1.5e-09 against a
  quoted error of 1.57e-09 — so that row should be read as "error ≤ ~2e-09, at the oracle's own
  floor", not as a measurement.
- **Memory** (4.2): peak 10.4 / 17.1 / 84.1 / 194.5 MB at N = 10³ / 10⁴ / 10⁵ / 10⁶, against
  outputs of 0.03 / 0.3 / 3.2 / 32 MB. Bounded and modest, but **not flat in N** — the
  non-output part grows 4.8× from N=10⁴ to 10⁵. The "O(block) + O(result)" claim is right in
  shape (the grid necessarily carries the N requested baselines) but "flat at ~56–59 MB" does not
  extrapolate past N = 4000.

**Battery 5 — the constructions that did *not* break it.**

- **The frozen-grid mode could not be reproduced.** A piecewise-constant profile with the
  interior breakpoints supplied but `l_ini`/`l_fin` missing — the exact
  `NOTES_ADAPTIVE_REFINEMENT.md` §4b construction — gives **6.658e-11** and raises
  `ToleranceNotAchievedWarning`. The strict probe runs to `max_n_slabs = 20000` rather than
  freezing, and `CUMULATIVE_N_ACC_SAFETY = 4` carries the rest. **Commit 051a3c7's always-strict
  probe works as designed**, and this is direct evidence for it.
- **Determinism is exact** (5.6): repeat call, shuffled input order, and `n_jobs` = 2 and 4 all
  give `max|diff| = 0.000e+00`.
- **Adversarial `n_slabs`** of 1, 2 and 10⁶ (× `n_tpts_per_slab` 2 and 500) all return the same
  correct answer: with a tolerance requested `n_slabs` is only a floor, and 10⁶ is clipped at the
  cap rather than allocated.

**Battery 7 — oracle diversity.** `expm` for constant `H`: 1.5e-13 (d=2) to 1.8e-11 (d=3) over 60
baselines, and the cumulative scan called directly agrees with `expm` to 1.6e-13 (d=2) and
6.4e-12 (d=3) — time-ordering and indexing are right at every dimension. Analytic 2ν vacuum
formula: 9.6e-14 over 200 baselines. Composition law `U(0→L₂)` vs `U(L₁→L₂)·U(0→L₁)`: 5.4e-13 to
6.4e-11. `average=True` against a brute-force 25-oscillation window average: **4.09e-06**
(better than the 1.5e-04 previously recorded).

**Unitarity** (measured on the package's own output via probability row/column sums, not on the
oracle's operator): worst **1.643e-11**, at d=2, N=10⁵. Across d = 2…5 and N = 25 … 10⁵ it
degrades only from ~3e-12 to ~1.6e-11 — four decades of N for half a decade of unitarity.
**Suspicion #5 is not borne out**; there is a long way to the 1e-9 the suite asserts.

---

## 5. Attribution — the most important table here

Every failing configuration, run under `main` and under the branch, dumping P_ee, the window
count and the certification flag:

| case | `main` P_ee / win / cert | branch P_ee / win / cert | |
|---|---|---|---|
| sub-threshold, w=3e-2 span | 0.8071837 / 0 / True | 0.8071837 / 0 / True | identical |
| sub-threshold, w=1e-2 span | 0.9418934 / 0 / True | 0.9418934 / 0 / True | identical |
| detection miss, w=3e-5 span | 0.0968947 / 0 / True | 0.0968947 / 0 / True | identical |
| detection miss, w=1e-5 span | 0.0968947 / 0 / True | 0.0968947 / 0 / True | identical |
| sinusoid, period span/7 | 0.9014671 / 0 / True | 0.9014671 / 0 / True | identical |
| **step function, edge unmarked** | 0.4975934 / 0 / True | 0.4975934 / 0 / True | identical |

**Bit-identical in every case.** Two conclusions, and both matter:

1. **The branch causes none of these.** They are `magnus.adiabatic` defects that predate it, and
   `DECISION_CUMULATIVE_DEFAULT.md` §4g already logs the class ("it certifies while wrong") as
   pre-existing and out of scope. That log is confirmed, and its scope widened: the mechanism is
   not confined to multi-resonance NSI profiles.
2. **Commit 9c7945a's fix is narrower than its message.** "Look for non-adiabaticity along the
   path, not only at gap extrema" is accurate about *where* γ is sampled and does close the
   4.3e-02 case it was written for. It does not address γ never crossing the initial threshold,
   which §3.2 shows is the larger exposure, nor a feature narrower than the probe spacing (§3.3).

---

## 6. Recommendations, in priority order

**All six were subsequently implemented — see §8 for what changed and what it measured.**
They are left here as originally written, so that §8 can be read against the reasoning that
produced it rather than against a tidied-up version of it.

1. **Guard the hybrid path against non-smooth profiles by *detecting* them, not by asking.**
   The dispatcher already computes the Hamiltonian on a 200-point probe grid; a cheap
   second-difference or total-variation test on `VCC` along that grid would catch §3.1 (a jump of
   15× across one probe interval) at essentially no cost, and decline to the general path, which
   is correct there. This is the single highest-value change available.
2. **Make certification depend on a window having been *considered*, not merely on two adiabatic
   answers agreeing.** Concretely: if no window opens at `threshold0`, require at least one
   further iteration at a lower threshold before the agreement test may certify. §3.2's table
   shows the answer is already correct once any window opens, so the fix is about *when* the loop
   is allowed to stop, not about the patching machinery.
3. **Amend the three documentation claims the measurements contradict:**
   - [adiabatic_strategy.rst:291](docs/source/adiabatic_strategy.rst:291) — "composes correctly
     with any number of simultaneous [or sequential resonances]". Already flagged false in
     `DECISION_CUMULATIVE_DEFAULT.md` §4g; still present.
   - [adiabatic_strategy.rst:125](docs/source/adiabatic_strategy.rst:125) — scanning every pair
     "so any number of simultaneous or sequential resonances are all found". The per-*pair*
     half of this is true and now measured (Battery 6.3, zero uncovered exceedances at d=2…5);
     the "all found" half is false for any feature narrower than the probe spacing (§3.3). The
     limit is `n_probe`, and saying so costs one sentence.
   - `DECISION_CUMULATIVE_DEFAULT.md` §4e — "without a separate warning" for lowered
     `max_n_slabs` (§4 above).
4. **Let `cumulative` through the wrapper layer** (§3.6) — a one-line `kwargs.pop`.
5. **Record the N = 25 accuracy step** in `HYBRID_YIELDS_TO_CUMULATIVE_MIN_POINTS`'s docstring.
6. **Then split the branch** as the handover suggests: commit 9c7945a is a different module from
   the other ten. Nothing found here argues against the split, and item 2 would land naturally on
   the `magnus.adiabatic` side of it.

---

## 7. What this validation did *not* cover

Stated so that "we found nothing" can be weighed against what was actually probed.

- **Notebooks 02 and 03 were not re-executed** (Battery 1.4). The figure-hashing and
  per-figure `solve_ivp` checks recorded in `DECISION_DISPATCH_ORDER.md`'s appendix were taken as
  given; this session tested the library paths those cells exercise, not the cells.
- **Timing was not re-measured.** The brief's performance criterion (no configuration >10 % slower
  than `main`, by alternating) was not run: the machine's documented run-to-run variance makes it
  a multi-round exercise, and every defect found is an accuracy defect. The incidental timings
  quoted above are single-shot and should not be treated as measurements.
- **`trapezoid`/`simpson`** were left where `DECISION_CUMULATIVE_DEFAULT.md` §4f left them.
- **The fuzzer drew smooth profiles only** (random Fourier sums). Given that §3.1 is a
  discontinuity defect, a fuzzer over *piecewise* profiles is the obvious next thing to build and
  would likely be productive.

**Every item in this section has since been closed — see §9.**

---

## 8. What was implemented, and what it measured

All six recommendations of §6 were applied. Two are algorithmic changes to `magnus.adiabatic`;
the rest are one-line or documentation changes. Reproductions live in
`adversarial_batteries/`, and the new behaviour is pinned by five tests in
`tests/test_adiabatic.py`.

### 8.1 Two changes to `magnus.adiabatic`

**A probe-scale resolution test** (`_profile_is_resolved`, rec 1). The guard against a non-smooth
profile was the proxy "did the caller pass `t_breakpoints`", which fails open. It is now a
measurement: halve the probe spacing and ask how much the largest adjacent change in `H` shrinks.
For a `C¹` Hamiltonian it halves (ratio → 0.5); across a jump the finer grid still straddles the
jump, so it does not (ratio → 1.0). `RESOLUTION_RATIO = 0.75` sits between the two limits, and
every profile measured falls decisively on one side.

*The first version of this was wrong, and the batteries caught it.* Testing only at `n_probe0`
cannot distinguish a genuine jump from a feature that is merely sharp at that density — and a
Gaussian of width 1e-3 of the domain, which the module answers **exactly** (1.09e-11) once
refinement resolves it, was abandoned as though it were a step function. The test is now
two-stage: the cheap check at `n_probe0`, and only if that fails, a confirmation at
`max_n_probe`. A jump is unresolved at every density; a sharp feature is not. Ordinary calls
never pay for the second stage.

**γ-aware certification** (rec 2). When no window opens, successive iterations differ only in the
adiabatic-transport grid, so they converge to the same adiabatic limit and agree with each other
whether or not that limit is right — the agreement test carries no information about the thing
that went wrong. Certifying an *empty* window list now additionally requires γ itself to be small
enough for the requested tolerance; otherwise the threshold keeps dropping until a window opens,
which is guaranteed because γ_max is measured on the same grid the threshold is compared against.
`GAMMA_TO_ERROR = 1.0` comes from the measured ratio |ΔP|/γ_max ∈ [0.30, 0.55] over two decades,
rounded conservatively. `find_nonadiabatic_windows` reports γ_max through a new optional `info`
dict, following the `convergence_info` convention already used in `osc_prob`.

### 8.2 What it fixed, measured through the public entry point

`osc_prob_matter_std_potential(..., strategy='auto')` — the default — against `solve_ivp`, or
against `expm` where the profile makes `expm` exact:

| case | before | after | |
|---|---|---|---|
| **step function, edge unmarked** | **5.395e-01, silent** | **3.610e-04, warns** | falls back to the general path |
| kink, C⁰ but not C¹ | 1.448e-02, silent | **5.845e-10** | window now opens |
| singularity approached | 8.625e-03, silent | **3.057e-05** | window now opens |
| sinusoid, period span/7 | 1.672e-02, silent | **4.444e-10** | window now opens |
| sub-threshold bump, w=1e-2 span | 7.701e-03, silent | **2.833e-09** | window now opens |
| sub-threshold bump, w=3e-2 span | 4.388e-03, silent | **1.104e-09** | window now opens |
| narrow bump, w=3e-5 span | 2.907e-02, silent | 2.907e-02, silent | **unchanged — see below** |

Battery 2, re-run in full over the same 42 configurations: **12 certified-but-wrong → 4**, and
all four survivors are the same narrow-feature case (w = 1e-5 and 3e-5 of the domain, at
d = 2, 4 and 5). Sub-test by sub-test:

| sub-test | before | after | |
|---|---|---|---|
| 2.1 narrow-resonance width sweep | 4 | 2 | both survivors w ≤ 3e-5 |
| 2.2 aliasing | 1 | **0** | the span/7 sinusoid, 1.672e-02 → 4.444e-10 |
| 2.3 edge crossings | 2 | **0** | 5.365e-03 → 7.772e-10 |
| 2.4 many crossings | 1 | **0** | 3.907e-02 → **4.782e-09**; this was the worst error found anywhere |
| 2.5 clustered crossings | 0 | 0 | already clean |
| 2.7 d = 4 and 5 | 4 | 2 | w=1e-2 cases fixed (7.7e-03 → 1.7e-09, 7.2e-03 → 1.1e-09) |

### 8.3 What was deliberately *not* fixed, because it cannot be

**A feature narrower than the probe spacing.** The two survivors above are a Gaussian of width
3e-5 and 1e-5 of the domain. Neither the probe grid (spacing 5e-3 of the domain) nor its
refinement ceiling (1.6e-4) samples it at all, so `H` looks perfectly smooth, γ looks small, and
the resolution test correctly reports "resolved" — there is genuinely no information at that
sampling density. The general Magnus path misses it too (2.907e-02), differing only in that it
warns. This is a property of any fixed grid, not of any of these tests, and the honest response
is the documentation change in rec 3 rather than a detector that pretends otherwise. The cure
available to a caller is `t_breakpoints` at the feature, or a larger `n_probe`.

### 8.4 The three other recommendations

**Rec 4** — `_resolve_cumulative_kwarg` pops a caller-supplied `cumulative` before the wrappers
set their own, so `cumulative=False` is now reachable from the wrapper layer instead of raising
`TypeError`. An explicit value wins over the `strategy='magnus'` opt-out, on the grounds that
naming `cumulative` is a specific request and `cumulative=True` raises rather than silently
falling back. Verified through all three wrappers: `cumulative=False` reproduces the general-path
answer exactly (4.859e-03, bit-identical to `strategy='magnus'`), `cumulative=True` and `'auto'`
give the cumulative answer (3.349e-08).

**Rec 5** — the N = 25 accuracy step is now recorded in
`HYBRID_YIELDS_TO_CUMULATIVE_MIN_POINTS`'s docstring, with the measured table and the note that
the step is always toward the truth.

**Rec 3** — `adiabatic_strategy.rst` now separates the claim that survived measurement (every
level pair is scanned, so simultaneous resonances do not mask each other — verified at d = 2…5
against a dense per-pair γ scan) from the two genuine limits, which are stated in a `warning::`
block. `DECISION_CUMULATIVE_DEFAULT.md` §4e's "without a separate warning" is amended with the
measured table showing `ToleranceNotAchievedWarning` does fire at every capped `max_n_slabs`, and
§4g's forward reference is marked resolved.

### 8.5 Regression coverage

Five tests in `tests/test_adiabatic.py` pin the new behaviour, each on the construction that
produced the defect:

| test | pins |
|---|---|
| `test_unmarked_discontinuity_is_detected_and_not_certified` | the 0.54 case; oracle is `expm`, exact here |
| `test_a_smooth_profile_is_still_reported_as_resolved` | the mirror — ordinary solar profiles, and a constant `H`, must keep certifying |
| `test_a_sharp_but_smooth_feature_is_not_mistaken_for_a_discontinuity` | the two-stage test; asserts unresolved at 200 **and** resolved at 6400 |
| `test_subthreshold_nonadiabaticity_is_not_certified_on_agreement_alone` | γ under `threshold0` but over the tolerance; asserts a window opens and the answer is right |
| `test_find_nonadiabatic_windows_reports_gamma_max_via_info` | the `info` out-parameter, including that omitting it still works |

---

## 9. Closing the gaps §7 recorded

### 9.1 Notebooks 02 and 03, re-executed

Both run cleanly under `nbclient` against the fixed code — 101 cells in 482 s and 108 cells in
478 s. Hashing every embedded PNG:

| | stored | re-run | changed | added | removed |
|---|---|---|---|---|---|
| notebook 02 | 14 images | 14 images | **none** | none | none |
| notebook 03 | 20 images | 20 images | **none** | none | none |

Warning output is identical too — cells 36, 41, 55, 69, 97 in notebook 02 and 43, 48, 62, 76, 104
in notebook 03, with **none gained and none lost**. Since no figure moved, the per-figure
`solve_ivp` re-scoring that Battery 1.4 called for has nothing to score: the stored figures are
still exactly the ones the earlier work validated.

(Per-cell timing hooks `async_execute_cell`, not the sync `execute_cell` — the handover's trap,
which records nothing while reporting success.)

### 9.2 Piecewise-profile fuzzing

150 random piecewise-constant profiles — random segment count (2–12), random breakpoint
positions, random densities over two decades, random d ∈ {2,3,4,5}, random energy, ν/ν̄, and scan
size. `expm` composed across the segments is the **exact** operator for these, so the oracle is
not an approximation and cannot itself step over a jump.

| | median | p90 | max | outside 1e-3 | **silent** |
|---|---|---|---|---|---|
| edges **not** declared | 7.76e-04 | 2.96e-03 | 3.46e-02 | 59 / 150 | **2** |
| edges declared (`t_breakpoints`) | **1.34e-12** | 1.40e-11 | 7.19e-11 | **0** | **0** |

Declaring the edges is never worse (0 / 150 cases more than 10× worse) and is essentially exact.
Undeclared, the answer is often inaccurate — which is expected and is what the warnings are for —
but only **2 of 150** are silently so, at 2.1e-03 and 1.4e-03, both at N = 80 on the cumulative
path rather than the hybrid one. That residue is a fair characterisation of what is left: a
uniform grid that happens not to align with a jump, overrunning the tolerance by about 2×.

### 9.3 Performance, measured by alternating — and the criterion fails, on `main`-vs-branch

Seven rounds, three trees interleaved round-robin (never back to back), two workloads the branch
cannot touch carried as controls. Medians, in ms, with the ratio against `main`:

| workload | main | branch (pre-fix) | branch + fixes | vs main |
|---|---|---|---|---|
| single point, solar (hybrid) | 22.6 | 34.0 | 32.0 | **1.42×** |
| single point, 3ν solar | 28.7 | 38.9 | 40.2 | **1.40×** |
| single point, multi-resonance | 1255.5 | 1202.6 | 1154.9 | 0.92× |
| solar scan N=400 (cumulative) | 8116.2 | 105.5 | 103.2 | **0.01×** |
| solar scan N=8 (hybrid) | 163.1 | 246.3 | 247.8 | **1.52×** |
| CONTROL vacuum scan N=300 | 31.0 | 32.1 | 31.7 | 1.02× |
| CONTROL const-density scan N=300 | 30.8 | 32.3 | 31.5 | 1.02× |

The controls come back at 1.02×, and the machine's own scatter on them is 1.05–1.15× within a
single tree, so anything under about 1.15× is noise. Two things follow, and they are different
things:

**The fixes in §8 are free.** Column 3 against column 2 is 0.94×, 1.03×, 0.96×, 0.98×, 1.01× —
inside the control scatter everywhere. The first version was not: `_profile_is_resolved` called
`H_func` one position at a time, costing 1.2× on an ordinary single-point call. These
Hamiltonians are array-capable — it is the same fast path `magnus.magnus` relies on — so
`_H_on_grid` now makes one vectorized call and falls back to the loop only for a scalar-only
`H_func`. Measuring rather than assuming is what caught this.

**The branch fails the 10% performance criterion against `main`, and did so before these
fixes.** Hybrid-answered single points and sub-threshold scans are **1.40–1.52×** slower, from
22.6 → 32.0 ms and 163 → 248 ms. §7 recorded timing as not covered, so this had never been
quantified at the entry point. The cause is commit 9c7945a's γ sweep, whose own commit message
claims 1.31× measured inside the detector; 1.4–1.5× at the wrapper is consistent with that.

This is a real criterion failure and is reported as one. It is also, on the evidence here, the
right trade: that sweep is what makes the difference between certifying a 4.3e-02 error and
detecting it, and the same code path now answers the multi-resonance point *faster* than `main`
(0.92×) and a 400-point scan **100× faster**. But it should be a decision taken knowingly rather
than an unmeasured side effect, which is what it was until now.

### 9.4 Quadrature methods

`gl`, `trapezoid` and `simpson` × `strategy` ∈ {auto, magnus}, on solar profiles: **no silent
miss**. The one row outside tolerance (`trapezoid`/`magnus`, 4.893e-03) warns. Under `auto` all
three integration methods return the same answer, because no window opens and the adiabatic
transport uses its own Simpson quadrature regardless — worth knowing, and not previously stated.
This leaves `DECISION_CUMULATIVE_DEFAULT.md` §4f's conclusion intact: these methods are
*inaccurate* on hard profiles, but they are not silently inaccurate.

---

## 10. The split (rec 6), and why it is stacked rather than parallel

Done, but not in the shape §6 assumed. The recommendation reads "commit 9c7945a is a different
module from the other ten", and in *source* that is exactly right: after the split,

| branch | touches |
|---|---|
| `cumulative-and-notebooks` | `oscprob.py`, notebooks 02/03, `test_oscprob.py`, the decision docs |
| `adiabatic-certification-honesty` | `adiabatic.py`, `test_adiabatic.py`, `adiabatic_strategy.rst` |

with no overlap in either direction — except for fifteen lines of `test_oscprob.py`.

Those fifteen lines are why the split is **stacked** (the adiabatic branch branches *from* the
cumulative one) rather than two independent branches off `main`. Cherry-picking 9c7945a straight
onto `main` conflicts, because the test it modifies —
`test_baseline_scan_across_many_resonances_matches_solve_ivp` — was added by 2624bbe on the
cumulative side and is built on `HYBRID_YIELDS_TO_CUMULATIVE_MIN_POINTS`, which does not exist on
`main`. The multi-resonance regression test that motivated the γ sweep is therefore genuinely
downstream of the cumulative routing, and pretending otherwise would mean rewriting a test to
make a branch topology look tidier than the work was.

Verification that the split preserves the work: `git diff adiabatic-certification-honesty
notebook-breakpoints-and-cumulative` is **empty** — the stacked pair reproduces the validated
tree exactly, commit for commit in content. The original branch is left in place unchanged as the
record of what was actually validated; nothing is pushed.
