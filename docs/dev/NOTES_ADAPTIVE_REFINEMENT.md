# Notes: what the `n_slabs` floor fixed, and what it did not

**Written:** 2026-08-02, executing §6 of `DECISION_OSCPROB_BATCHING.md`.
**Machine:** the same 12-core box the decision doc used; every number below is a
min-of-2-or-3.

---

## 1. The mechanism is not aliasing

The decision doc (§2) diagnoses the false convergence as aliasing: "the coarse grid
*aliases* the periodic profile: consecutive refinements agree with each other and
disagree with reality." That framing predicts a systematic near-agreement between
neighbouring refinement levels. The ladder does not behave that way. Fixed-grid
answers on notebook 03's castle wall at L = 5242.4 km, truth 0.986024:

| `n_slabs` | 1 | 2 | 3 | 4 | 5 | 6 | 8 | 12 | 18 | 27 | 41 | 62 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| P | 0.478 | 0.428 | 0.131 | 0.131 | 0.645 | 0.117 | 0.680 | 0.956 | 0.747 | 0.991 | 0.980 | 0.988 |

Below about 40 slabs — roughly 1.5 per wall, with 26 walls crossed — the sequence
**thrashes**. It does not converge slowly and it does not alias consistently. The
3-vs-4 agreement that `np.allclose` fired on (2.3e-4, comfortably inside the default
1e-3) is a **coincidence** between two adjacent samples of a wildly varying sequence;
its neighbours are 0.645 and 0.117.

This matters for three reasons.

1. It explains why a tighter `rtol` cannot help, and does so more sharply than the
   aliasing story: there is no systematic error being under-resolved, there is a
   coincidence being trusted.
2. It says what *would* have caught it: requiring agreement **twice in a row**. See §3.
3. It means the failure is not restricted to periodic profiles. Any profile whose
   feature scale the seed underestimates puts the ladder in this regime.

## 2. Why the "make `suggest_n_slabs` feature-scale-sensitive" option is not the cheap one

§6.1 of the decision doc offers two alternatives and calls the second the cheap one.
Having measured both, the first is not merely more expensive — the obvious cheap
version of it does not work at all.

`suggest_n_slabs` estimates the accumulated phase from 17 probe points. The natural
cheap fix is to check that estimate against a doubled probe and refine until it
stabilises — which would detect a profile the probe cannot see. Measured on the
castle wall, the estimate is *already stable*:

| probe points | 17 | 33 | 65 | 129 | 257 | 513 | 1025 | 2049 |
|---|---|---|---|---|---|---|---|---|
| phase [rad] | 9.353 | 9.411 | 9.411 | 9.395 | 9.366 | 9.366 | 9.366 | 9.368 |

The seed is not wrong because it is under-sampled. It is wrong because **the phase
criterion is the wrong criterion for this profile**: `suggest_n_slabs` measures the
*integral* of the Hamiltonian along the path, and an integral is blind to structure
that averages out. Fifty walls oscillating about their mean accumulate ~9 radians
net, so "2π radians per slab" honestly returns 2 slabs. The estimate is accurate;
the inference from it is not.

A criterion that would work has to measure *variation* rather than the integral —
total variation along the path, or the deviation of the profile from its running
mean — and then calibrate how much variation demands how many slabs. That is a new
accuracy heuristic on the package's highest-traffic entry point, not a tweak. It was
not built.

## 3. What the floor does not cover, and the cost of covering it

The floor makes a caller-supplied feature scale authoritative. A caller who supplies
none still gets the old ladder. `test_unfloored_ladder_is_what_the_floor_protects_against`
pins this explicitly.

The measured option for closing that gap is **two consecutive agreements** before
declaring convergence — the safeguard `_osc_prob_ip_exp_core` already applies, for
the same class of failure (see its comment at `oscprob.py:2795`). Against the ladder
in §1 it works: 3-vs-4 agrees, 4-vs-6 does not (1.4e-2), so the loop would keep
climbing.

It was not built here because the cost is broad and this brief did not ask for it.
Every adaptive call would run one extra refinement loop, and since loop *n* costs
about 1.5× loop *n−1*, the total lands near 1.6× for **every** tolerance-driven call
in the package — against a benefit that only materialises on profiles whose feature
scale the seed underestimates. That trade deserves its own decision, with the same
scrutiny this one got. It is the natural next item.

## 4. Numbers this work produced

Castle-wall baseline scan, 6000 points over 100–10⁴ km, sampled at 200–300 of them,
against a converged (`n_slabs=24000`, `rtol=None`) reference itself checked against
`solve_ivp`/DOP853 at the worst points:

| | worst \|ΔP\| | points > 1e-3 | cost |
|---|---|---|---|
| as shipped (`n_slabs=150` discarded) | **0.855** | 7.2% > **1e-2** | 2.6 ms/pt |
| `n_slabs=150` honoured as a floor | 1.5e-2 | 10.0% | 3.8 ms/pt |
| floor + `t_breakpoints` at the walls | **1.9e-3** | 1.5% | **1.4 ms/pt** |

The third row is the one the notebooks now use. Breakpoints are worth more than
slabs here because the profile is a step function: a slab straddling a wall degrades
the quadrature to low order regardless of `magnus_exp_order`, so an edge placed *on*
the wall buys accuracy that uniform refinement has to pay for many times over. Across
the four castle-wall configurations in notebooks 02 and 03 (narrow/wide × baseline/
energy) breakpoints improve the worst case by 3.6–7.5× and run 1.9–4.7× faster.

### Every figure the floor moved, it moved toward the truth

The floor changes more than the castle-wall figures: every notebook cell that passes
`n_slabs > 1` with a tolerance is now refined at least that finely. Nine stored figures
changed in each notebook. Each affected profile family was checked against a
`solve_ivp`/DOP853 reference, comparing the old ladder (equivalent to passing
`n_slabs=1`) with the floored one:

| profile (notebook 03 cell) | old error | new error | |
|---|---|---|---|
| exponential density, 200–1000 km (43) | 1.0e-4 – 3.9e-4 | 5.9e-9 – 2.9e-6 | better |
| Gaussian density, 200–1000 km (43) | 1.3e-4 – 4.3e-4 | 4.3e-9 – 8.1e-8 | better |
| PREM, costhz −0.1/−0.5/−1.0 (85) | 2.6e-4 – 1.9e-3 | 3.1e-5 – 7.7e-4 | better |
| solar exponential (99) | 1.3e-5 – 2.1e-5 | unchanged | its ladder already exceeded 200 |

Nothing measured got worse. This is worth stating explicitly because error is *not*
monotone in `n_slabs` (see the ladder in §1), so "more slabs is better" is an empirical
claim here, not an automatic one.

### The notebooks got *faster*, which settles decision-doc §5.3

§5.3 named the scenario that could flip the batching decision on its own: "the accuracy
fix lands and turns out to make the per-point path much more expensive … a correct
per-point loop might cost 5-10× what today's wrong one does, at which point 115 s
becomes 600-1000 s and proposal (3) is suddenly the majority of the notebook. I would
re-run §3.1 immediately after the accuracy fix."

Re-run, on the same machine:

| | before | after |
|---|---|---|
| notebook 03 | 652 s | **445 s** (−32%) |
| notebook 02 | — | 394 s |

The honest fix made the notebooks *cheaper*, not 5-10× dearer, because the breakpoints
buy back more than the floor costs. §3.1's split moves accordingly — notebook 03's
energy scanning goes 504.8 s → 320.9 s and its baseline scanning 115.0 s → 93.9 s, so
baseline scanning is 21% of a much shorter notebook rather than 17.6% of a long one.

**This is the scenario §5.3 flagged, and it resolved against proposal (3), not for it.**
The condition that would have made the case was the per-point path getting much more
expensive; it got cheaper. Nothing here argues for building the cumulative scan sooner.

Both notebooks re-ran with no errors, and the only change to their warning output is
that five cells stopped emitting `MagnusConvergenceWarning` — none started.

### The stale performance comments (§6.4)

Eight cells asserted ~4.5 ms per `osc_prob` call, "almost entirely fixed entry-path
cost", evidenced by "n_slabs=1 and n_slabs=150 both measure 4.5 ms per call". All
three claims are false on this machine today:

```
osc_prob on a varying profile, fixed grid:
  n_slabs=1      0.343 ms/call
  n_slabs=150    1.141 ms/call
  n_slabs=2000  10.364 ms/call
```

Loop-vs-array at each cell's real point count, measured: 1.35×, 1.60×, 1.66×, 1.99×,
2.21×, 2.25×, 2.27×, 2.28× — against the ~22× the comments predicted. The array call's
cost is genuinely per-point (0.13 ms/pt, flat from N=100 to N=40000), so the entry
path is not a fixed cost waiting to be amortised. Passing the array is still the right
advice, for a smaller and correctly-stated reason.

## 4b. The two-consecutive-agreements safeguard, re-derived (2026-08-03)

§3 estimated the cost at ~1.6× and left the trade undecided; `HANDOVER_DISPATCH_AND_ADOPTION.md`
Task 3 then recommended against building it. Both predate the measurement below, which was
prompted by finding the general path silently outside its requested tolerance at the **default**
tolerance while measuring the dispatch reorder (`DECISION_DISPATCH_ORDER.md` §5).

Method: replay the ladder rung by rung (seed from `suggest_n_slabs`, then
`n_slabs <- round(1.5*n_slabs)`), scoring every rung against `solve_ivp`/DOP853 at `rtol=1e-12`.
The replay was validated by checking that the rung where a *single* agreement first fires
reproduces the shipped adaptive answer **bit-for-bit**; it does, on both points below.

### The mechanism is confirmed, and it is a coincidence between two wrong rungs

10 MeV over R_sun, `rtol=atol=1e-3`, seed 977:

| k | 0 | 1 | 2 | **3** | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|---|---|
| `n_slabs` | 977 | 1466 | 2199 | **3298** | 4947 | 7420 | 11130 | 16695 | 20000 |
| error | 5.9e-02 | 3.8e-03 | 1.6e-02 | **1.7e-02** | 8.1e-03 | 4.5e-03 | 3.5e-06 | 1.6e-07 | 6.8e-08 |

`np.allclose` fires at k=3 because rungs 2 and 3 agree to 1.088e-03 — while both are wrong by
~1.6e-02, and the next rung moves by 2.5e-02. Exactly the §1 mechanism, at the default tolerance
on the package's most ordinary profile.

### It works on the thrashing mode, and the ~1.6× estimate was right — for healthy calls

15 configurations (solar 2nu at nine (E, L) points, castle wall narrow/wide with breakpoints,
noisy high/low amplitude, 3nu exponential and Gaussian):

| | |
|---|---|
| silently outside 1e-3 today | **6 / 15** |
| fixed by two consecutive agreements | **4 / 6** |
| time cost on calls that were already correct | median **1.53×**, max 1.75× |
| time cost on the calls it fixes | 3.5× – 8.3× |

The two costs are different quantities and §3 conflated them. ~1.6× is the *tax* on a call whose
first agreement was genuine — one extra rung. On a call the safeguard actually rescues, the
second agreement is several rungs away, because convergence genuinely has not happened: 8.26× at
10 MeV, 8.25× on the low-amplitude noise. That asymmetry is a feature (the cost is paid where the
answer was wrong), but it must not be quoted as 1.6×.

### The second failure mode is an incomplete breakpoint list, not an algorithmic blind spot

§3 claimed the safeguard "works" on the castle wall: "3-vs-4 agrees, 4-vs-6 does not (1.4e-2), so
the loop would keep climbing." That was measured **without** `t_breakpoints`. With them, the
failure looks completely different — the answer is **frozen**, not thrashing. Castle wall narrow,
49 breakpoints over 50 walls, seed 3:

| `n_slabs` | 4 | 6 | 9 | 14 | 21 | 32 | 48 |
|---|---|---|---|---|---|---|---|
| error | 1.564e-02 | 1.564e-02 | 1.564e-02 | 1.564e-02 | 1.564e-02 | 1.564e-02 | 1.135e-03 |
| \|ΔP\| vs previous | — | 4.7e-16 | 4.4e-16 | 5.6e-16 | 2.0e-15 | 5.1e-15 | 1.5e-02 |

Six consecutive rungs bit-identical and all wrong by 1.6e-02, so **no consecutive-agreement count
escapes it**. That much stands. But the *diagnosis* first recorded here — "a stable fixed point of
the discretisation", implying an algorithmic blind spot — was wrong, and the correction matters
because it dissolves the problem rather than deferring it.

The profile is exactly piecewise constant between consecutive breakpoints (sampled `H_00` spread:
**0.000e+00**) in every interval but the first. The castle wall is defined on
`[l_ini, l_fin] = [100, 10000] km`, and `castle_wall_breakpoints` returns the 49 *interior* walls
— but the trajectory starts at **0**, so the profile switches on at `l_ini`, a genuine
discontinuity with no breakpoint on it. Every level integrates across it; Magnus is exact on all
the other intervals, so refinement changes nothing and the shared answer stays wrong.

Adding `l_ini` and `l_fin` to the list:

| `n_slabs` | 4 | 32 | 150 |
|---|---|---|---|
| as supplied | 1.564e-02 | 1.564e-02 | 1.572e-04 |
| with `l_ini`, `l_fin` | **3.583e-12** | **3.585e-12** | **3.607e-12** |

Nine orders of magnitude, at every slab count. `castle_wall_breakpoints` is a notebook-local
helper, so this is a **notebook bug**, not a library one, and it is not evidence for any change to
the refinement algorithm.

So the failure has two modes, and only the first is the library's problem:

- **thrashing** (§1) — successive rungs vary wildly, `allclose` fires on a coincidence.
  Two agreements catches this: 4/4 of the measured cases.
- **unmarked discontinuity** — the caller's `t_breakpoints` misses a place where the Hamiltonian
  is non-smooth. No refinement strategy can help; the cure is the missing edge. Worth stating in
  the user-facing documentation, which it now is.

**This revises the tally.** Of the 6 configurations recorded above as silently outside tolerance,
2 are this notebook bug. Of the remaining 4, three are solar, where the observable is a
phase-average that removes most of the error (§4c). The genuine, unmitigated library-level
exposure is narrower than the raw count suggests.

## 4c. What the observable is, which changes the cost/benefit again

The errors above are all **pointwise** — P at one (E, L). For solar neutrinos that is not the
observable: at 10 MeV the survival probability oscillates ~2100 times over a solar radius, and
what an experiment measures is the average. Measured on the 10 MeV/R_sun case, averaging over a
window of 25 oscillations near R_sun:

| | pointwise max \|ΔP_ee\| | \|⟨P_ee⟩ − ⟨P_ee⟩_true\| |
|---|---|---|
| ladder stop (`n_slabs` = 3298) | 2.52e-02 | **1.86e-04** |
| converged (`n_slabs` = 20000) | 2.83e-06 | 1.11e-09 |

The error is oscillatory, not a bias (`|mean|/std = 0.014`), so averaging removes ~135× of it and
what survives is *inside* the requested 1e-3. For an averaged solar observable the ladder's
"wrong" answer is fine, and tightening the pointwise criterion buys nothing.

The package already has the right tool: `average=True` answers this case in **36 ms** via
`avgprob.averaged_probabilities_adiabatic`, bypassing the ladder entirely, giving 0.311806
against a brute-force 25-oscillation window average of 0.311956.

Where pointwise accuracy *is* the observable — a probability-versus-baseline or versus-energy
curve, an oscillogram, a fixed (E, L) — none of this applies and the safeguard's value is real.

### Recommendation, revised twice

**Build it, opt-in, off by default** — shipped as `osc_prob(..., strict_convergence=True)`.

The case for a default flipped twice under measurement, and both reversals were right:

1. The raw count was 6 of 15 configurations silently outside tolerance. But 2 of those are the
   notebook's incomplete breakpoint list (§4b), not a library failure.
2. Of the remaining 4, three are solar, where the observable is a phase-average that removes
   most of the error (§4c above).

A flat 1.53× tax on every adaptive call in the package, to fix errors that are unobservable in
one of its two main regimes and a notebook bug in another, is not a trade worth making by
default. As an opt-in for callers who need the oscillating probability itself, it is.

Its practical reach is narrower still than that: since the dispatch reorder the solar points are
answered by the hybrid strategy rather than the general path, and baseline scans can opt into
`cumulative=True`. What remains genuinely exposed is a raw `osc_prob` call, on a smooth varying
profile, where the oscillating probability is what the caller wants.

## 5. Still open

- **Two-consecutive-agreements in `osc_prob`'s adaptive loop.** Cost/benefit now derived in
  §4b: recommended, at 1.53× median on healthy calls, fixing 4 of 6 measured silent misses.
  Not yet built. The frozen-grid mode it does *not* fix needs a separate idea.
- **The energy axis** (decision doc §6.3). Adding `t_breakpoints` to the castle-wall
  energy scans cut them 4.7× (33.7 → 7.2 ms/pt), which is a real bite out of the
  77.4%, but it is profile-specific rather than the general improvement §6.3 has in
  mind. Extending `_osc_prob_scan_separable`'s detection to more profile shapes is
  untouched.
- ~~**PREM cells could take breakpoints too.**~~ **Done** (2026-08-03): notebook 03 cells
  85/90 and notebook 02 cells 78/83 now pass `earth.prem_layer_edges_along_chord`. Measured
  against `solve_ivp` at the grid those cells use, at no cost in time: 1.28e-04 → 6.78e-11,
  5.12e-04 → 1.42e-07, and 4.16e-03 → 3.30e-06 for the three directions — the through-the-core
  one having been outside the default tolerance without them.
- ~~**Proposal (3)**~~ **Done**: adopted in the notebooks, and `cumulative` now defaults to
  `'auto'`. See `DECISION_CUMULATIVE_DEFAULT.md`, which also records that the accompanying
  dispatch change makes the *hybrid* strategy the thing it must beat, not the per-point path.
- **`strict_convergence` is load-bearing beyond its own flag.** The cumulative scan's grid is
  sized by a probe that is always strict, because one early-stopping probe would misplace an
  entire scan rather than one point. Removing or weakening the flag would silently degrade every
  baseline scan in the package.
