# Decision: try the hybrid strategy before the interaction-picture fast path

**Written:** 2026-08-03, executing Task 1 of `HANDOVER_DISPATCH_AND_ADOPTION.md`.
**Machine:** the same 12-core box the earlier documents used; every timing is min-of-2-or-3.
**Status: SHIPPED**, as an ordering change at the three call sites in `oscprob.py`
(`osc_prob_matter_std_potential`, `osc_prob_matter_nsi`, `osc_prob_liv`). Nothing was removed.

---

## 1. The call

**Reorder.** `_osc_prob_hybrid_dispatch` now runs before `_osc_prob_ip_exp_dispatch`.

The handover brief framed this as a performance change, and the performance case is large
(median **397×**). But the stronger reason is that **the old order contradicted the package's
own published contract**, in three places at once:

> ``'auto'`` tries the hybrid strategy first, under the same conditions, but falls back
> silently to the ``'magnus'`` strategies above

— the `strategy` docstring of all three wrappers, and `docs/source/adiabatic_strategy.rst`,
which both explicitly list the interaction-picture integrator as one of the `'magnus'`
strategies. The code did the reverse. A second, quieter consequence: a user who passed
`strategy='hybrid'` on an exponential profile did **not** get the hybrid strategy, because
the fast path answered first and `strategy` was never consulted.

So this is a conformance fix that happens to be worth 397×, rather than an optimisation that
happens to be safe.

## 2. What was measured

50 (energy, baseline) configurations at the default `rtol=atol=1e-3`, spanning the standard,
NSI and LIV families and 0.5–100 MeV over 0.028–1.0 R_sun. Each dispatch path was isolated by
monkeypatching the other two to `NotImplemented`; a spy recorded whether the path under test
*answered* or *declined*, so a decline is never mistaken for a fast answer from the fallback.

**Accuracy is scored against `solve_ivp`/DOP853 at `rtol=1e-12, atol=1e-14` and nothing else.**
The reference was checked for convergence by tightening to `rtol=1e-13` (it moves by 5e-11) and
for unitarity (8e-11). No Magnus path was ever compared against another Magnus path.

| | hybrid | `ip_exp` (as shipped) |
|---|---|---|
| certified | **50 / 50** | 22 / 50 |
| warnings emitted | **0 / 50** | — |
| worst error | 1.8e-04, against a requested 1e-3 | — |
| speedup where both answer | **28×–594×**, median **397×** | — |
| cost of its 28 declines | — | 369 s total, mean **13.2 s** each |

### 2.1 The energy range, which is where the brief expected trouble

The brief's principal worry was that hybrid's certification would fail below 10 MeV — "exactly
where `ip_exp` was originally built to help" — and that the win would shrink there. **The
opposite is true.** 2nu standard potential, L = R_sun:

| E [MeV] | hybrid | `ip_exp` | general |
|---|---|---|---|
| 0.5 | 0.030 s, 3.8e-06 | **declines**, 13.3 s | 0.094 s, 5.4e-04 |
| 1 | 0.031 s, 4.0e-06 † | **declines**, 11.7 s | 0.470 s, 7.8e-05 |
| 2 | 0.031 s, 2.3e-05 | **declines**, 12.7 s | 0.182 s, 1.8e-05 |
| 3 | 0.028 s, 3.2e-06 | **declines**, 12.2 s | 0.128 s, 3.5e-06 |
| 5 | 0.027 s, 7.8e-06 | **declines**, 12.1 s | 0.125 s, 2.3e-05 |
| 8 | 0.028 s, 2.3e-06 | **declines**, 12.8 s | 0.170 s, 1.6e-06 |
| 10 | 0.033 s, 2.2e-07 | 12.4 s, 6.3e-05 | 0.031 s, **1.7e-02** |
| 20 | 0.026 s, 8.8e-06 | 11.6 s, 1.1e-04 | 0.124 s, 3.2e-06 |
| 40 | 0.029 s, 5.4e-06 | 11.7 s, **4.6e-06** | 0.111 s, 2.9e-06 |
| 100 | 0.031 s, 7.3e-06 | 12.9 s, **2.4e-07** | 0.034 s, 5.7e-04 |

† the R_sun/1 MeV row was measured while other work contended for the CPU (0.218 s); the value
quoted is the uncontended 0.9 R_sun/1 MeV row, which agrees with every other hybrid row.

Below 10 MeV the fast path does not certify **at all**, and charges ~12 s for the privilege
before falling through. Above it, the fast path certifies but is 379–457× slower.

### 2.2 The sibling families

The brief warned not to assume the 2nu result transfers, citing the neglected-`Omega_2`
coefficient jumping three orders from 2 to 3 flavors. Measured across NSI (`eps_aa`/`eps_ab` at
0.05/0.02 and 0.5/0.2) and LIV (`n_liv` = 0 and 1, eigenvalues at 10% and 100% of the vacuum
splitting at 10 MeV), 25 configurations: **hybrid certified 25/25**, `ip_exp` 10/25, and the
pattern is identical — `ip_exp` declines below ~10 MeV, and is 380–430× slower where it answers.

### 2.3 Tolerance: the fast path has a narrower window than anyone thought

Every measurement above is at the default 1e-3. Tightening it, at L = R_sun:

| tolerance | 1e-3 | 1e-5 | 1e-6 | 1e-8 |
|---|---|---|---|---|
| `ip_exp` certifies (10/40/100 MeV) | **yes** | no | no | no |
| hybrid certifies | yes | yes | yes | yes |

`ip_exp` certifies **only at the default tolerance, and only above ~10 MeV**. At anything
tighter its early-refusal bound (Fix B of `BUG_IP_EXP_MEMORY.md`) fires and it declines in
0.14–0.34 s — cheaply and correctly, exactly as that document predicted for tighter tolerances.
Hybrid reaches 1.6e-08 at 40 MeV and 1.6e-09 at 100 MeV when asked for 1e-8.

## 3. What the reorder does *not* do, contrary to the brief

**It does not speed up notebooks 02 or 03.** The brief's stated reason to sequence this before
notebook adoption was that notebook 03 cell 104 (117.4 s, the most expensive cell in either
notebook) is a solar energy scan. It is — but it is written as a raw `oscprob.osc_prob` loop on
a *fixed* grid (`n_slabs=200, n_tpts_per_slab=100, magnus_exp_order=2`, no tolerance), and
`osc_prob` never reaches either dispatcher. An AST walk over every function in `oscprob.py`
confirms only three functions call both: `osc_prob_matter_std_potential`, `osc_prob_matter_nsi`
and `osc_prob_liv`. Across all twelve notebooks, the only call reaching the reordered code is
notebook 12 cell 5 (`osc_prob_2nu_sun`, 18 MeV, L = 4 l_scale), which goes from 12.089 s to
0.026 s — **465×**, on a cell that is not a bottleneck.

Capturing cell 104 requires *rewriting* it onto the wrapper layer, which also swaps an
author-chosen fixed grid for an adaptive one. That is a notebook-content decision with its own
accuracy question, not a consequence of this change. **Task 2 therefore does not depend on this
one**, and the brief's sequencing rationale does not hold.

## 4. Cost, and what argues against

**`ip_exp` is more accurate at 40–100 MeV**, in 8 of the 22 points where it answers: 2.4e-07
against hybrid's 7.3e-06 at 100 MeV is a real 30×. Every one of those points is at the default
tolerance, where both are already 100× inside what was asked for, and it costs 400× to collect.
For solar neutrinos this is also above where the physics lives.

**Hybrid over-claims at tight tolerances.** In 4 of 12 tolerance-sweep rows it reported
*certified* while missing the requested tolerance — 2× over at 40 MeV/1e-8, up to **26× over at
10 MeV/1e-8** (2.6e-07 for a requested 1e-8). Its self-certification compares successive
refinement levels, not truth, so it can converge to its own floor and call that success. This is
not an argument for the old order — `ip_exp` declines at every one of those rows — but it is a
real limitation, and it is the natural next thing to measure on this path.

**Hybrid has a ~0.026 s floor per point** that the general path does not. On short, easy
baselines the general path answers the same query in 0.002–0.006 s. This is pre-existing to
`strategy='auto'` and unchanged by the reorder; it only becomes visible now that hybrid is
reached more often.

## 5. A defect found while measuring, not fixed here

The **general** path returned an answer outside the tolerance it was given, without saying so,
in **9 of the 50 configurations** (18%): worst 2.0e-02 against a requested 1e-3 (LIV, 40 MeV),
and 1.7e-02 on the plain 2nu solar case at 10 MeV. In every instance the only warning raised is
`MagnusConvergenceWarning` — which also fires on rows accurate to 1.6e-06, so it does not
discriminate. `ToleranceNotAchievedWarning` never fires: the ladder reaches two successive
agreements and declares victory.

This is the false-convergence signature of `NOTES_ADAPTIVE_REFINEMENT.md` §1, on a third code
path, at the *default* tolerance. It is recorded here rather than fixed because it is
independent of dispatch order — but note that the reorder now *masks* these points (hybrid
answers 10 MeV/R_sun at 2.2e-07 instead of 1.7e-02), which is a reason to have written it down.

It also bears on `HANDOVER_DISPATCH_AND_ADOPTION.md` Task 3, which recommends *against* building
the two-consecutive-agreements safeguard. That recommendation predates these measurements; the
cost/benefit should be re-derived rather than inherited.

## 6. Why the change is safe in shape

Both dispatchers already return `NotImplemented` when they decline, and the caller falls
through. Putting hybrid first cannot lose an answer: where hybrid does not certify, `ip_exp`
runs exactly as before. `strategy='magnus'` makes the hybrid dispatcher return
`NotImplemented` without doing any work, so that route to the fast path is bit-for-bit
unchanged — pinned by `test_magnus_strategy_still_reaches_the_interaction_picture_fast_path`.

The ordering itself is pinned by
`test_hybrid_strategy_precedes_the_interaction_picture_fast_path`, which was **mutation-checked**:
with the old order restored in `osc_prob_matter_std_potential` alone, it fails
(`assert None is True`). `test_sun_2nu_fast_path_matches_solve_ivp` would otherwise have kept
passing while silently testing the hybrid strategy instead of the fast path it is named for, so
it now disables hybrid explicitly and asserts, via a spy, that the fast path is what answered.

## 7. Deliberately not done

- **Narrowing or retiring the `ip_exp` gate.** §2.3 makes the case discussable — the fast path
  certifies only at the default tolerance and only above ~10 MeV — but the reorder already
  removes almost all of the harm, because `ip_exp`'s ~13 s is now paid only when hybrid has
  already declined, which did not happen in 50 configurations. Retiring it would also give up
  the fallback that makes this change safe. Separate decision, separate evidence, per the
  project's `integration_method='gl'` precedent.
- **The batched multi-energy axis.** `_osc_prob_ip_exp_dispatch` is genuinely batched over
  `(nE, n_slabs, d, d)`; `_hybrid_propagator_scan` is a per-point Python loop. Every measurement
  here is single-point. `BUG_IP_EXP_MEMORY.md` makes the batched fast path look unattractive
  anyway (~1.3 GB and ~10 s per energy), but it is unmeasured against hybrid.
- **Antineutrinos.** All measurements are `nubar=False`.
