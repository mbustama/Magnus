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

## 5. Still open

- **Two-consecutive-agreements in `osc_prob`'s adaptive loop.** §3 above. Needs its
  own cost/benefit; ~1.6× on every adaptive call is the price.
- **The energy axis** (decision doc §6.3). Adding `t_breakpoints` to the castle-wall
  energy scans cut them 4.7× (33.7 → 7.2 ms/pt), which is a real bite out of the
  77.4%, but it is profile-specific rather than the general improvement §6.3 has in
  mind. Extending `_osc_prob_scan_separable`'s detection to more profile shapes is
  untouched.
- **PREM cells could take breakpoints too.** Notebook 03 cells 85/90 scan PREM
  per-point without `t_breakpoints`, though `earth.prem_layer_edges_along_chord`
  exists and `tests/test_oscprob.py:303` already shows it improves accuracy there.
  Left alone as out of scope for this brief.
- **Proposal (3)**, held per the decision doc §6.5 and unchanged by any of the above
  except that §3.1's cost split should be re-measured against the notebooks as they
  now stand.
