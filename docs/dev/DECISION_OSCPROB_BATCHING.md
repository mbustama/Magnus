# Decision: `osc_prob` baseline batching — do not build it (yet)

**Judges:** `docs/dev/HANDOVER_OSCPROB_BATCHING.md`
**Written:** 2026-08-02, against `main` at `155e01e`.

> **Status (executed 2026-08-02, branch `fix-adaptive-n-slabs-floor`).** §6 items 1,
> 2 and 4 are done; see the CHANGELOG "Fixed" entries and `NOTES_ADAPTIVE_REFINEMENT.md`
> for what the work changed about the diagnosis in §2 — the mechanism is not aliasing,
> and the "make `suggest_n_slabs` feature-scale-sensitive" alternative offered in §6.1
> is *not* the cheap one. §6.3 (energy axis) is measured but not built. §6.5 (proposal
> (3)) is deliberately still held.
**Machine:** 12 cores, `uptime` load 1.3–1.7 throughout (one other interactive
session idling; every number below is a min-of-2-or-3, so load noise is bounded).

---

## 1. The call

**Do not build proposal (2) — it already exists and buys nothing. Do not build
proposal (3) yet — it is correct and its speedup ratio is enormous, but it
addresses 17.6% of the only consumer's runtime while the untouched energy axis
is 77.4%. Before either, fix the accuracy bug that this investigation turned up:
the per-point baseline loops proposal (3) would replace are currently returning
answers wrong by up to 0.855 in probability.**

That last clause is the real result. The performance question is a wash; the
correctness finding is not.

---

## 2. The headline finding, which is not about performance

Notebook 03's castle-wall loop asks for `n_slabs=150`. Because `rtol`/`atol` are
left at their defaults, **`n_slabs` is ignored** (documented behaviour — see
`osc_prob`'s docstring) and the adaptive loop refines up from `min_n_slabs=1`.
On a profile with 50 density oscillations across the baseline it stops at
**`n_slabs=4`, `n_tpts_per_slab=2`** and declares convergence.

```python
# castle wall, 50 slabs of alternating density, E = 50 MeV, L = 5242.4 km
info = {}
oscprob.osc_prob(H_l, 0, 5242.4*gd.CONV_KM_TO_INV_EV, n_slabs=150,
                 n_tpts_per_slab=100, magnus_exp_order=3, convergence_info=info)
# -> info: {'n_slabs': 4, 'n_tpts_per_slab': 2}
```

Four slabs cannot resolve fifty density steps. The successive-iterate
convergence test passes anyway because the coarse grid *aliases* the periodic
profile: consecutive refinements agree with each other and disagree with reality.

| L [km] | truth (`solve_ivp` DOP853, rtol 1e-11) | notebook's call | `min_n_slabs=64` | `min_n_slabs=512` | `rtol=atol=1e-6` |
|---|---|---|---|---|---|
| 5242.4 | 0.986024 | **0.130993** | 0.985948 | 0.986285 | 0.986044 |
| 4780.4 | 0.365250 | **0.651263** | 0.360847 | 0.365278 | **0.651263** |

Three things to note:

1. The error is up to **0.855 in an oscillation probability**, not a tolerance-level
   discrepancy.
2. **Tightening the tolerance does not fix it.** At `rtol=atol=1e-6` the second
   point is still wrong, and wrong by exactly the same amount. Only raising
   `min_n_slabs` above the profile's feature scale fixes it. The criterion is
   under-sampling, and a stricter threshold on an aliased comparison is still aliased.
3. The user *asked* for 150 slabs, which would have been fine, and the library
   silently discarded that number in favour of 4.

Across the 600-point subsample of that scan, **7.2% of points differ by more than
1e-2** from the warm-started path, which is accurate everywhere I checked.

### How this nearly fooled me, and the lesson for proposal (3)

My first reference run used `n_slabs=4000, magnus_exp_order=4, rtol=atol=1e-10,
max_n_slabs=400000`. At L = 4780.4 km it returned 0.651291 — agreeing with the
notebook loop to *exactly* 0.00e+00. I nearly reported "the loop is fine here."

Both were wrong. The truth is 0.365. Two runs of the same aliasing scheme agree
with each other; that agreement is not evidence. Only `solve_ivp` with a
different integrator, plus a non-adaptive grid ladder shown to actually settle
(2000 → 128000 slabs, changes 9.0e-4 → 4.0e-5), separated them.

This is the third instance of the pattern already in the project's memory
(`hybrid_propagator`'s false certification; "agreement between two capped runs is
not convergence"). **It bears directly on proposal (3)'s validation plan**, which
proposes "agreement with the per-point path on non-commuting profiles" as a
primary oracle. That oracle is unsound: I have measured the per-point path being
wrong by 0.855 on exactly the non-commuting profile class the brief names. A
correct cumulative scan would *fail* that test.

The sound oracle is the one the project already uses on the energy axis
(`tests/test_oscprob.py:332`): pin `rtol=atol=None` so both paths use identical
grids, then demand agreement to 1e-12. That tests bookkeeping, which is what a
cumulative product can get wrong, and it does not pretend to test accuracy.

---

## 3. Measurements

Every script is in `~/Downloads/magnus-batching-bench/` (`b1`–`b11`, plus
`bench_setup.py` and the raw notebook timing in `nb03_times.json` / `nb03.log`),
deliberately outside the repo so nothing here needs committing. Run any of them
with `python3 b<N>_*.py` from that directory. The essential ones are inlined
below as well.

All start from array-capable `H_func`s copied verbatim from notebook 03, with
`ScalarHamiltonianWarning` promoted to an error so a scalar fallback cannot go
unnoticed (trap 2.2), and with `rtol=atol=None` pinned wherever integration
rather than refinement is being timed (trap 2.1).

### 3.1 Where notebook 03's 11 minutes actually goes

Executed the real notebook with per-cell timing (`nbclient`, `on_cell_executed`
hook), then classified every cell over 1 s by AST — looking for comprehensions
whose loop variable is a baseline vs. an energy and whose body calls `osc_prob`.

```
notebook 03 total: 652 s  (10.9 min -- confirms the ~11 min figure)

ENERGY  scanning via osc_prob   504.8 s  (77.4%)   cells 62, 76, 104, 48, 90
BASELINE scanning via osc_prob  115.0 s  (17.6%)   cells 57, 71, 99, 85, 43
plotting / other                 10.7 s  ( 1.6%)
everything under 1 s             13.3 s  ( 2.0%)
```

**Proposal (3) touches the 17.6% and nothing else.** The three most expensive
cells in the notebook (62, 76, 104 — 435 s between them) are energy scans at
fixed baseline. The brief's framing — "every baseline scan in the notebooks is a
Python loop" — is true but selects the minority of the cost.

*(A regex first told me baseline scanning was 0.4%; it broke on the `[nu_i][nu_f]`
indexing between the call and the `for`. Re-audited with an AST walk, per trap 5.
17.6% is the AST number.)*

### 3.2 The cumulative pass: measured floor, and it is real

A cumulative scan's traversal can be measured today without touching `src/`:
`t_slab_edges` pinned at the N requested baselines with `rtol=atol=None` performs
exactly the traversal, then discards the intermediates instead of recording them.

```python
edges = np.concatenate([[0.0], distances*gd.CONV_KM_TO_INV_EV])
pairs = np.stack([edges[:-1], edges[1:]], axis=1)
oscprob.osc_prob(H, 0, distances[-1]*gd.CONV_KM_TO_INV_EV, t_slab_edges=pairs,
                 n_tpts_per_slab=100, magnus_exp_order=3, rtol=None, atol=None)
```

| scan (6000 baselines) | per-point loop | cumulative floor |
|---|---|---|
| castle wall narrow | 15.8 s | **0.01 s** |
| castle wall wide | 23.6 s | **0.01 s** |
| PREM, 3000 pts × 3 directions | 2.1 / 1.9 / 2.6 s | 0.01 s each |

And it is genuinely correct, not just fast:

- constant `H`, 6000 segments vs. `expm`: **max |ΔP| = 1.5e-12** (the ordering
  test the brief asks for, and it passes)
- castle wall at L = 1e4 km: floor is 3.15e-3 from a tight reference; the
  notebook's own per-point answer is 2.95e-3 from it. **Same accuracy class.**

So the brief's asymptotic claim is vindicated: ~1500× on the traversal, and the
mathematics works. The problem is what it is 1500× *of*.

**Honest upper bound: 652 s → ~537 s. Notebook 03 goes from 10.9 min to 9.0 min.**
That is an upper bound in the strict sense — it credits proposal (3) with
reducing all baseline scanning to zero and charges it nothing for the
per-segment refinement it would actually need.

### 3.3 Proposal (2) already exists, and it delivers 0.92×

This is the flat contradiction of the brief. `osc_prob_energy_baseline` — public,
one layer above `osc_prob`, documented in `architecture.rst`, accepting **an
arbitrary `H_func` and an array of baselines** — is proposal (2), shipped.

```python
oscprob.osc_prob_energy_baseline(H_func, energy, L_array, magnus_exp_order=3)
# -> (npts, d, d);  H_func may be H(l), H(E), or H(E, l)
```

Measured on the castle-wall profile, 600 baselines:

```
osc_prob_energy_baseline(L=array) :    3.55 s
notebook-style manual loop        :    3.28 s
-> the existing API is 0.92x the manual loop
```

The brief estimates proposal (2) at "roughly 1.5–2×" from amortized entry cost.
An existing implementation of exactly that idea — accept the array, loop
internally — achieves **no speedup at all**. Entry cost is not where the time is.
Anyone who wants proposal (2) can have it today by calling one function up, and
they will find it does not help.

### 3.4 The 26× figure does not generalise

Both briefs lean on "wrappers batch, measured 26×". That number is specific to
vacuum and constant density, where `_osc_prob_scan_separable` engages. On a
spatially varying profile:

```
PREM scan, costhz=-1.0, 3000 baselines, E=20 MeV
  osc_prob_2nu_earth(L=array)      :  2.113 s
  same wrapper, per-point          :  3.345 s     -> 2x, not 26x
```

```
energy scan, 500 pts, exponential profile
  manual loop                      :  1.16 s
  osc_prob_energy_baseline(E=array):  0.57 s      -> 2.0x, not 26x
```

So the "wrappers already batch, the gap is only for custom `H_func`" framing in
the assignment brief is too generous to the status quo: **the wrappers do not
batch the baseline axis on a varying profile either.** The gap is wider than
stated — which is an argument *for* proposal (3), and I am recording it against
my own recommendation.

### 3.5 Two things that turned out not to matter

- **`n_jobs`.** The notebooks pass `n_jobs=10` on these loops. Measured 13.1 s
  (`n_jobs=10`) vs 14.1 s (`n_jobs=1`) for the same 6000-point scan — no
  meaningful difference. Not a lever.
- **Entry-path overhead.** Per 3.3, amortizing it buys nothing.

---

## 4. Cost/benefit, stated plainly

**Benefit:** at most 115 s per run of notebook 03 (18%), for the one consumer in
the repo. Similar in notebook 02. For an external user with a custom profile who
scans baselines, the ratio is large (1500× on the traversal) but the absolute
numbers are seconds — the largest scan in the repo costs 24 s today.

**Cost:** the mathematical core is genuinely small. `src/magnus/oscprob.py:2048`
already builds the full per-slab chain `U_chain` and collapses it with
`reduce(np.matmul, U_chain[::-1])`; the cumulative scan is that `reduce` replaced
by a running product that records intermediates. I want to be fair to the brief
here — it over-states the difficulty of the composition itself.

What is *not* small is everything around it, and the brief is right about all of it:
per-segment refinement, a tolerance contract for a product of N separately
converged factors, breakpoint placement under `logspace`, sort/unsort, and no
free composition with the energy-batched engine. Call it a week of careful work
plus a validation suite whose primary proposed oracle (§2 above) has to be
redesigned before it is written.

**The ratio:** roughly a week of design-heavy work on the package's single
highest-traffic entry point, to remove ~2 minutes from an 11-minute notebook that
runs when someone regenerates the docs. That is a bad trade *today*.

**Risk (§4.4):** the keyword-gated design does contain the blast radius, and the
identical-grid oracle would catch bookkeeping errors to 1e-12. I do not think a
transposed matrix product would ship. What would ship — because it already has,
in the per-point path — is a plausible-looking accuracy contract that is wrong on
oscillatory profiles. Proposal (3) would build a new tolerance story on top of a
convergence criterion that §2 shows is not currently delivering its nominal
tolerance. **Fixing the foundation is a prerequisite, not a follow-up.**

---

## 5. What would change my mind

In rough order of how likely each is to happen:

1. **A real user scans baselines with a custom profile at N ≫ 10⁴, or in a fit
   loop.** Everything above is about a 6000-point scan costing 24 s. At 10⁶
   points, or at 6000 points inside a minimizer called 10⁴ times, the 1500×
   stops being decorative. This is the condition I would actually watch for.
2. **The energy axis gets fixed first and baseline scanning becomes the
   majority of the cost.** If the 504.8 s of energy scanning drops to ~250 s
   (measured 2× is available today; see §6), baseline scanning goes from 17.6%
   to ~30% of a much shorter notebook. Proposal (3) gets proportionally more
   attractive as the thing it does not address shrinks.
3. **The accuracy fix lands and turns out to make the per-point path much more
   expensive.** This is the interesting one. Forcing `min_n_slabs` high enough to
   resolve the profile is the correct fix, and it will make every per-point call
   substantially slower — the honest cost of a correct answer. A correct
   per-point loop might cost 5–10× what today's wrong one does, at which point
   115 s becomes 600–1000 s and proposal (3) is suddenly the majority of the
   notebook. **I would re-run §3.1 immediately after the accuracy fix**, because
   this specific scenario could flip the answer on its own.
4. `osc_prob` acquiring array baselines for API-consistency reasons rather than
   performance ones — a legitimate motivation the brief does not claim.

What would *not* change my mind: a faster machine, or a better constant factor
on the per-point path.

---

## 6. What I would do instead, in order

1. **Fix the false convergence.** The smallest useful increment: make
   `suggest_n_slabs` (or `min_n_slabs`) sensitive to the profile's feature scale,
   or — much cheaper and arguably more honest — stop silently discarding a
   user-supplied `n_slabs` when tolerances are on. Treat it as a floor rather
   than ignoring it. The notebook asked for 150 and got 4; had 150 been honoured
   as a minimum, every number in §2 would have been right. This is a small change
   with a large correctness payoff, and it needs a test built on the
   `solve_ivp` oracle, not on path-vs-path agreement.
2. **Re-run the notebooks and check the figures.** 7.2% of the castle-wall scan
   is materially wrong. Notebooks 02 and 03 both ship affected figures. This is
   a docs-correctness issue, not a performance one.
3. **Then look at the energy axis, not the baseline axis.** It is 77.4% of the
   runtime, a 2× is already reachable through `osc_prob_energy_baseline` with no
   new algorithm, and the separable engine that gets 26× on constant density is
   already written — extending its detection to more profile shapes is a smaller,
   better-understood job than proposal (3).
4. **Delete the stale performance comments.** Eight sites across notebooks 02 and
   03 (cells 9/14/22/27 and 8/17/29/34) assert that "the ~4.5 ms is almost
   entirely fixed entry-path cost … n_slabs=1 and n_slabs=150 both measure 4.5 ms
   per call, which is how you can tell the physics is not what costs." That is
   the reasoning `HANDOVER_OSCPROB_BATCHING.md` §1.2 identifies as wrong, and
   §3.3 above independently refutes it. It is currently on `main`, in the
   teaching material, presented as a lesson.
5. **Revisit proposal (3) after (1)–(3)**, with §3.1 re-measured.

If proposal (3) is built anyway, the smallest useful increment is: cumulative
recording behind a keyword, sorted input only, `rtol=atol=None` only (no adaptive
story at all in v1), validated by the identical-grid/1e-12 oracle and against
`expm` for constant `H`. That version is genuinely small — it is the `reduce` at
`oscprob.py:2058` plus a sort — and it would let the notebooks opt in explicitly
while the accuracy contract is designed separately. Everything expensive about
the proposal lives in the adaptive story; ship it without one first.

---

## 7. Where I disagree with the brief being judged

| `HANDOVER_OSCPROB_BATCHING.md` says | this document finds |
|---|---|
| (2) is worth "roughly 1.5–2×" | The same idea is already shipped as `osc_prob_energy_baseline` and measures **0.92×** |
| (3) is "the difference between minutes and seconds" | True per-scan (15.8 s → 0.01 s); **18% of the notebook**, because 77% is energy scanning |
| "Do (3), and let (2) fall out of it" | Do neither yet; the pullback in its own final paragraph is the more accurate sentence |
| Validate by "agreement with the per-point path on non-commuting profiles" | Unsound — the per-point path is **wrong by 0.855** on that exact profile class |
| Composition is settled, the work is bookkeeping and accuracy control | Agreed, and the bookkeeping is smaller than implied (one `reduce`); accuracy control is the whole job |

And where the assignment brief itself is too generous to the status quo:
"wrappers batch both axes / anyone using the wrappers already batches" is not
true on a varying profile — `osc_prob_2nu_earth(L=array)` gets 2×, not 26×.
