# Decision: the cumulative baseline scan (proposal (3)) — revisited after the accuracy fix

**Supersedes the verdict in:** `DECISION_OSCPROB_BATCHING.md` §1 and §4.
**Written:** 2026-08-02, against `main` at `bdf3074` (i.e. *after* the `n_slabs` floor and the
notebook `t_breakpoints` landed).
**Machine:** the same 12-core box the earlier documents used; every timing is a min-of-2-or-3.
**Prototype:** ~40 lines, in the scratchpad, deliberately outside the repo. Nothing in `src/`
was touched to produce any number here.

---

## 1. The call

**Build it — but for accuracy, not for speed.**

The performance case on its own is weak and would not justify the work: the honest
end-to-end effect is **−21% on notebook 03 and −24% on notebook 02**, and it is capped
there by Amdahl, because 72% of both notebooks is energy scanning that a cumulative
baseline scan does not touch.

The case that changed is accuracy. On every profile measured, the cumulative scan is
**2 to 7 orders of magnitude more accurate than the per-point path at 140–585× less cost**.
`DECISION_OSCPROB_BATCHING.md` §4.4 listed accuracy as proposal (3)'s principal *risk*;
measured, it is its principal *benefit*. For a package whose recent history contains two
false-convergence bugs, a baseline scan that is structurally more accurate and two orders
of magnitude cheaper is worth more than the 21%.

**The one blocking condition:** the tolerance story is mandatory, not deferrable. §4 below
has a measured counterexample where the cumulative scan is silently wrong by 1.6e-2. The
"ship v1 with no adaptive story" shortcut that `DECISION_OSCPROB_BATCHING.md` §6 proposed
is **not** safe.

---

## 2. What the strategy is

For a fixed Hamiltonian `H(l)` and N requested baselines `L_1 < L_2 < ... < L_N`, the
per-point path computes N independent answers, each re-traversing the profile from the
origin:

```
P(L_1):  [====]
P(L_2):  [=========]
P(L_3):  [==============]              N traversals, total work ~ sum_i n_i slabs
```

But the evolution operator is a time-ordered product, so each answer is a *prefix* of the
next: `U(0->L_2) = U(L_1->L_2) @ U(0->L_1)`. One traversal that records its running
product yields all N:

```
         [==|====|=====|====|===...]   1 traversal, total work ~ (n_acc + N) slabs
            ^     ^     ^     ^        record here
```

Concretely: build one slab grid, compute every slab's evolution operator in a single
batched call, then accumulate `U <- U_k @ U` left to right and snapshot `U` at each
requested baseline. That is exactly the `reduce(np.matmul, U_chain[::-1])` already at
`oscprob.py:2058`, with its intermediates recorded instead of discarded — which is why the
mathematical core is small.

**The grid is the whole design.** It must be the union of three things:

```
grid = (the N output baselines) u (a uniform accuracy grid of n_acc slabs) u (breakpoints)
```

- the **output baselines**, so every requested answer lands exactly on a slab edge and is
  read off rather than interpolated;
- a **uniform accuracy grid**, because the output baselines are typically logspaced and so
  are dense where accuracy is cheap and sparse where it is expensive — the opposite of what
  is needed;
- the **breakpoints**, for the usual reason (a slab straddling a density discontinuity
  degrades the quadrature regardless of `magnus_exp_order`).

`DECISION_OSCPROB_BATCHING.md` §3.2 measured a "floor" of 0.01 s by pinning `t_slab_edges`
at exactly the N output baselines — that is the union with `n_acc = 0`. It is fast, and on
that particular scan it happened to be accurate, but it ties resolution to *how many points
the user asked for* rather than to what accuracy requires. It is not a floor for a correct
answer, and §4 shows where it breaks.

---

## 3. Measurements

### 3.1 Where the time goes now

Re-derived on current `main` from the post-fix notebook runs, classified by AST walk (a
regex got this wrong once already):

| | notebook 02 | notebook 03 |
|---|---|---|
| total | 394.3 s | 445.0 s |
| **energy scanning** | 280.9 s (71.2%) | 320.9 s (72.1%) |
| **baseline scanning** | 96.2 s (24.4%) | 95.9 s (21.5%) |
| everything else | 17.2 s | 28.2 s |

Baseline scanning is a larger *fraction* than the 17.6% recorded before the accuracy fix,
but a smaller *absolute* cost (115.0 s → 95.9 s), because the notebooks got faster.

### 3.2 Per-cell, against a `solve_ivp` oracle

Every row is the same configuration on both sides, timed the same way. The error column is
against DOP853 at `rtol=1e-11`, sampled across the scan — never per-point against
cumulative, which is the unsound comparison this project has been burned by.

| cell | profile | per-point | cumulative | speedup | error: per-point → cumulative |
|---|---|---|---|---|---|
| 03/57 | castle wall, narrow | 17.9 s | 0.034 s | **525×** | 8.7e-4 → 1.1e-10 |
| 03/57 | castle wall, wide | 19.0 s | 0.032 s | **585×** | 1.6e-3 → 1.3e-10 |
| 03/71 | noisy, high-amplitude | 23.5 s | 0.048 s | **489×** | 1.1e-5 → 7.4e-9 |
| 03/71 | noisy, low-amplitude | 23.9 s | 0.056 s | **430×** | 1.2e-4 → 3.0e-8 |
| 03/99 | solar exponential | 29.9 s | 0.108 s | **278×** | 1.8e-5 → 2.8e-5 |
| 03/85 | PREM, 2 directions | 7.3 s | 0.02–0.05 s | **140–380×** | 2.8e-7 → 1.4e-10 |
| 03/43 | exponential / Gaussian | 1.8 s | 0.007 s | **172–297×** | equal or far better |

The castle-wall rows reach 1e-10 for a reason specific to that profile: it is piecewise
*constant*, so once breakpoints sit on the walls each slab has a constant `H` and the Magnus
expansion is exact inside it. The smooth and PREM rows are the representative ones.

### 3.3 How the speedup scales with N

Castle wall, `n_acc = 2000`, everything else fixed:

| N | 20 | 100 | 500 | 2000 | 6000 | 20000 |
|---|---|---|---|---|---|---|
| speedup | 6× | 21× | 93× | 261× | 387× | 334× |

It plateaus at **300–400×**, not the 1500× quoted in `DECISION_OSCPROB_BATCHING.md` §3.2.
The asymptote is the mean slab count per point: per-point costs `sum_i n_i ~ N * mean(n)`,
cumulative costs `n_acc + N`, so the ratio tends to `mean(n)`, a few hundred. Below
N ≈ 500 the win is under 100×; below N ≈ 100 it is ~20×.

### 3.4 Correctness of the bookkeeping

Three separate oracles, because a cumulative product's characteristic failure is
transposition and off-by-one, not inaccuracy:

- **Constant `H`, 500 outputs over 699 slabs, against `expm`:** max |ΔP| = **2.05e-13**.
  This is the ordering test the original brief asked for.
- **Identical-grid oracle** (the sound one, per `tests/test_oscprob.py:332`): for 118
  output baselines, `osc_prob` run on the *exact* grid prefix with `rtol=atol=None` versus
  the cumulative snapshot: max |ΔP| = **1.49e-14**.
- **Unitarity** of every cumulative output: **5.9e-13**.

### 3.5 The comparison is fair

The obvious objection is that the per-point baseline is handicapped — the notebooks call
`osc_prob` in a raw loop, which gets no warm starts, while `osc_prob_energy_baseline` seeds
each point from the previous point's converged grid. Measured on the castle wall at
N = 2000:

```
raw per-point loop                    5.24 s   err 3.0e-4
warm-started (osc_prob_energy_baseline)  5.50 s   err 3.5e-4     -> warm starts buy 0.95x
cumulative                            0.0201 s  err 8.2e-11
```

Warm starts buy nothing, consistent with the 0.92× that
`DECISION_OSCPROB_BATCHING.md` §3.3 measured for the same API. The 260× stands against the
better of the two per-point paths.

---

## 4. The blocking condition: `n_acc` cannot be guessed

The solar scan is the counterexample, and it is why "ship v1 with no adaptive story" is not
an option:

| solar exponential, 1000 baselines | time | error vs `solve_ivp` |
|---|---|---|
| per-point (what the notebook does) | 29.9 s | 1.8e-5 |
| cumulative, `n_acc = 2000` | 0.016 s | **1.65e-2** ← silently wrong |
| cumulative, `n_acc = 20000` | 0.108 s | 2.8e-5 |

At `n_acc = 2000` the answer is wrong by 1.6e-2 and nothing warns. The reason is measurable
directly — the per-point path's own converged slab count along that trajectory:

```
per-point at   1% R_sun : converges at n_slabs = 300
per-point at  30% R_sun : converges at n_slabs = 4770
per-point at 100% R_sun : converges at n_slabs = 14883
```

So the rule is **`n_acc` must be at least what the per-point path needs at the longest
baseline**, and 2996 total slabs is not. This is the same failure class as the
`n_slabs`-floor bug: a plausible-looking answer, exactly unitary, no warning. Shipping the
cumulative scan without a criterion that discovers `n_acc` would reintroduce it on a new
code path.

Note also what this implies for the criterion's *design*: successive refinement in `n_acc`
is the natural candidate, but the lesson recorded in `NOTES_ADAPTIVE_REFINEMENT.md` §1 —
that a ladder can thrash rather than converge, and that `np.allclose` fires on coincidences
— applies here unchanged, so agreement should be required twice in a row.

---

## 5. The energy axis: measured, and there is no comparable win

Since the energy axis is 72% of both notebooks, it is the obvious question. Two things were
measured and one was verified structurally.

**The existing batched API buys almost nothing.** `osc_prob_energy_baseline` on the
castle-wall energy scan (400 points, fixed L): **1.24×**.

**Even the separable engine, called directly with a perfect decomposition, buys 2.3×.**

```
castle-wall energy scan, 400 pts, fixed L = 1e4 km
  per-point loop              3.01 s   (7.5 ms/pt)
  _osc_prob_scan_separable    1.32 s   (3.3 ms/pt)   -> 2.3x   (agreement 2.7e-4)
```

This is the structural asymmetry between the two axes, and it is worth stating plainly:

> **Baselines nest; energies do not.** `P(0->L_1)` is a *prefix* of `P(0->L_2)`, so N
> baselines share one traversal. `P(E_1)` shares nothing with `P(E_2)` — each energy needs
> its own propagation through the whole profile. Separability saves re-evaluating `VCC(l)`
> across energies and improves batching, which is a constant factor; it removes no work.

The 26× recorded elsewhere for the separable engine is specific to **constant density**,
where the profile needs exactly one slab, so there is nothing to scale.

**The real energy-axis lead is elsewhere, and is not this proposal.** `osc_prob` — the entry
point the notebooks loop over — never dispatches to `_osc_prob_ip_exp_dispatch` or
`_osc_prob_hybrid_dispatch`; both are reachable only from the wrapper layer (`oscprob.py`
lines 4088/4098, 4426/4433, 4780/4785, 7819). Notebook 03 cell 104 — at 117.4 s the single
most expensive cell in either notebook — is a **solar** energy scan written as a raw
`osc_prob` loop, so it runs the general Magnus path while the package holds two dedicated
fast paths for exactly that case. **This is structurally confirmed but its magnitude is
unmeasured**, and it should be measured before anyone relies on it.

---

## 6. Cost, stated plainly

**Small, and demonstrated:** the traversal, the running product, the snapshotting and the
grid union are ~40 lines, written and validated to 1.5e-14 (§3.4).

**Not small:** the `n_acc` criterion (§4), which is the whole job, as
`DECISION_OSCPROB_BATCHING.md` §4 said. Beyond it: sort/unsort of user input, validation,
the scalar-`H` fallback (the prototype assumes an array-capable Hamiltonian), and plumbing
for `convergence_info`, logging and `n_jobs`.

One defect in the prototype worth carrying forward: it stores a matrix at *every* grid edge,
not just the N output points. Harmless at these sizes, but at large `n_acc` it wastes memory
for no reason — a shipped version should snapshot only at the output indices.

**No composition with the energy axis.** A cumulative scan is per-energy by construction.
That is a documented limitation, not a defect.

---

## 7. What would change this call

- **The `n_acc` criterion turns out to be as hard as the per-point one.** §4's ladder could
  thrash the same way §1 of `NOTES_ADAPTIVE_REFINEMENT.md` describes. If a sound criterion
  cannot be found in about a day of work, the accuracy argument that carries this decision
  collapses, and the answer reverts to "not for 21%".
- **The solar/hybrid lead in §5 is measured and is large.** If routing the notebooks'
  energy scans through the wrapper layer takes a 117 s cell to single digits, that is a
  bigger prize than this proposal on a bigger share of the runtime, and it should go first.
  It needs no new algorithm — only reaching code that already exists.
