# Handover: dispatch order, notebook adoption, and what to leave alone

> **Task 1: DONE** (2026-08-03), see `DECISION_DISPATCH_ORDER.md`. Two things this brief
> assumed turned out differently:
>
> * **The sequencing rationale does not hold.** "Task 1 before Task 2" rests on the reorder
>   speeding up notebook 03 cell 104. It cannot: cell 104 is a raw `osc_prob` loop on a fixed
>   grid, and `osc_prob` reaches neither dispatcher. Across all twelve notebooks only notebook
>   12 cell 5 touches the reordered code. **Task 2 does not depend on Task 1.**
> * **Hybrid does not weaken at low energy**, which was this brief's main worry. It certified
>   50/50 across the standard, NSI and LIV families, 0.5-100 MeV; it is `ip_exp` that declines
>   below ~10 MeV, 28 times out of 50, at a mean of 13.2 s each.
>
> Task 3's recommendation ("I recommend against doing it next") should be **re-derived, not
> inherited**: measuring for Task 1 found the general path silently exceeding its requested
> tolerance in 9 of 50 configurations, worst 2.0e-02 against 1e-3, which is new evidence for
> exactly the safeguard Task 3 describes. See `DECISION_DISPATCH_ORDER.md` §5.

**Written:** 2026-08-03, at the close of the session that produced PR #21 and PR #22.
**Starting point:** `main` once PR #22 merges. Read `docs/dev/BUG_IP_EXP_MEMORY.md` and
`docs/dev/DECISION_OSCPROB_CUMULATIVE.md` first; this brief assumes both.
**Machine for every number quoted:** 12-core box, min-of-2-or-3 timings.

---

## Where things stand

| | state |
|---|---|
| PR #21 | merged — `n_slabs` is a floor; notebooks 02/03 figures corrected |
| PR #22 | open — `ip_exp` memory fix, early refusal, output guard, cumulative scan |
| tests | 659 passing, ruff clean, docs clean under `-W` |
| notebooks | 02 (394 s), 03 (445 s), 07 regenerated; **not** yet using `cumulative=True` |

**Not code, but blocking:** GitHub Pages is disabled on the repository, so the
"Documentation Deployment" workflow has failed on every commit since well before this work
(`Get Pages site failed... verify that the repository has Pages enabled`). It needs
Settings → Pages → source "GitHub Actions". Nothing in the codebase can fix it.

---

## Task 1 (do first): reorder the dispatch so hybrid precedes `ip_exp`

### The finding

`_osc_prob_ip_exp_dispatch` is tried before `_osc_prob_hybrid_dispatch`, unconditionally and
regardless of `strategy`. On solar configurations it is **400–500× slower for the same
answer**. Measured, 2nu, L = R_sun, default tolerance, one energy per call, error against
`solve_ivp` (DOP853, rtol 1e-10):

| E [MeV] | `ip_exp` (as shipped) | hybrid | general per-point |
|---|---|---|---|
| 100 | 9.47 s, 1.7e-06 | **0.02 s, 1.5e-05** | 0.04 s, 7.1e-04 |
| 40 | 9.35 s, 9.4e-06 | **0.02 s, 9.3e-06** | 0.08 s, 5.0e-06 |
| 10 | 11.20 s, 7.4e-06 | **0.03 s, 7.4e-06** | 0.16 s, 5.2e-06 |

At 40 and 10 MeV hybrid matches `ip_exp` to two significant figures. Because `ip_exp` runs
first, hybrid never gets the chance.

Reproduce with the script sketched in "Measurement recipes" below — it monkeypatches the two
dispatchers to `NotImplemented` to isolate each path.

### Why the change is low-risk in shape

The dispatchers already return `NotImplemented` when they decline, and the caller falls
through. Putting hybrid first therefore cannot lose an answer: if hybrid does not certify,
`ip_exp` runs exactly as it does today. The change is an ordering, not a removal.

### What must be validated before it ships

The evidence above is **three energies, one profile, 2nu**. Do not ship on it. Needed:

1. **Down the energy range**, especially below 10 MeV and into the 1 MeV region where solar
   physics actually lives. `hybrid`'s certification is least certain there, and that is
   exactly where `ip_exp` was originally built to help. If hybrid declines below some
   energy, the reorder still works (fallback), but the win shrinks — quantify where.
2. **Across the sibling families.** `_osc_prob_ip_exp_dispatch` is reached from
   `osc_prob_matter_std_potential`, `osc_prob_matter_nsi` and `osc_prob_liv`
   (`oscprob.py` ~4088/4098, ~4426/4433, ~4780/4785). NSI and LIV add terms to `H_E`; the
   project's history records the neglected-`Omega_2` coefficient jumping three orders from 2
   to 3 flavors, so do not assume the 2nu result transfers.
3. **Accuracy, against `solve_ivp` only.** Never `ip_exp` against hybrid — agreement between
   two of this package's own paths is the trap that produced two false-convergence bugs
   already. `solve_ivp` is affordable above ~20 MeV and prohibitive below ~5 MeV (the project
   once measured 523 s for a single 2 MeV reference); plan the grid around that.
4. **The existing tests that pin current dispatch.** `test_sun_2nu_fast_path_matches_solve_ivp`
   and `test_generic_osc_prob_sun_hybrid_strategy_resolves_hard_case` both encode assumptions
   about which path answers. Expect to update them deliberately, not to silence them.

### Expected payoff

Notebook 03 cell 104 — at 117.4 s the most expensive cell in either notebook — is a solar
energy scan. If hybrid answers it at ~0.02 s/point instead of ~0.12 s/point, that cell drops
substantially. **This is the reason to do this task before Task 2.**

---

## Task 2: re-measure the notebooks, then decide adoption and the default flip together

Deliberately sequenced *after* the reorder, because the reorder moves the numbers this
decision rests on.

1. **Re-derive the cost split** (`docs/dev/DECISION_OSCPROB_CUMULATIVE.md` §3.1 has the AST
   recipe — use an AST walk, not a regex; a regex got this wrong once already). Current
   split: energy scanning 72%, baseline scanning 21.5% (nb03) / 24.4% (nb02).
2. **Adopt `cumulative=True`** in the baseline-scan cells where it still matters. Expected
   −21% (nb03) and −24% (nb02) *before* the reorder; recompute after.
3. **Then decide the default flip** for `cumulative`. Evidence is in
   `DECISION_OSCPROB_CUMULATIVE.md` §8: 124× faster and 11× more accurate on a solar scan,
   195× and seven orders on the castle wall. The argument against is blast radius — it moves
   every baseline scan in the package. The project's precedent is `integration_method='gl'`:
   implement, measure, flip in a separate decision with its own evidence.
4. **Fold in PREM breakpoints** while there: notebook 03 cells 85/90 scan PREM without
   `t_breakpoints`, though `earth.prem_layer_edges_along_chord` exists and
   `tests/test_oscprob.py::test_prem_breakpoints_improve_accuracy` already shows it helps.

Re-run notebooks with `nbclient` and per-cell timing; compare stored figures by hashing the
embedded PNGs, and check every changed figure against `solve_ivp` rather than assuming a
finer grid is better — error is **not** monotone in `n_slabs` (see
`NOTES_ADAPTIVE_REFINEMENT.md` §1).

---

## Task 3 (lowest priority): two consecutive agreements in `osc_prob`'s ladder

Would close the case the `n_slabs` floor does not cover — a caller who states no feature
scale still gets a ladder that can fire on a coincidental agreement
(`NOTES_ADAPTIVE_REFINEMENT.md` §1 and §3). `_osc_prob_ip_exp_core` already uses exactly this
safeguard.

**I recommend against doing it next.** It costs ~1.6× on *every* adaptive call in the package
(one extra refinement level, each ~1.5× the last), for a benefit that only materialises on
profiles whose feature scale the phase seed underestimates — and the cumulative scan now
sidesteps that case for baseline scans entirely. Keep it recorded; spend elsewhere.

---

## Measurement recipes and traps

**Isolating a dispatch path.** Monkeypatch the dispatcher to a function returning
`NotImplemented`, restore afterwards:

```python
real_ip = op._osc_prob_ip_exp_dispatch
op._osc_prob_ip_exp_dispatch = lambda *a, **k: NotImplemented
```

**Always run measurements under a memory cap.** `(ulimit -v 6000000; python3 ...)`. This
session OOM-killed the desktop application three times before the cause was found; the cap
turns a machine-wide kill into a clean `MemoryError` naming the shape.

**Traps that cost real time here, all of them repeatable:**

- A **module constant consumed as a default argument** binds at import and cannot be
  monkeypatched. A test that varies it passes trivially, comparing two identical runs. Read
  such constants at call time (`max_entries=None`, resolve inside).
- **Check that a passing test can fail.** After finding the above, the fixed test was checked
  by mutation: correct block order gives exactly 0, reversed gives O(1).
- **`gh` is snap-confined** and has a private `/tmp` — it cannot read a `--body-file` under
  `/tmp`. Pipe through `--body-file -`.
- **`pgrep -f <script>` matches the waiting shell's own command line**, so
  `until ! pgrep -f foo.py; do sleep; done` never exits. Use `run_in_background` on the real
  command instead. Several orphaned waiters were leaked this way.
- **Piping through `tail` buffers everything** until EOF, defeating `python3 -u`. Write to a
  file and read it.
- `/tmp` is cleared on reboot; a scratchpad there does not survive one.
- **`solve_ivp` without `t_eval` stores `y` at every accepted step.** Harmless at high energy
  (~20k steps), not at low.

**The oracles, and which is which.** `solve_ivp` (a genuinely different integrator) is the
only accuracy oracle. Identical-grid comparison with `rtol=atol=None` is the *bookkeeping*
oracle — it tests ordering and indexing, not accuracy. `expm` for constant `H` tests time
ordering. Never certify one configuration of the Magnus path with another; that specific
mistake produced the 0.855 probability error in `DECISION_OSCPROB_BATCHING.md` §2.

---

## Deliberately not done, and why

- **The "real cure" for `ip_exp`** — a second interaction picture absorbing the diagonal
  matter phase exactly — destroys the closed-form slab integral the method rests on
  (`∫V` for an exponential profile is itself exponential, giving
  `exp(i(Δ_jk l + c(1-e^{-l/λ})))`, not elementary). That is research, and it may be moot: if
  the reorder shows hybrid dominating, the right change is to **narrow the `ip_exp` gate**,
  not to make it cleverer.
- **Energy-axis batching.** Measured and closed: the existing batched API buys 1.24×, and the
  separable engine called directly with a perfect decomposition buys 2.3×. Baselines nest and
  energies do not — `P(E1)` shares nothing with `P(E2)`. There is no cumulative-scan
  counterpart for energies, and looking for one is wasted effort.
