# Handover brief: batching `osc_prob` over baselines

**Repo:** `/home/mbustamante/Research/magnus` (GitHub `mbustama/Magnus`, currently private)
**Written:** 2026-08-01, from a session that profiled the notebooks and found where the time actually goes.
**Status of this work:** *not started*. This brief exists so it can be picked up cold.

This is a **performance and API exploration**, deliberately separated from the
notebook work that produced it. Nothing here is required for correctness. Read
section 1 before deciding whether to do any of it — one of the two proposals is
much more valuable than the other, and the cheap win has already been taken.

---

## 0. The one-paragraph version

`osc_prob` takes a single baseline. Every baseline scan in the notebooks is
therefore a Python loop that re-integrates from `L = 0` for each point. Two
things could change: **(2)** let `osc_prob` accept an array of baselines, which
amortizes the fixed per-call cost; and **(3)** compute a whole baseline scan in
*one* pass, by partitioning `[0, L_max]` at the requested points and taking
cumulative products of the per-segment evolution operators. (2) is a constant
factor. (3) is asymptotic, and is the one worth doing.

---

## 1. What has already been done, and what the measurements actually say

**Do not redo this.** A prior session measured all of it on this machine
(12 cores, load < 1 when timing).

### 1.1 The cheap win is already taken

The single largest factor was never `osc_prob` — it was whether the
user's `H_func` accepts an array of positions. `magnus.magnus._evaluate_A`
probes for this and falls back to a Python loop if it fails.

| `H_func` | per `osc_prob` call |
|---|---|
| scalar-only (how every notebook was written) | 7.80 ms |
| array-capable, same physics, identical output | **1.70 ms** |

That 4.6× has been taken: the notebooks now use array-capable Hamiltonians,
`ScalarHamiltonianWarning` fires once per session on the fallback, and
`README.md` plus `docs/source/methodology.rst`
(`.. _array-capable-hamiltonians:`) document how to write one.

**Any benchmark for this work must start from an array-capable `H_func`**, or it
will measure the thing that is already fixed.

### 1.2 A correction worth carrying, because it is easy to re-derive wrongly

An earlier reading of the data claimed *"the cost is fixed entry-path overhead,
not the physics"*, on the evidence that `n_slabs=1` and `n_slabs=150` both timed
at 4.5 ms. **That conclusion was wrong**, and the reasoning behind the error is
instructive:

```
n_slabs=1,   default tolerances     192 H evals    20.4 ms
n_slabs=150, default tolerances     192 H evals    10.2 ms
n_slabs=1,   rtol=atol=None           9 H evals     0.96 ms
n_slabs=150, rtol=atol=None         158 H evals     8.11 ms
```

`n_slabs` is a *starting point for the adaptive refinement, not a cap*. With
tolerances left at their defaults the loop refines regardless, so both cases
converged to the same 192 evaluations — which is why they timed alike. It was
never evidence that integration is free.

**When benchmarking anything here, pin `rtol=atol=None`** or you are timing the
refinement loop's convergence, not the thing you changed.

### 1.3 Where the time goes (cProfile, 200 calls, scalar `H_func`)

```
_evaluate_A                      63% cumulative
  the user's H_func              38,400 calls for 200 osc_prob calls (192 each)
  hamiltonian_3nu_matter         0.437 s
  np.diag                        0.344 s   (38,400 calls)
_expm_stack / eigh               0.40 s / 0.28 s
suggest_n_slabs                  0.19 s
probe_eval_mode                  0.08 s
```

With an array-capable `H_func` the first block collapses and the profile
re-centres on `eigh` and the slab products — **re-profile before optimizing;
this table is the "before" picture.**

### 1.4 What already batches, and what does not

- `osc_prob_{2,3,4,5}nu_*` wrappers **accept arrays for `energy` and `L`** and
  return `(npts, d, d)`. Measured 0.174 ms/point against 4.54 ms/point for the
  loop — 26×. Notebooks 02 and 03 now use this for their vacuum and
  constant-density sections.
- `_osc_prob_scan_separable` batches the **energy** axis when many energies
  share one baseline, for the separable form `H(E,l) = H_E(E) + V_CC(l)·M`.
  See `docs/source/methodology.rst`, "Silent vectorization and the
  energy-batched scan engine".
- **Nothing batches the baseline axis for a general user `H_func`.** That is
  the gap this brief is about.
- `oscprobstd.osc_prob_3nu_vacuum_std` does **not** accept arrays (it reshapes
  to `(3,3)`); `osc_prob_2nu_vacuum_std` and `osc_prob_2nu_matter_std` do. Do
  not assume uniformity across `oscprobstd`.

---

## 2. Proposal (2): accept an array of baselines in `osc_prob`

### What

`osc_prob(H_func, t_ini, t_fin, ...)` currently requires scalar `t_fin`; an
array raises `ValueError: The truth value of an array with more than one
element is ambiguous`. Make it accept an array and return `(npts, d, d)`,
matching the wrappers' existing convention.

### Why it is only worth so much

The naive implementation — loop internally over `t_fin` — moves the loop inside
the library and saves only the repeated entry cost (validation, probing,
dispatch). Measured entry cost is ~2.4 ms of a ~6.9 ms call with default
validation, so expect **roughly 1.5–2×**, not 26×. The 26× the wrappers get
comes from the *separable* batching in 1.4, which a general `H_func` cannot use.

If you implement (2) naively, say so in the docstring, or users will reasonably
expect wrapper-like speedups and not get them.

### Design notes

- **Return shape.** Follow the wrappers: scalar in, `(d, d)` out; array in,
  `(npts, d, d)` out. `np.ndim(t_fin) == 0` is the discriminator used elsewhere.
- **`t_ini` stays scalar.** Every use in this codebase scans the far end.
- **Sort, compute, unsort.** Do not assume the input is monotonic; the
  cumulative scheme in (3) requires sorted baselines and users will pass
  `np.logspace` (fine) and hand-built lists (not fine).
- **The refinement loop is per-point today.** Batching means deciding whether
  convergence is judged per point or for the whole array. Per point is more
  faithful but complicates the loop; whole-array is simpler but makes one
  hard point dictate the cost for all. The energy-batched engine already
  solves exactly this with per-energy convergence masking — read
  `_osc_prob_scan_separable` before inventing a third answer.
- **Validation.** `validate_input` currently checks a scalar `t_fin`. Extend
  rather than skip: `tests/test_validation.py` has the conventions
  (`ValueError` naming the offending parameter, message via
  `gd.ERROR_MSG_NO_COLOR`).

---

## 3. Proposal (3): one cumulative pass for a whole baseline scan

**This is the one worth doing.**

### The idea

A baseline scan asks for `U(0, L_1), U(0, L_2), …, U(0, L_N)` with
`L_1 < L_2 < … < L_N`. These are nested. Instead of `N` independent
integrations from zero, partition `[0, L_N]` at the requested points and
accumulate:

```
U(0, L_k) = U(L_{k-1}, L_k) · U(0, L_{k-1})
```

Each segment is integrated once, and the running product is recorded at every
requested baseline. Total work is one traversal of `[0, L_N]` rather than `N` of
them.

### Why it is asymptotically better

Present cost is `sum_k (slabs needed for [0, L_k])`, which for a scan with
comparable slab density everywhere is `O(N · n_slabs)` and in the worst case
`O(N²)` in the number of requested points. The cumulative scheme is
`O(N + total_slabs)` — the segments partition the interval, so their slab counts
add up to what a *single* `L_N` calculation already costs.

**Expected win on the notebooks' remaining loops:** notebook 02's PREM section
is 3,000 points at 100 slabs; the cumulative version integrates roughly the same
total interval once. That is the difference between minutes and seconds, not a
constant factor.

### This is exactly the composition the engine already relies on

`U(l_c, l_a) = U(l_c, l_b) · U(l_b, l_a)` is exact, and the multi-slab chain is
already built on it (`compute_evolution_operator_multiple_slabs`). A prior
session confirmed the same principle when stitching adiabatic and Magnus
segments — see the seventh/eighth-pass notes in the project memory: *"Composition
of pieces from different methods needs no special smoothing."* So the
mathematics is settled; the work is bookkeeping and accuracy control.

### The hard parts, in the order they will bite

1. **Accuracy is not composable the way the operators are.** Each segment
   carries its own truncation error, and the product accumulates them. `N`
   segments each converged to `rtol` does **not** give a product converged to
   `rtol`. Decide the contract: per-segment tolerance `rtol/N`? A global
   estimate from a halved-partition comparison? Whatever you choose, state it
   in the docstring — the existing tolerance semantics are already "difference
   between successive refinements", not a global bound, and this makes that
   looser still.

2. **Refinement must be per segment.** A segment crossing the Earth's core
   needs more slabs than one in the crust. Refining the whole partition
   uniformly throws away most of the benefit. The existing
   `t_breakpoints`/PREM-layer-aligned machinery in `oscprob.py` is the closest
   precedent.

3. **Requested baselines are not good segment boundaries.** With `np.logspace`
   the first segments are tiny and the last enormous. You almost certainly want
   an internal slab grid chosen on physics grounds, with the requested `L_k`
   inserted as mandatory breakpoints — which is precisely what
   `t_breakpoints` already does for PREM layers. Reuse it.

4. **Unitarity is free but is not accuracy.** The product of unitary factors is
   unitary regardless of how wrong it is, so unitarity checks will pass while
   the answer drifts. **Do not use unitarity as the convergence signal** — this
   codebase has been caught by a self-consistent-but-wrong convergence flag
   before (see the project memory's note on `hybrid_propagator`'s false
   certification, where agreement between two identical runs passed a check).

5. **It does not compose with the energy-batched scan for free.** That engine
   batches energies at fixed baseline; this batches baselines at fixed energy.
   Combining them is a two-axis problem — scope it out of v1.

### Validation plan (do not skip)

- **Exactness first, no quadrature.** Constant `H`: the cumulative product must
  equal `expm(-i·H·L_k)` for every `k`, to machine precision. This isolates the
  bookkeeping from every numerical question and is the test that catches an
  ordering bug. `tests/test_magnus_expansion.py` has a two-constant-slab
  ordering test to model it on.
- **Against the current path.** Same `H_func`, same tolerances, per-point
  `osc_prob` versus the cumulative scan, over: a smooth exponential profile, an
  asymmetric profile with `dCP != 0` (non-commuting slabs — this is where
  ordering errors show), and a full PREM crossing.
- **Against `solve_ivp`** at a handful of baselines, tolerance 1e-10, as the
  independent arbiter. Note `solve_ivp` gets very slow at low energy — see the
  memory's sixth-pass note — so choose points, do not sweep.
- **Reversed and unsorted input**, to prove the sort/unsort round-trips.
- **Prove each new test bites** by breaking the thing it protects and watching
  it fail, and assert the mutation actually applied. A no-op `str.replace`
  looks exactly like a working guard; this project has been bitten.

---

## 4. Working notes for whoever picks this up

- **Repo conventions:** branch off `main`, never commit to it; numpydoc
  docstrings with `r"""`, `.. versionadded:: 1.0.0` on new public functions and
  **no** `versionchanged` before 1.0.0 ships; update `CHANGELOG.md` (one
  `[1.0.0rc1]` section, no `[Unreleased]`); why-focused commit prose.
- **Checks:** `python3 -m pytest tests/ -q` (~4 min, 500+ tests);
  `--cov` enforces a 90% floor; `python3 -m ruff check src/magnus/ tests/ docs/`
  (rule set pinned on purpose — do not widen); `python3 docs/check_doc_snippets.py --rst-only`;
  fast docs build recipe in the notebook handover brief (keep autoapi **on**).
- **A static pre-flight exists and is worth reusing.** The notebook session
  wrote a script that walks every notebook cell's AST and binds each
  `magnus.*` call against the installed signature. Discovering an API mismatch
  by running a notebook costs minutes and surfaces only the first one. If you
  change `osc_prob`'s signature, run that before re-executing anything.
- **Benchmark hygiene:** check `uptime` before quoting a number (a concurrent
  job on this machine once inflated a wall clock by 60%), pin
  `rtol=atol=None` when timing integration, and start from an array-capable
  `H_func`.
- **Persistent notes:** `/home/mbustamante/.claude/projects/-home-mbustamante-Research-magnus/memory/`
  — twenty-four passes of history, including several traps in this exact area.

## 5. Recommendation

Do **(3)**, and let **(2)** fall out of it: the array-baseline API is the natural
interface to a cumulative scan, and implementing (2) on its own buys ~2× for
work that (3) would then largely replace. Build it behind a keyword
(`cumulative=True` or a `strategy`-style flag) so the per-point path stays as
the reference implementation the tests compare against — that is how the
adiabatic hybrid strategy was introduced, and it worked well.

If time is short, the honest answer is that **nothing here is needed**: the
4.6× from array-capable Hamiltonians is already banked, and the wrappers already
batch the cases most users hit. This is an optimization for people who supply
their own `H_func` and scan baselines — real, but a narrower audience than the
effort might suggest.
