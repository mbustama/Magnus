# Plan: exploiting palindromic density profiles in Magνs

**Written:** 2026-08-05, branch `dev-palindrome`, based on `main` at the PR #32 merge.

**Goal.** A chord through a spherically symmetric Earth meets every radius twice, so its density
profile reads the same from either end. NuOscProbExact exploits that to compose such a chord at
roughly two thirds of the cost (measured 1.4x-1.8x from about fifteen slabs upward). This plan
asks what the equivalent is for Magνs — and the answer is **not** the same algorithm.

---

## 0. The measurement that constrains the whole design

NuOscProbExact's saving rests on one fact: for a palindromic sequence of **constant-H** slabs,
`U_j = U_{n-1-j}` exactly, so each distinct operator is built once and used twice
(`fastkernels._slab_product_3nu_mirrored`). Two accumulators, `acc_a` growing on the right and
`acc_b` on the left, and the answer is `acc_a @ acc_b`.

**That identity does not hold in Magνs.** Measured on a perfectly symmetric profile
(`vcc(l) = vcc(L-l)` to machine precision), 12 slabs, `n_tpts_per_slab=41`, trapezoid:

| `magnus_exp_order` | `max_j |U_j - U_{n-1-j}|` | |
|---|---|---|
| 1 | 9.87e-13 | mirrors, to quadrature noise |
| **2** (the default) | **7.58e-03** | **broken** |
| 4 | **1.450** | unrelated operators |

**Why.** NuOscProbExact's slabs are constant in H, so a mirrored slab has *identical* inputs.
Magνs integrates H *across* each slab, and reversing an interval reverses the time-ordering
inside it. Term by term, with `s = a + b - t`:

* `Ω₁ = ∫H dt` is invariant — a plain integral does not care which way it is traversed;
* `Ω₂ = -(i/2)∫dt₁∫^{t₁}dt₂ [H(t₁),H(t₂)]` has its inner and outer limits exchanged, which
  swaps the commutator's arguments and **flips its sign**;
* **`Ω_k → (-1)^{k+1} Ω_k` holds for k <= 2 ONLY.** For k >= 3 it is false: reversing a nested
  commutator such as `[H₁,[H₂,H₃]]` gives `[H₃,[H₂,H₁]]`, which is not `±` the original. Measured
  in §3d — the sign rule leaves a 3.7e-08 residual at order 4, and re-deriving the terms from
  reversed samples is exact instead. *(An earlier draft of this plan asserted the sign rule "in
  general". That is wrong.)*

So the mirrored slab operator is not `exp(-iΩ)`. Copying the
NuOscProbExact composer would return a different answer at every order above 1 — silently, which
is the one thing an optimisation must never do, and which that project's own `palindromic`
docstring is emphatic about.

**This is the finding that makes the port a design problem rather than a translation.**

---

## 1. What the saving actually is here

Three distinct opportunities, in increasing order of value and difficulty.

### (A) Halve the Ω construction, not the exponentials — the general case

Compute the Magnus terms `Ω_k` once for each slab in the first half. The mirror slab's terms are
the same numbers with the even ones negated, so it needs **no Hamiltonian evaluations, no
quadrature and no commutator algebra** — only its own `expm`.

What that saves, per mirrored pair:

| stage | saved? |
|---|---|
| `H_func` evaluations (`n_tpts_per_slab` per slab, or 1-3 for `'gl'`) | **yes** |
| quadrature (`_cumulative_integral` / `_full_integral`) | **yes** |
| commutator algebra (`_omega_integrand`, nested chains — the dominant cost at order >= 4) | **yes** |
| matrix exponential (`_expm_stack`) | **no** — two different Ω, two exponentials |

This is worth most exactly where Magνs is slowest: a user `H_func` that is expensive to evaluate,
and high `magnus_exp_order` where the commutator chains dominate. It is worth least at order 2
with `'gl'`, where three Hamiltonian evaluations and one commutator are cheap against an `expm`.

**It must be measured before it is built.** The cost split between "build Ω" and "exponentiate Ω"
is the whole question, and it is not obvious.

### (B) The global transpose identity — cheaper, but conditional

If `H(l) = H(L-l)` **and** `H` is complex *symmetric* (`Hᵀ = H`, not merely Hermitian), then
reversing the time-ordering of a product is the same as transposing it, and

```
U(L, 0) = U(L/2, 0)ᵀ  @  U(L/2, 0)
```

which halves everything, exponentials included. **The condition is real**: for 3ν with a CP
phase `δ ≠ 0` the vacuum Hamiltonian is Hermitian but not symmetric, so this does *not* apply to
the general case. It does apply to 2ν, to 3ν with `δ = 0`, and to any case where the imaginary
part of `H` vanishes. Worth having as a special case precisely because Earth chords at `δ = 0`
are a common benchmark — but it must be gated on a property of `H`, not on the profile alone.

### (C) Reuse across a scan

Orthogonal to both, and possibly the largest win in practice: `osc_prob_earth` computes the same
chord at many energies. The *profile* is energy-independent, so the palindrome test is too. It
should be decided once per call, not once per point — the same reasoning that put
`_scan_for_hidden_features` at the entry point rather than inside an engine.

---

## 2. Where it applies in Magνs

| caller | palindromic? | notes |
|---|---|---|
| `osc_prob_earth`, `osc_prob_{2,3,4,5}nu_earth` | **yes, by construction** | a chord meets every radius twice; `earth.prem_layer_edges_along_chord` already returns edges symmetric about the midpoint |
| user profile via `osc_prob_matter_std_potential` | sometimes | must be detected, never assumed |
| solar / supernova profiles | **no** | monotonic or shock-structured; the test will correctly decline |
| `t_breakpoints` supplied by the user | depends | the test must consider the *slab edges actually used*, not the requested profile |

The Earth path is the target. It is also the path that already carries `t_breakpoints` from PREM,
so the slab grid is not free-form — which makes exact symmetry achievable rather than accidental.

---

## 3. Exactness: the part to copy verbatim

NuOscProbExact gets one thing exactly right and it should be carried over unchanged:

* **The test is exact equality, never a tolerance.** `np.array_equal`, not `np.allclose`. A
  tolerance would return a different answer for a nearly-symmetric profile depending on how
  nearly symmetric it is, which is a silent accuracy change keyed on an invisible property.
* **The producer is responsible for making it exactly symmetric**, not the consumer for tolerating
  near-symmetry. `earth._earth_slabs_cached` averages each width and density with its mirror —
  `w = (w + w[::-1])/2` — because floating-point addition is commutative, so the two ends of a
  pair come out bitwise identical. Their comment is worth quoting: *"This is not housekeeping."*
* **Two separate predicates**: "is this a palindrome" and "is the saving worth having at this
  size". NuOscProbExact keeps `palindromic()` and `worthwhile_mirror()` apart, with a per-flavour
  slab-count threshold, because the mirrored composer costs an extra matrix multiply per slab and
  only leads from about fifteen slabs up.

**And one thing to do differently:** their test is on the slab *inputs* (`h_stack`, `widths`).
Ours must be too — testing the *outputs* `U_j` would be both more expensive and, as §0 shows,
false at order >= 2 even when the profile is perfectly symmetric.

---

## 3b. MEASURED, 2026-08-05 — the numbers that set the criteria

Everything below is measured on this branch, on a real PREM chord (`costhz = -0.9`, 3nu,
`dCP = 3.70`) with a working prototype, not modelled. The prototype reproduces the shipped path
to **1e-15** wherever the slab inputs are exactly symmetric, so the timings are like-for-like.

**Where the time goes** (64 slabs). The exponential — the one thing the mirror cannot halve —
is 0.01-0.04 ms out of 0.8-31.6 ms. It is free. Everything else is halvable:

| method | order | total | `expm` |
|---|---|---|---|
| gl | 2 | 0.80 ms | 0.04 ms |
| trapezoid | 4 | 9.46 ms | 0.01 ms |
| trapezoid | 6 | 31.64 ms | 0.01 ms |

**Speed-up, cheap analytic Hamiltonian:**

| method | order 2 | order 4 |
|---|---|---|
| gl, 128 slabs | 1.38x | 1.39x |
| gl, 512 slabs | 1.09x | 1.23x |
| trapezoid, 256 slabs | 2.20x | 2.08x |
| simpson, 256 slabs | 2.77x | 2.09x |

**Speed-up against `H_func` cost** (gl, 128 slabs), which is the axis that decides it:

| spectral modes per evaluation | order 2 | order 4 |
|---|---|---|
| 0 (cheap analytic) | 1.19x | 1.28x |
| 100 | 1.28x | 1.45x |
| 400 | **1.61x** | **2.00x** |

The mirror halves every evaluation of `H_func`, so with `f` the share of slab time spent in
`H_func`, the speed-up is **`1/(1 - f/2)`** — 1.11x at `f = 0.2`, 1.43x at 0.6, 2.00x at 1.0.
That formula, not a slab count, is the thing to measure for any given Hamiltonian.

**The Earth chord is palindromic to rounding, and that is easily fixed.** On a uniform grid over
a `costhz = -0.9` chord, widths differ from their mirror by 3.9e-03 in natural units — which
sounds large but is **4.3e-15 relative**, i.e. floating-point noise, exactly what NuOscProbExact
describes ("about 1e-12 km on a 100 km slab"). `w = (w + w[::-1])/2` makes it exactly palindromic.
*(An earlier draft of this plan called that "a real asymmetry, not rounding". That was wrong: it
is rounding, and the fix is one line.)*

**But symmetrising the grid does NOT make the mirror exact for the quadrature methods.** Measured
directly — order 4, trapezoid, 64 slabs, before and after symmetrisation:

| order | plain grid | symmetrised grid |
|---|---|---|
| 2 | 2.01e-15 | 1.02e-15 |
| **4** | **3.72e-08** | **3.72e-08 — unchanged** |

So the 3.7e-08 is not grid asymmetry. It is the **discrete cumulative quadrature breaking the
sign rule**: `Ω_k → (-1)^{k+1} Ω_k` is exact for the continuum integrals, and for the closed-form
Gauss-Legendre scheme at every order, but the cumulative trapezoid/Simpson rule used for
`Ω_{k>=3}` integrates *from the left*, and that direction does not mirror. The identity survives
at order <= 2 and degrades at order >= 4.

**Consequences, and they differ by method:**

* **`gl` — exact at every order** (1e-15 in all configurations measured). Symmetrising the Earth
  grid is therefore all that is needed to unlock the saving on the Earth path, with no accuracy
  change beyond the 4.3e-15 grid nudge itself.
* **`trapezoid`/`simpson` — exact at order <= 2**, and at order >= 4 the mirror introduces an
  error of ~1e-8. Harmless against a 1e-3 tolerance, but it is an approximation *introduced by an
  optimisation*, which is a different thing from an approximation the caller asked for. Either
  restrict the mirror to order <= 2 on those methods, or make the cumulative integral direction-
  symmetric, or accept and document the 1e-8. **Not decided here.**

---

## 3c. Decision criteria, set in advance

**Ship the mechanism, gated on exact symmetry of the slab inputs actually used.** This is the
staging that removes almost all of the risk, and it was not obvious until §3b:

* a grid that is not exactly palindromic — which includes today's Earth chord — **takes the
  existing path, unchanged**. So the feature can ship with **zero bit-identity movement**, and
  `bitident.py` must show 0 of 11 moved. If it shows anything else, the gate is wrong.
* the beneficiaries are callers who supply their own symmetric profile or explicit
  `t_slab_edges`, where symmetry is the caller's responsibility and already exact.

**Symmetrising the Earth grid is a separate, later decision.** It is the only part that moves
bit-identity on the package's most-used path, and it is worth making only for someone whose Earth
workload is `H_func`-expensive. Decoupling it means the risky change is never bundled with the
useful one.

**Enable/disable threshold.** The mirror does strictly less work — half the `H_func` evaluations,
comparable Ω arithmetic — and was never slower in any configuration measured (worst 1.06x).
So unlike NuOscProbExact's `worthwhile_mirror`, no slab-count floor is obviously needed. That
must still be checked at small counts, where fixed overheads could invert it; if a floor is
needed it goes in as a measured constant with its population, like every other constant here.

**Worth building if**, on the caller's own Hamiltonian, the H-evaluation share `f` gives
`1/(1 - f/2) >= 1.3` — i.e. `f >= 0.46`. Below that the gain does not pay for a second code path
through a dispatch layer that produced three silent-wrongness defects in one session. A one-line
timing of `H_func` against the slab total settles it for any given user.

---

## 4. Proposed phases

Each phase ends in a decision, and any of them can end the work.

**Phase 1 — DONE, see §3b.** The cost split is measured and the answer is that `expm` is free,
so (A) is worth building wherever `H_func` is not trivially cheap. What is still unmeasured is
the small-slab-count end, where fixed overheads might need a floor.

**Phase 2 — the predicate.**
`magnus.palindromic(*arrays)` mirroring NuOscProbExact's semantics exactly, plus
`_slab_inputs_are_palindromic(h_samples, widths)`. Pure, cheap, exact. Tests: symmetric and
asymmetric profiles, odd and even slab counts, single slab, empty, and a **near**-symmetric
profile that must return `False`.

**Phase 3 — make the Earth chord exactly symmetric.**
Audit `earth.prem_layer_edges_along_chord` and the slab grid Magνs derives from it. If widths are
symmetric only to ~1e-12, apply the mirror-average. **This is a bit-identity change** and must be
justified per workload with `bitident.py`, as every such change in this repo has been.

**Phase 4 — the mirrored Ω path (opportunity A), if Phase 1 justifies it.**
In `magnus_expansion_multislab`: when the slab inputs are palindromic and the count is above
threshold, evaluate and integrate only the first half, form the mirror's terms by negating even
orders, and exponentiate the full stack as now. `_expm_stack` is already batched, so the slab axis
stays intact.

**Phase 5 — the transpose identity (opportunity B), gated on `Hᵀ = H`.**
Only if Phase 1 shows `expm` is a large share, since this is the only route that halves it.

**Phase 6 — acceptance.**
* bit-identity across the 11 workloads: **any** movement justified in writing;
* the physical-profile population re-run — the palindrome test must decline on solar, supernova
  and shock profiles, and a false positive there would be a correctness bug, not a missed
  optimisation;
* a symmetric-profile accuracy test against `solve_ivp` at orders 1, 2, 4 and 6;
* the timing re-measured on the same harness as Phase 1, so the claimed speed-up is the measured
  one.

---

## 5. Risks, in the order I would worry about them

1. **The order >= 2 sign flip (§0).** Already measured, and it is the whole reason this is not a
   port. Any implementation that skips the sign alternation will be wrong at the default order —
   and wrong *silently*, on the Earth path, which is the most-used path in the package.
2. **A near-symmetric profile taking the fast path.** Mitigated by exact equality, and by making
   the producer exact rather than the test tolerant.
3. **Bit-identity movement on the Earth workloads** (Phase 3). Expected, must be justified rather
   than absorbed.
4. **The saving may not exist at default settings.** Phase 1 exists to find that out before
   anything is built. `worthwhile_mirror`'s ~15-slab floor in NuOscProbExact is a warning: their
   composer trades a matrix multiply for half an expansion, and below that size the trade loses.
   Ours trades differently, so the crossover has to be measured, not inherited.
5. **`n_jobs > 1` and the batched slab axis.** The mirrored path changes which slabs are evaluated
   together; the parallel path must not silently diverge from the serial one. The existing
   invariant test (`n_jobs > 1` vs 1, asserted as exactly 0.0) will catch it.

---

## 6. What this plan does not propose

* **Porting `fastkernels`.** Magνs has no compiled kernels and this plan does not add any; the
  saving here is in avoiding work, not in doing it faster.
* **Touching the adiabatic/hybrid path.** It does not compose a fixed slab chain, so the
  palindrome has nothing to act on there. The target is the general Magnus ladder and the
  cumulative scan.
* **Any tolerance-based symmetry detection**, for the reason in §3.
