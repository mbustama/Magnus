# Plan: exploiting palindromic density profiles in Magνs

**Written:** 2026-08-05, branch `dev-palindrome`, based on `main` at the PR #32 merge.

**Goal.** A chord through a spherically symmetric Earth meets every radius twice, so its density
profile reads the same from either end. NuOscProbExact exploits that to compose such a chord at
roughly two thirds of the cost (measured 1.4x-1.8x from about fifteen slabs upward). This plan
asks what the equivalent is for Magνs — and the answer is **not** the same algorithm.

---

## 0. The measurement that constrains the whole design

> **CORRECTED — see §3d(i).** The claim below that the sign rule is false for k >= 3 is wrong;
> it is exact for every k, and the measured 3.7e-08 is a convergent O(h^2) quadrature error.
> The section is kept as written because the *rest* of it — that `U_j != U_{n-1-j}` at order
> >= 2, so the port does not transfer — is correct and is the finding that matters.

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

> **CORRECTED — see §3d(ii).** The gate described below cannot be implemented as specified:
> testing "the slab inputs actually used" requires evaluating them, which is the cost the
> feature exists to avoid. The `1/(1 - f/2)` speed-up applies only to the declare-it route.

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

## 3d. CORRECTIONS, 2026-08-05 (later session) — two claims above are wrong

Everything in this section was measured on `dev-palindrome` against the shipped code, on the
same PREM chord (`costhz = -0.9`, 3nu, `dCP = 3.70`). Scripts under the session scratchpad.

**Confirmed unchanged:** §0's `U_j != U_{n-1-j}` at order >= 2 (so the port really does not
transfer); §3b's grid palindromy (`3.9062e-03` absolute, **4.30e-15 relative**, `array_equal`
False, mirror-average fixes it); §1(B)'s transpose identity (at `dCP = 0`, `H - H^T` is exactly
zero and `U = F^T F` to 1.2e-15; at `dCP = 3.70`, `H` is asymmetric at **17% relative** and the
identity fails at 7.4e-01, while `U(d) = F(-d)^T F(d)` holds to 2.1e-15).

### (i) The sign rule is exact for every k, not only k <= 2

§0 and §3b assert that `Omega_k -> (-1)^{k+1} Omega_k` is algebraically false for k >= 3. **It is
not.** The argument given — that reversing `[H1,[H2,H3]]` yields `[H3,[H2,H1]]`, which is not
`+/-` the original — inspects one term of a commutator group in isolation. Summed over the whole
group the integrand is invariant: `Omega_3`'s integrand over the ordered simplex is
`([A1,[A2,A3]] + [[A1,A2],A3])`, which reversal maps to `([A3,[A2,A1]] + [[A3,A2],A1])`, and
antisymmetry gives `[A3,[A2,A1]] = [[A1,A2],A3]` and `[[A3,A2],A1] = [A1,[A2,A3]]` — the same
pair. There is also an operator-level derivation valid at all k: with `W(s) = U(b, a+b-s)` one
has `W' = W Atilde`, `W(0) = I`, `W(b-a) = U(b,a)`, whence `U(b,a) = exp(-Omega[-Atilde])`, and
since `Omega_k[-A] = (-1)^k Omega_k[A]` the rule follows for every k.

The **measurements in §3b are reproducible and correct** — 3.7366e-08 (trapezoid, order 4,
`n_tpts=41`) and 8.4360e-09 (simpson, order 4) — but 3.7e-08 is not an algebraic residual. It is
the error of the discrete cumulative quadrature, and it converges away:

| `n_tpts_per_slab` | 11 | 41 | 161 | 641 |
|---|---|---|---|---|
| sign-rule residual | 6.08e-07 | 3.74e-08 | 2.33e-09 | 1.45e-10 |

Exactly **4.0x per doubling** — O(h^2), the trapezoid's own rate. That rate holds on the PREM
chord, on a generic random complex `A(t)`, and at orders 2 through 8. Against `solve_ivp`,
sign-flipped terms converge to the *reversed* propagator through order 6 exactly as the unflipped
terms converge to the forward one; were the rule false at k >= 3, that sequence would stall at
order 3, and it does not. The residual also sits **4-5 orders below the error the same grid
already commits**: at `n_tpts=41`, order 4, it is 3.7e-08 against a quadrature error of 1.2e-03.

§3b's dead end (averaging the two sweep directions leaves the residual "exactly 3.72e-08 —
identical") is consistent with this rather than evidence against it: the cumulative trapezoid rule
is *already* direction-symmetric, so that patch was a no-op.

**Consequence.** The choice on `trapezoid`/`simpson` at order >= 4 is not "wrong versus exact". It
is "a 1e-8 discretisation difference for 2.2-2.8x" versus "exact for 1.02x". **Decision taken:
resample.** Reverse the retained samples and recompute, which is bitwise-exact at every order and
every method, saves only the `H_func` evaluations, and — importantly — removes the sign rule from
the shipped code path altogether, so none of the above needs to be relied upon in production.

### (ii) The gate in §3c cannot be implemented as specified

§3c says to gate on "exact symmetry of the slab inputs actually used" while evaluating only the
first half. Those two requirements are incompatible; you cannot test what you have not evaluated.
Measured:

* **Widths alone are not sufficient, and failing this is catastrophic.** A monotonic
  (solar-like) profile on a symmetrised uniform grid passes `np.array_equal(w, w[::-1])`, and the
  mirror is then wrong by **3.340e-01** — against 9.3e-16 for a genuinely symmetric profile on the
  same grid. This is the silent-wrongness mode the whole plan exists to avoid, and a widths-only
  gate walks straight into it on the *solar* family §2 predicts will "correctly decline".
* **The sampled `A` is never bitwise palindromic**, so an exact test on it never fires. Even with
  a perfectly symmetric profile on a symmetrised grid, `array_equal(At, At[::-1,::-1])` is False,
  differing by 4.9e-16 / 6.5e-16 / 5.7e-16 relative at gl orders 2 / 4 / 6 — because the shipped
  sampling computes the mirror slab's nodes as `(L-b) + h*s` and the forward slab's as
  `a + h*s`, which are different floating-point expressions for the same real number.

So the design has a genuine fork, and it is not the one §3c anticipated:

* **(a) Construct the mirror's sample positions exactly** (`t_mirror = L - t_forward`, one
  subtraction) and evaluate `A` at all of them, then test `array_equal`. Sound and exact — the
  test now means "is the profile symmetric", which is the property being exploited. But it
  forfeits the halved `H_func` evaluations, i.e. the entire `1/(1 - f/2)` saving that motivated
  the feature. What remains is only the halved Omega construction.
* **(b) Have the producer declare it** — the Earth entry points know by construction that a chord
  meets every radius twice. This delivers the saving, and it is what §3 already quotes
  NuOscProbExact as prescribing ("the producer is responsible for making it exactly symmetric,
  not the consumer for tolerating near-symmetry"). But it is a declaration, not a test, so a
  caller who declares wrongly gets the 3.34e-01 failure above.

**Not decided here.** (b) is the only route that delivers the advertised speed-up, and it is
the route with a silent-wrongness mode. That trade is the real decision this plan owes, and §3b's
`1/(1 - f/2)` figures should be read as belonging to (b) alone.

### (iii) The §2.7 prototype is wrong for odd slab counts

Transcribed verbatim, `gl_mirror` allocates `Om = np.empty((n,3,3))` and writes `Om[:m]` and
`Om[n-m:]` with `m = n//2`. For odd `n` those two slices skip index `m`, so the middle slab is
returned as uninitialised memory. Measured against the shipped path, with the worst slab being the
middle one in every case:

| `n_slabs` | 31 | 63 | 129 |
|---|---|---|---|
| max abs error | 7.13e-01 | 3.01e-01 | 1.48e-01 |

Even counts are unaffected (1e-15 at orders 2, 4 and 6), which is why the prototype looked sound.
Any implementation must handle the unpaired middle slab explicitly.

---

## 4. Proposed phases

Each phase ends in a decision, and any of them can end the work.

**Phase 1 — DONE, see §3b.** The cost split is measured and the answer is that `expm` is free,
so (A) is worth building wherever `H_func` is not trivially cheap. What is still unmeasured is
the small-slab-count end, where fixed overheads might need a floor.

> **Phases 2, 4 and 5 are superseded in part by §3d.** Phase 2's predicate is still wanted, but
> §3d(ii) shows what it can and cannot be applied to. Phase 4 must not "form the mirror's terms
> by negating even orders" — the decision recorded in §3d(i) is to resample instead, which is
> bitwise-exact and keeps the sign rule out of the shipped path. Phase 5 is dead for 3nu with
> `dCP != 0` per §1(B), confirmed in §3d.

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

> **Reordered by §3d.** The top risk is no longer the sign flip (which is exact, and which the
> resample decision removes from the shipped path anyway). It is **a caller, or an entry point,
> declaring a profile symmetric when it is not** — measured at 3.34e-01 on a monotonic profile
> that passes a widths-only gate. Every route that actually delivers the speed-up rests on such a
> declaration, so that is where the correctness argument has to be made.

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
