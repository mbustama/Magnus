# Handover: session log 2026-08-04/05, and the palindromic-profile investigation

**Written:** 2026-08-05, at the close of a long session. **Read §§0-2 before touching anything.**

**Where things are:** `main` at the PR #32 merge (`6362bf3`), all CI green. Working branch
**`dev-palindrome`**, based on that, containing **one uncommitted file**:
`docs/dev/PLAN_PALINDROMIC_PROFILES.md`. Nothing else is uncommitted.

---

## 0. Verify the base before starting

```bash
git -C ~/Research/magnus branch --show-current      # dev-palindrome
git -C ~/Research/magnus log --oneline -1           # 6362bf3 Merge pull request #32
git -C ~/Research/magnus status --porcelain         # only PLAN_PALINDROMIC_PROFILES.md
python -c "import sys; sys.path.insert(0,'src'); import magnus.oscprob as o, magnus.adiabatic as a; \
  print(o.HYBRID_YIELDS_TO_CUMULATIVE_MIN_POINTS, o.CUMULATIVE_AUTO_MIN_POINTS, \
        hasattr(a,'oscillation_sampling'))"
# must print: 8 2 True
```

If `oscillation_sampling` is missing you are on a pre-#28 tree and none of §1 applies.

**Suite baseline:** 758 tests, all passing. `ruff check src/ tests/` clean. Docs build clean under
`-W --keep-going`. GitHub Actions **works now** (the repo went public mid-session, which restored
free Actions minutes) and GitHub Pages **is enabled and deploying**.

---

## 1. What landed this session — five merged PRs

The session executed `HANDOVER_PHYSICAL_PROFILES.md` and then three follow-ups. Full write-up is
`FINDINGS_ROBUSTNESS_PROGRAMME.md` §13 (§§13.1-13.21). Summary:

| PR | branch | what |
|---|---|---|
| #28 | `dev-robustness` | the physically-motivated profile population + findings; two shipped bug fixes; `oscillation_sampling` |
| #29 | `routing-regression-test` | hybrid no longer stands aside for a *disabled* engine |
| #30 | `example-notebooks` | notebooks 13 (tabulated solar) and 14 (supernova shock) |
| #31 | `docs-reconcile` | reconciled stale documents; documented the turbulence blind spot |
| #32 | `data-provenance` | `SOURCES.md` for the one third-party data file |

All five branches deleted, local and remote. `origin/dev`,
`origin/notebook-breakpoints-and-cumulative` and `origin/reorder-hybrid-before-ip-exp` were
**kept** — the last has a commit that exists **only on that laptop** ("ahead 1"), so do not
delete it without pushing first.

### 1.1 The three runtime changes (everything else was docs/tests/diagnostics)

* **Missing scalar squeeze.** `osc_prob_3nu_earth(E, ..., cumulative=True)` with scalar energy and
  baseline returned `(1, d, d)` instead of `(d, d)`, so `P[nu_i][nu_f]` silently indexed the wrong
  array. Pre-existing on `main`.
* **`convergence_info` + `cumulative=True`** raised `TypeError` instead of returning a
  probability. Pre-existing on `main`.
* **The hybrid path stood aside for a disabled engine.** With `cumulative=False` above the seam it
  declined into the chain and landed on the general ladder: **1.157e-05 -> 2.966e-03**, a factor of
  256, on one extra baseline. Fixed. Semantics clarified: `strategy='magnus'` is the exact route to
  pre-1.0.0 numbers at every N; `cumulative=False` guarantees only that the cumulative scan is not
  used.

### 1.2 The scientific result, and the methodological lesson behind it

**Whether an error is observable depends on the family, and the same test decides it.** Averaging
an instantaneous scan over ~6 oscillation lengths and comparing against `solve_ivp`:

| configuration | instantaneous | averaged | reduction |
|---|---|---|---|
| BS05 solar model, d = 2, 5 MeV | 1.380e-03 | **2.603e-05** | 53x |
| SN turbulence C\* = 0.1, 45 MeV | 1.701e-03 | **1.565e-04** | 23x |
| SN shock, 70 km front | 1.095e-03 | 9.773e-04 | 2x |
| SN shock, 0.07 km front | 2.033e-01 | **2.135e-01** | 3x |

Errors that shrink by >20x are **phase** and no observable sees them; errors that barely move are
**envelope** — a shock changes the adiabaticity of the level crossing, which moves the conversion
probability itself. **Of everything measured, nothing reaches a user silently.**

**Three claims made during the session were retracted on measurement**, and the retractions are
kept in the findings rather than edited away:

1. the solar "silent miss" headline (it is phase, not envelope);
2. "the code handles the solar model" (the d = 2 failures cluster *entirely* at 5 MeV);
3. "`t_breakpoints` cures the shock" (true on scans; on single points it improved 7 of 18
   configurations and **worsened 11**).

**The lesson to carry forward:** state, per family, *which observable* a pass criterion applies to
**before** measuring. The brief's P1 ("inside 1e-3 or warns") was executed exactly as written and
measured a quantity solar physics never observes.

---

## 2. THE LIVE WORK: palindromic profiles

**Goal.** A chord through a spherically symmetric Earth meets every radius twice, so its density
profile reads the same from either end. NuOscProbExact
(`/home/mbustamante/Research/NuOscProb/NuOscProbExact`) exploits that to compose such a chord at
~2/3 cost (`src/fastkernels.py`: `palindromic`, `worthwhile_mirror`,
`_slab_product_3nu_mirrored`). The user asked for the Magνs equivalent.

**Status: investigated in depth, NOT implemented. Plan written, uncommitted.** Every number below
was measured on `dev-palindrome` with working prototypes, on a real PREM chord
(`costhz = -0.9`, 3nu, `dCP = 3.70`, 2 GeV).

### 2.1 It is not a port — the algorithm does not transfer

NuOscProbExact's saving rests on `U_j = U_{n-1-j}` for palindromic slabs. That holds because
**their slabs have constant H**, so a mirrored slab has bitwise identical inputs. Magνs integrates
H *across* each slab. Measured, on a perfectly symmetric profile:

| `magnus_exp_order` | `max_j |U_j - U_{n-1-j}|` |
|---|---|
| 1 | 9.87e-13 (mirrors) |
| **2 — the default** | **7.58e-03 (broken)** |
| 4 | **1.450 (unrelated)** |

**A direct port would be silently wrong at the default order, on the Earth path.**

### 2.2 The three opportunities, and which survive

**(A) Halve the Ω construction.** Evaluate H and build the Magnus terms for the first half only;
derive the mirror's Ω from them. Saves H evaluations, quadrature and commutator algebra; does
**not** save the exponential. **Survives, and is the one to build.**

**(B) The global transpose identity** `U(L,0) = U(L/2,0)ᵀ U(L/2,0)` halves everything including
the exponential — but requires **`Hᵀ = H`**, i.e. *symmetric*, not merely Hermitian. **Dead for
this user's work**, measured:

| | `max\|H - Hᵀ\|` | `max\|U - FᵀF\|` |
|---|---|---|
| δ = 0 | 0 exactly | 3.4e-15 ✓ |
| δ = 3.70 | 1.4e-11 | **1.3e-01** ✗ |

It *does* generalise — `U(δ) = F(-δ)ᵀ F(δ)`, verified to **6.3e-15** — but that needs the half
evolution at both `+δ` and `-δ`, which is two half-computations, i.e. one full computation. **No
saving.** Losing (B) costs almost nothing, because the exponential it uniquely halves is
0.01-0.04 ms out of 0.8-31.6 ms — free.

**(C) Decide the palindrome once per call, not per energy.** The profile is energy-independent.
Same reasoning that put `_scan_for_hidden_features` at the entry point.

### 2.3 THE CRITICAL CORRECTION — the sign rule is only valid for k <= 2

An intermediate conclusion in this session was that `Ω_k → (-1)^{k+1} Ω_k` under interval
reversal, in general. **That is wrong for k >= 3**, and the plan has been corrected. Reversing a
nested commutator `[H₁,[H₂,H₃]]` gives `[H₃,[H₂,H₁]]`, which is not `±` the original.

Measured — two ways of building the mirror's Ω for the quadrature methods, on an exactly
symmetrised grid, 64 slabs:

| method | order | **sign rule** | **reverse samples & recompute** | speedup (resample) |
|---|---|---|---|---|
| trapezoid | 2 | 1.02e-15 ✓ | 1.11e-15 ✓ | 1.41x |
| trapezoid | 4 | **3.72e-08** ✗ | **8.95e-16** ✓ | 1.02x |
| trapezoid | 6 | **3.74e-08** ✗ | 1.42e-15 ✓ | 1.02x |
| simpson | 4 | 8.39e-09 ✗ | 1.00e-15 ✓ | 1.03x |
| simpson | 6 | 8.41e-09 ✗ | 9.49e-16 ✓ | 1.11x |

So for the quadrature methods: **sign rule is fast but wrong at order >= 4; resampling is exact
but the saving collapses to ~1.02x.** The mirror is worth having there only at order <= 2.

**`gl` is exact at every order** because it is a closed-form expression in the node values, so the
mirror is derived exactly (swap the Gauss-Legendre nodes) rather than via a series identity:

* order 2 — one node, the midpoint, is its own mirror: `Ω_mirror = Ω`;
* order 4 — two nodes swap: `Ω = S + K`, `Ω_mirror = S - K` where
  `S = (h/2)(A₁+A₂)`, `K = (√3/12)h²[A₂,A₁]`;
* order 6 — three nodes, `A₁ ↔ A₃`: recompute `_magnus_gl(An[:, ::-1], h, order)`.

Prototype code for all three is in §2.7 below.

### 2.4 A dead end, measured so nobody repeats it

**Making the cumulative integral direction-symmetric does not fix the order >= 4 error, and costs
1.6-1.7x.** Patching `_cumulative_integral` to average the left-to-right sweep with the reversed
right-to-left one left the error at **exactly 3.72e-08** — identical — while making trapezoid
order 4 run 1.71x slower. The integration *direction* was never the cause; §2.3 is.

### 2.5 Measured speed-ups

**Where the time goes** (64 slabs). The exponential is free; everything else is halvable:

| method | order | total | `expm` |
|---|---|---|---|
| gl | 2 | 0.80 ms | 0.04 ms |
| trapezoid | 4 | 9.46 ms | 0.01 ms |
| trapezoid | 6 | 31.64 ms | 0.01 ms |

**`gl`, cheap analytic Hamiltonian** (exact to 1e-15 throughout):

| slabs | order 2 | order 4 | order 6 |
|---|---|---|---|
| 32 | 2.71x | 1.85x | 1.35x |
| 128 | 1.38x | 1.39x | 1.19x |
| 512 | 1.09x | 1.23x | 1.06x |

*Treat the 32-slab row as noise* — sub-millisecond timings, and one row had order 4 faster than
order 2, which is incoherent.

**`gl` against `H_func` cost** (128 slabs) — the axis that decides the whole question:

| spectral modes / evaluation | order 2 | order 4 |
|---|---|---|
| 0 (cheap analytic) | 1.19x | 1.28x |
| 100 | 1.28x | 1.45x |
| 400 (expensive) | **1.61x** | **2.00x** |

With `f` = share of slab time in `H_func`, **speed-up = `1/(1 - f/2)`** — 1.11x at `f = 0.2`,
1.43x at 0.6, 2.00x at 1.0. This matched the sweep and is the formula to use for any given
Hamiltonian.

**Quadrature methods, sign rule** (fast but only valid at order <= 2): trapezoid 2.20-2.51x,
simpson 1.81-2.77x.

### 2.6 The Earth chord — palindromic to rounding, trivially fixable

On a uniform grid over a `costhz = -0.9` chord:

```
slab width               = 9.0806e+11
max|w - w[::-1]|         = 3.9062e-03   ->  RELATIVE 4.30e-15
max|mid - (L-mid)[::-1]| = 7.8125e-03   ->  relative 1.34e-16
```

**That is rounding, not physics** — exactly what NuOscProbExact describes ("about 1e-12 km on a
100 km slab"). *(An intermediate conclusion this session called it "a real asymmetry, not
rounding". Wrong.)* `w = (w + w[::-1])/2` makes it exactly palindromic — verified.

**Why it matters:** the gate must test **exact** equality, so today the Earth path would decline
and gain nothing. Symmetrising the grid is what unlocks it. **It also moves bit-identity** on
every Earth workload (by ~1e-15 relative), which needs justifying per workload.

Note the sampled potential is still not bitwise palindromic on either grid, because
`cumsum` positions do not mirror exactly — this did **not** matter for `gl` (still 1e-15) but is
worth knowing.

### 2.7 Working prototype (reproduces the shipped path to 1e-15 on `gl`)

```python
def gl_mirror(A, edges, order):
    """Evaluate A on the FIRST HALF only; the mirror slab's Omega follows by symmetry."""
    n = edges.shape[0]; m = n//2
    nodes = mg.gl_nodes(order)
    a, b = edges[:m, 0], edges[:m, 1]; h = b - a
    ts = a[:, None] + np.outer(h, nodes)
    An = A(ts.ravel()).reshape(m, len(nodes), 3, 3)     # half the H evaluations
    hh = h[:, None, None]
    if order <= 2:                       # midpoint is its own mirror
        Om_f = hh*An[:, 0]; Om_b = Om_f
    elif order <= 4:                     # nodes swap: symmetric part kept, commutator flips
        A1, A2 = An[:, 0], An[:, 1]
        S = 0.5*hh*(A1 + A2)
        K = (np.sqrt(3.)/12.)*(hh**2)*mg.commutator(A2, A1)
        Om_f, Om_b = S + K, S - K
    else:                                # order 6: reversal swaps A1<->A3, A2 fixed
        Om_f = mg._magnus_gl(An, h, order)
        Om_b = mg._magnus_gl(An[:, ::-1], h, order)
    Om = np.empty((n, 3, 3), dtype=complex)
    Om[:m] = Om_f; Om[n-m:] = Om_b[::-1]
    return mg._expm_stack(Om)
```

For the quadrature methods, replace the middle with `mg._magnus_terms_quadrature(Bt, order,
method)` and build the mirror by **reversing `Bt` along its sample axis and recomputing** (exact),
*not* by flipping term signs (wrong at order >= 4). Note `_magnus_terms_quadrature` returns shape
**`(order, m, d, d)`** — term axis **leading**.

### 2.8 Decision criteria already agreed with the user

* **Gate on exact symmetry of the slab inputs actually used.** A non-palindromic grid — including
  today's Earth chord — takes the existing path unchanged, so the feature ships with **zero
  bit-identity movement**. `bitident.py` must show 0 of 11 moved; anything else means the gate is
  wrong.
* **Symmetrising the Earth grid is a separate, later commit.** It is the only piece that moves
  bit-identity on the most-used path. Decoupling keeps risk and benefit apart.
* **Build when `f >= 0.46`** (speed-up >= 1.3x). Below that the gain does not pay for a second path
  through a dispatch layer that produced three silent-wrongness defects in one session.
* **Exact equality, never a tolerance** — `np.array_equal`. Copy NuOscProbExact's discipline
  verbatim: *"a tolerance here would silently return a different answer for a nearly-symmetric
  profile, which is the one thing an optimisation must never do."*
* **No slab-count floor is obviously needed** (unlike NuOscProbExact, whose composer adds a matrix
  multiply per slab and loses below ~15 slabs). Ours does strictly less work and was never slower
  — worst 1.06x. Untested below 32 slabs; if a floor is needed it ships as a measured constant.

### 2.9 The user's situation, and the open question

**The user mostly uses `gl`, the default.** That is the path where the mirror is exact at every
order, and where the saving is 1.1-1.4x on a cheap analytic Hamiltonian rising to 2.0x on an
expensive one. The user said they **may** end up with expensive Hamiltonians but has none in mind,
and asked for criteria to be set rather than supplied — hence §2.8.

**Immediately open:** the user asked whether it would help to run **Fable 5 with Ultracode**
pointedly on the §2.1-2.4 conclusions, to be "triple-sure", and asked to be answered *before* it
is activated. **That question was never answered** — the session moved to this handover. It is a
reasonable request: the conclusions are subtle (a sign rule valid only for k <= 2, a transpose
identity valid only for symmetric H, and one dead end) and every one of them is now backed by a
reproducible measurement in this document, which makes them cheap to re-verify independently.

---

## 3. Suggested next steps, in order

1. **Answer the Fable-5/Ultracode question** (§2.9) before anything else — the user asked to be
   consulted first.
2. **Commit the plan** (`PLAN_PALINDROMIC_PROFILES.md`, currently uncommitted on
   `dev-palindrome`). It carries the §2.1 and §2.6 findings, which are worth keeping regardless of
   whether the feature is built.
3. **Build the gated mechanism for `gl`**: `magnus.palindromic(*arrays)`, a slab-input test, and
   the `gl_mirror` path of §2.7. Gate on exact equality; assert `bitident` shows 0 of 11.
4. **Decline the mirror on quadrature methods at order >= 4** — §2.3 shows the exact version saves
   1.02x, which is not worth a code path.
5. **Only then** consider symmetrising the Earth grid, as its own commit with its own
   bit-identity justification.

---

## 4. Traps from this session, all paid for once already

* **A `SIGSTOP` corrupts a `perf_counter` stopwatch.** Pausing a timing battery for 75 min made
  the in-flight call report the whole pause as its execution time. Kill and re-run timing
  measurements; accuracy measurements resume fine.
* **Read minima, not medians, for timing on a shared machine.** Interference only ever adds time.
* **Sub-millisecond timings are noise.** If order 4 comes out faster than order 2, stop and use a
  bigger workload.
* **`pgrep -f foo` matches the shell issuing it.** Kill by PID from a targeted `ps`. Also check for
  `sleep` children reparented to `systemd` after their shell dies.
* **Never leave two waiters polling the same condition.** One `until … sleep … done` per condition.
* **Ask "does *any* row satisfy this", not "does the *worst* row satisfy it".**
  `bs05_energy_band.py`'s first verdict tested only the largest miss and concluded the opposite of
  the truth.
* **A monotone step is not a "hidden feature."** `find_hidden_features` computes
  `per_interval - endpoints`; the two nodes bracketing a step see its full height, so the statistic
  is zero at any width.
* **Verify a new diagnostic measures what you think.** The first aliasing probe took the *smallest*
  adjacent eigenvalue gap (the slowest oscillation) instead of the largest spread (the fastest,
  which is what aliases) — and agreed with a 4096-point reference to 1.000x, which looked like
  validation and was only evidence that the wrong quantity is smooth.
* **GitHub Actions was billing-blocked** for most of the session; every job "failed" in 2-9 s
  without starting. Making the repo public restored it. Do not diagnose a CI failure without
  reading why the job did not start.
