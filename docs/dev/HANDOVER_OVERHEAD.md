# Handover: reducing Magνs's per-call overhead

**Written:** 2026-08-09; §§9–10 added 2026-08-10 at the close of the backend work, §11 later the
same day at the close of the second review.
**If you are here to run a code review, read §10 and §11 first** — between them they list
everything two max-effort passes have already established, refuted or deliberately deferred, and
a third pass that rediscovers any of it has spent its budget on nothing. §11's last paragraph
names the two claims neither pass verified, which is where to start instead. Otherwise read
§§0–2 before touching anything — but note §§2–3 are superseded, see the DONE block below.

**The task for the next session.** Implement the Cayley–Hamilton matrix-exponential backend as a
selectable model with a numba implementation as the default, validated against degeneracies and
resonances, for 2ν through 5ν. Everything else in this brief is either already done or context
for that.

> **DONE, 2026-08-09 (same day, later session).** `src/magnus/expmkernels.py`, switch
> `magnus.magnus.EXPM_BACKEND` (`'auto'`/`'numba'`/`'eigh'`), 91 tests in
> `tests/test_expm_backend.py`. **6.8× on the exponential at N=108, 2.11× end to end** on a
> 60-energy PREM scan. §§2–3 below are superseded on four points, each of which cost something
> to find:
>
> 1. **No confluent (Hermite) form is needed.** §2 and §3 both assume coincident eigenvalues
>    require one. They do not: a Hermitian matrix is never defective, so a polynomial matching
>    exp on the *distinct* eigenvalues is already exact. With eigenvalues **sorted** and the
>    spectrum shifted to put the median at zero, the ill-conditioned coefficient multiplies a
>    matrix whose norm shrinks with the same gap, so its error is bounded by ε·gap and vanishes
>    as the gap closes. **No tolerance, no crossover, no near-degenerate branch** — 1e-16 at
>    splittings of 1e-2, 1e-6, 1e-10, 1e-14 and exactly 0 alike. **That holds at ‖K‖ ~ 1 only**,
>    which this note originally failed to say: the eigenvalue error scales with the norm, so a
>    clustered spectrum at large norm reaches 2.7e-07 against `eigh`'s 3.0e-11. See §9 and
>    `expmkernels.SEV_TOL`.
> 2. **`np.linalg.eigh` reads the LOWER triangle** (`UPLO='L'`), and `_expm_stack` admits input
>    anti-Hermitian only to 1e-12. The first version of the kernel read the upper triangle and
>    the two backends diverged by ~2e-12 — big enough to matter, small enough to read as
>    rounding. The kernel now reads the lower one.
> 3. **"19.1× at N=1" does not reproduce** through `_expm_stack`: it is 1.9×. `eigh` on one 3×3
>    costs 3.6 µs, and reaching it through `_expm_stack` costs 14.2 µs — the ~10 µs of framing
>    (the anti-Hermiticity test and its temporaries) does not shrink with the stack and is now
>    the single-point cost. **That is the next bottleneck, not the exponential.**
> 4. **The closed-form eigenvalues degrade to ~4e-9 at an exact degeneracy** and the exponential
>    stays at 2.5e-16 anyway, because interpolation error is *second* order in the displacement
>    of a coalescing node. Both halves are asserted; do not "fix" the eigenvalues.
>
> §6's memory corrections are applied. Still open: the ladder (§5), and the notebook-25
> re-measurement.

Your memory of this project is in
`~/.claude/projects/-home-mbustamante-Research-magnus/memory/` — start with `MEMORY.md`, then
`magnus-project-state.md` and `magnus-per-call-overhead.md`. The latter is now partly
superseded by this document; §6 says how.

---

## 0. Verify the base before starting

```bash
git -C ~/Research/magnus branch --show-current      # dev-overhead
git -C ~/Research/magnus log --oneline -3
git -C ~/Research/magnus status --porcelain         # empty
python -m pytest tests/ -q -n auto                  # 908 passed, ~10 min
```

`dev-overhead` is cut from `main` at `6e3251c` (PR #36, merged). **It has not been pushed** —
decide whether to push before or after the numba work.

**Measure on an idle machine.** This session produced three contaminated measurements before
the pattern was recognised: a parallel notebook build or a background `pytest -n auto` inflates
everything by 3–5×, and the *ratios between codes* shift too, so nothing is safe. Check
`uptime` first. All numbers in this brief were taken with load below ~1.5.

---

## 1. What is already done

Committed on `dev-overhead`:

* **`a7cdd07` — constant Hamiltonians, 73 µs → 19 µs.** The pipeline is skipped rather than
  tuned: for constant *A* every Magnus term past Ω₁ is a commutator of *A* with itself and
  vanishes, so Ω = −iH·Δ is the whole expansion. Also: the probe was rediscovering its own
  precondition (`if not callable(H_func)` *is* the proof of constancy); `eigh` reads one
  triangle so the explicit symmetrisation did nothing; `np.linspace` is 15× the cost of writing
  two endpoints down.

In the working tree, tested, awaiting the suite (see §5):

* **`magnus.ordered_product`** — the time-ordered slab chain by pairwise tree reduction instead
  of `reduce(np.matmul, U[::-1])`. 3.2× at N=108, 4.8× at N=2048, relative agreement 1e-15.
  Associativity licenses it; commutativity is neither needed nor assumed.
* **Per-slab constancy in `_magnus_gl`** — where a slab's GL nodes are bit-identical the
  commutators are exactly zero and Ω = h·A. Order 4: 35.7 → 5.4 µs; order 6: 102.9 → 8.7 µs.
  Results bit-identical when it fires. **Detection is exact equality, deliberately** — a nearly
  constant *A* has commutators carrying real information, and a tolerance would silently lower
  the order.
* **`magnus.cached_eval_mode` / `_eval_mode_for`** — the evaluation mode is a property of
  `H_func`, not of the call, so it is remembered weakly against the function. Verified 5 → 2
  `H_func` evaluations per call. Worth ~7% on a PREM lookup and about a third on notebook 19's
  Yukawa quadrature.

Five new tests in `tests/test_magnus_expansion.py` cover the first two.

---

## 2. THE TASK: a Cayley–Hamilton backend, numba by default

### Why Cayley–Hamilton and not eigenvectors

The session first proposed an analytic *eigenvector* route to preserve `_expm_stack`'s
advertised exact unitarity. **That objection was wrong and the measurement is in §6.** The
`eigh` path is *not* exactly unitary — 4.4e-16 at N=1, 3.8e-15 at N=4096, never zero, degrading
with stack size exactly as Cayley–Hamilton does. `_expm_stack`'s docstring says "yields an
exactly unitary result (probabilities that sum to 1 by construction)"; **that claim is false and
should be corrected as part of this work.**

So take the Cayley–Hamilton route: no eigenvectors, no null-space extraction, no eigenvector
multiplicity branching. exp(−iK) = a₀I + a₁K + a₂K² with the coefficients from interpolating
exp on the eigenvalues. Simpler, fewer branches, fewer places for a silent bug — which matters
more here than 1e-16 against 1e-15. It is also what NuOscProbExact does, so it is proven in the
sibling project (`~/Research/NuOscProb/NuOscProbExact/src/oscprob3nu.py`, and its
`fastkernels.py` for the numba side).

### Why it is worth doing

`np.linalg.eigh` costs **~1.25 µs per 3×3 regardless of stack size** — 1, 108 or 4096, the
per-matrix cost is flat, because it loops over LAPACK internally rather than vectorising. It is
**24% of a single 108-slab pass**, and `_expm_stack` as a whole is 33%.

Three routes were measured:

| approach | N=1 | N=108 | N=4096 | accuracy |
|---|---|---|---|---|
| `eigh` (current) | — | — | — | 4e-16 … 4e-15 |
| pure-numpy analytic | **0.2×** | 1.5× | 2.3× | 3e-15 |
| **numba kernel** | **19.1×** | **6.3×** | **6.2×** | 1e-15, unitary to 1e-15 |

The pure-numpy version *loses* below N≈100: ~20 numpy calls each pay dispatch overhead, and the
arithmetic is trivial. Only a compiled kernel removes that. **numba 0.60 is installed but
unused; `import numba as nb` sits commented out at `oscprob.py:298`** and it is not in
`pyproject.toml` — someone tried this before and stopped.

### The shape of the work

1. **Backend selection.** An `expm_backend` parameter reaching `_expm_stack`, plus a module
   default: `'numba'` when available, `'eigh'` otherwise, and `'auto'`. `HAVE_NUMBA` detection
   with graceful fallback so an install without numba still works. numba as an **optional**
   dependency in `pyproject.toml`.
2. **The kernel.** dim 2 and dim 3 analytically; **dim ≥ 4 delegates to `eigh`** — there is no
   practical closed form for 4×4/5×5 Hermitian eigenproblems. 4ν and 5ν therefore keep working
   but do not get faster, and the docs must say so rather than imply otherwise.
3. **Degeneracy, which is the crux.** Use `cache=True` on the jitted kernel; first-call
   compilation is ~0.7 s.
4. **Tests first** — see §3.
5. **Docs** — the `expm_backend` parameter in docstrings, a section in
   `implementation_details.rst` with what each backend does, the dim ≥ 4 limitation, the
   measured numbers, and the corrected unitarity claim.

---

## 3. Degeneracy is the whole risk, and there is a failing case waiting

A prototype kernel (§7) reached **6× and was completely wrong** — errors of 0.5, unitarity off
by 2.1 — because of a single conjugate in the determinant: `conj(b12)*b02` where it should have
been `b12*conj(b02)`. It still produced plausible, nearly-unitary matrices. That is the
silent-wrong-answer shape this package has been burned by repeatedly; see
`magnus-hybrid-certifies-while-wrong` and `magnus-crosscheck-cannot-reach-shared-blindness` in
memory.

With that fixed, the eigenvector prototype then **crashed with `ZeroDivisionError` on the very
first degeneracy case (H ∝ I)**: coincident eigenvalues make the null space 2-dimensional and
every row cross-product vanishes. Cayley–Hamilton avoids that failure mode entirely — but it has
its own, since Lagrange interpolation divides by eigenvalue differences. **Coincident
eigenvalues need the confluent (Hermite) form**, and near-coincident ones need care about which
branch is taken and where the crossover sits.

**Write these tests before the kernel.** They are what the prototype failed:

* H ∝ I (three-fold degenerate)
* two-fold degenerate, e.g. `diag(1, 1, 2)`
* the zero matrix
* a near-degeneracy sweep with splittings 1e-2, 1e-6, 1e-10, 1e-14, 0 — **this is the MSW
  resonance case**, where two eigenvalues genuinely approach each other, and it is the reason
  the user asked for resonance robustness
* random Hermitian stacks across dims 2, 3, 4, 5 and sizes 1 … 4096
* backend equivalence: `'numba'` against `'eigh'` to ~1e-14 on all of the above
* unitarity and row-sum invariants
* physics level: `osc_prob` agreeing between backends on a solar profile, an Earth chord, and a
  genuine NSI-induced resonance (notebook 12 has one)

---

## 4. What was tried and did NOT work — do not redo these

* **"Cancellation in ‖∫A dt‖ makes `suggest_n_slabs` underestimate."** False: ‖∫A dt‖ = 9.74 and
  ∫‖A‖dt = 9.76 on the PREM chord. No cancellation.
* **A tolerance-aware seed, n ~ θ^((p+1)/p)·tol^(−1/p).** Matched the first configuration
  (98 predicted, 108 actual) and then failed on four of five, over-predicting by 10–20×.
  Required slabs depend on profile structure, not total phase. Do not resurrect this without a
  much better model.
* **Per-slab constancy on PREM.** Implemented and correct, but **0 of 108 slabs** have identical
  GL samples — PREM is piecewise-*smooth*, not piecewise-constant. It fires on castle-wall
  profiles and `t_breakpoints`-delimited uniform regions, not on the Earth. Do not expect it to
  move the PREM number.

---

## 5. State of the ladder, and why it is lower priority than it looks

On the PREM chord the adaptive call climbs `[2, 3, 4, 6, 9, 14, 21, 32, 48, 72, 108]` — 11
levels, 319 slab-evaluations to deliver an answer at 108. The ladder is **77% of the adaptive
call**.

But it is largely amortised where it matters: a 20-point loop runs at 0.63× of 20 independent
calls, and batched `osc_prob_earth` uses **3.2 magnus calls per point** rather than 11, because
`osc_prob_energy_baseline` sets `A_eval_mode` once and warm-starts neighbours. Scans are what
the code mostly does. The 3× waste is a single-point cost.

Still open if you want it: `growth_factor_n_slabs` 1.5 → 2.0 gives 7 levels instead of 11 (254
against 319 evaluations, ~20%), but `MIN_EFFECTIVE_REFINEMENT = 1.25` is tied to it by a
documented measurement at `oscprob.py:370` and following. Read that before touching either.

---

## 6. Corrections to earlier records

`magnus-per-call-overhead.md` in memory records "344 µs against ~40 µs of irreducible algebra".
**Those absolute numbers were taken under load and are inflated ~4–5×.** On an idle machine the
same measurements are:

```
numpy eigh of the 3x3          3.6 us
scipy expm of the 3x3          6.8 us
NuOscProbExact                10.2 us
Magnus osc_prob, before       72.9 us     (6.5x NuOscProbExact)
Magnus osc_prob, after a7cdd07 ~19 us     (~1.9x)
```

The qualitative conclusion in that memory — generality overhead rather than the exponential —
holds. Update the file with these figures.

Notebook 25's stored speed/accuracy figures were also measured under load. They are not wrong
about the *shape* of the comparison, but the absolute microseconds are high. Worth re-running
`--only 25_` on an idle machine at some point; not urgent.

---

## 7. Where the prototypes are

The working numba prototype — corrected determinant, 6.3× at N=108, 1e-15 accuracy — is in the
session scratchpad at
`/tmp/claude-1000/-home-mbustamante-Research-magnus/c89d33ca-380f-45f9-bff7-14a9c4475a7e/scratchpad/nbk.py`.
**Scratch space is not durable**; if it is gone, the algorithm is the trigonometric solution of
the characteristic cubic for the Hermitian part, and §3 lists everything it failed.

Note that prototype takes the *eigenvector* route, which §2 has now rejected. Use it for the
eigenvalue computation and the timing baseline, not for the structure.

---

## 8. Traps from this session, all paid for once

* **Measure on an idle machine.** Three separate contaminated measurements, one of which
  (344 µs) reached a memory file and a notebook.
* **`git checkout -- notebooks/` reverts the generator**, because `make_notebooks.py` lives in
  the same directory as its output. It silently undid a fix that the next commit then claimed to
  make. Restore artefacts with `git checkout -- 'notebooks/*.ipynb'`.
* **Stage before running the suite.** `test_tree_matches_git` compares `TREE` against
  `git ls-files`, so an untracked file is invisible and 903 tests pass locally while CI fails.
  `python tests/test_file_tree.py` does *not* run that assertion — only pytest does.
* **This machine's numpy is 1.26; CI installs 2.x**, where `np.trapz` is gone. Test
  version-sensitive code by patching numpy to look like 2.0.
* **A consistency check between two things cannot tell you when both moved together.** The
  `sources_match` guard passed while the generator had been silently reverted — the artefacts
  agreed with a generator that was itself wrong.
* **`make_notebooks.py --only 19_,24_`** rebuilds just those: 38 min → 9.9 s. Use it.
* **`pgrep -f` matches the shell issuing it.** A kill that appears to fail may have worked.
* **The user asks for brevity.** Lead with the answer; keep the measurements exhaustive and the
  prose short.

---

## 9. The four inherited findings — closed

A max-effort code review of this branch (2026-08-09) produced fifteen findings. Eleven were
fixed in the commits that introduced them. These four were **inherited from `a7cdd07`/`2478fd7`
rather than from the backend work**, were held back so the backend commits stayed reviewable,
and are now closed by one refactor of the evaluation-mode cache. All four were confirmed by
execution before being fixed, and the reproduction was written before the fix.

They turned out to be one defect wearing four faces: a public cache with a wrong key, and a
private second copy of it that was doing the actual work.

1. **The key omitted the interval — the largest silent error the review found, 5.8e-02.**
   `probe_eval_mode`'s `'constant'` verdict means "sampling A across [t0, t1] gave the same
   matrix every time", which a wider interval can falsify. A two-layer profile that
   short-circuits when all requested positions fall in one layer probes as `'constant'` on a
   short baseline and `'vector'` on a long one, *for the same function object*. Keyed on the
   function alone, a mode learned on the short interval was served for the long one, and
   `_evaluate_A`'s `'constant'` branch then broadcast one sample over a profile that varies —
   with no spot-check, unlike the `'vector'`/`'scalar'` hints, which self-heal. Unitary,
   unwarned, and wrong only for one order of the caller's loop. Now keyed on
   `(function, t0, t1)`: measured loop-order dependence 5.8e-02 → **0.00e+00**.

2. **The lookup was unguarded, so callable Hamiltonian *objects* raised.** `WeakKeyDictionary`
   raises `TypeError` for a key that cannot be weakly referenced (`__slots__`) or hashed (a
   `@dataclass`, which sets `__hash__ = None`, or anything defining `__eq__`). All are ordinary
   ways to write a Hamiltonian and all worked before the cache existed. The guarded public
   version wrapped its lookup; the private copy in `oscprob` wrapped only its *store*. Now
   guarded, and tested for all three shapes: a cache must not decide whether a call succeeds.

3. **The comment claimed identity keying; `WeakKeyDictionary` keys by equality.** Two profile
   objects comparing equal shared an entry (P_ee 0.860706 against 0.858081). The container's
   behaviour is not a choice made here, so the comment is corrected rather than the code.

4. **The public `cached_eval_mode` was dead code and the private copy was its unguarded twin.**
   Its only call site sat in the `not callable(H_func)` arm of a ternary, and `H_func` is
   unconditionally rebound to a closure ~40 lines earlier, so the guard was never true
   (instrumented: 0 calls, while the copy got 1). `oscprob._eval_mode_for` is **deleted** — it
   also reached across a module boundary into `magnus._EVAL_MODE_CACHE` — and the one call site
   now uses `cached_eval_mode(..., key=H_func)`. `key` exists because the caller must wrap H in
   `lambda t: -1j*H(t)`, which would miss every time; scaling by a constant cannot change
   whether a function accepts an array, so the two share a verdict.

**What the interval key costs: nothing measurable.** The case the cache exists for is a
refinement ladder calling repeatedly at *one* interval, and that still probes once. A scan over
twelve distinct baselines probes twelve times, and came in faster than twelve identical calls,
because the probe is swamped by propagation.

Still recorded and not acted on: `_samples_identical`'s docstring claims `np.array_equal`
short-circuits on the first differing element. It does not — NumPy materialises the full
comparison — and the check costs 17.6%/7.4%/4.8% of an order-4 `_magnus_gl` call at
n = 1/108/2048, on every refinement iteration, on smooth profiles where it can never fire.
Separately, `integration_method='nonsense'` and `min_n_slabs=-5` are accepted by both the
constant engine and the per-point path; that is pre-existing and unchanged.

---

## 10. Handover for the code-review session (written 2026-08-10)

**Read this section first if you are here to run `/code-review max`.** The review has already
run once on this branch at max effort. Everything it found is fixed. Re-running it is
deliberate — the fixes were written and verified by the same agent that wrote the code — but it
must not rediscover the same fifteen findings, and §10.3 lists things that were tried and
refuted, which cost real time to establish.

### 10.0 State

Branch `dev-overhead`, cut from `main` at `6e3251c`, **never pushed**. `main` and `upstream` are
the same GitHub repo (identical tips); push to `origin`. Do not trust a tip hash written here —
every edit to this brief invalidates it. `git log --oneline main..HEAD` is the authority; these
eight commits carry the code, and the rest are edits to this brief.

| commit | what |
|---|---|
| `a7cdd07` | constant Hamiltonians, 73 µs → 19 µs — see §1 |
| `2478fd7` | `ordered_product`, per-slab constancy, the evaluation-mode cache — see §1 |
| `2debd51` | Cayley–Hamilton numba backend, `EXPM_BACKEND`, `SEV_TOL` gate |
| `fd0a5d5` | the `'constant'` engine, wrapper-overhead caches, `h_matt` fix |
| `e8b1d05` | notebook 25 (batched NuOscProbExact, PREM 3ν and 3+1) + full 26-notebook rebuild |
| `b68fad6` | order-6 conjunction test; corrections to this brief |
| `98024e7` | evaluation-mode cache: interval in the key, duplicate deleted |
| *(this one)* | the second review's four findings, plus the pre-existing ceiling bug they uncovered — see §11 |

The review's scope is all eight — the first two predate the backend work but are where §9's four
findings came from, so do not scope to the middle five alone.

Gates, all green as of §11's commit: **1056 tests** (`pytest tests/ -q -n auto`, 8 min 41 s),
`ruff check src/ tests/ notebooks/make_notebooks.py`, and `make clean html` with
`SPHINXOPTS="-n -W --keep-going"` — **clean**, not incremental, which is the only form of that
gate worth anything; see §10.5. At `98024e7` the same gates read **1044 tests**, and all 26
notebooks executed with **no accuracy column changed** — that notebook run still stands, since
§11's changes alter no computed number. `2debd51` and `fd0a5d5` were each verified to pass in
isolation in a throwaway worktree.

**Run only the gate whose inputs moved.** `docs/dev/` is invisible to all of them: it sits
outside `docs/source/`, so Sphinx never reads it, and `test_file_tree.py` collapses it to a
single entry, so neither its contents nor its filenames reach `TREE`. An edit confined to this
brief cannot change any gate's result, and re-running the suite to "confirm the tip" is 8.5
minutes that cannot come back different. `src/` or `tests/` moving is what obliges it — and
`src/` also invalidates every notebook cache, so CI re-executes all 26.

**26 of the 49 changed files are regenerated `.ipynb` JSON** (3 351 insertions, about half the
branch by line count) from `e8b1d05`, and three more are regenerated figures (`fig/*.pdf`,
`img/gallery/*.png`). All are build artefacts of `make_notebooks.py`. Review the generator; skim
the notebooks only for the accuracy columns. The reviewable surface is the other twenty files,
3 816 insertions:

```bash
git diff main HEAD -- . ':(exclude)*.ipynb' ':(exclude)fig/*' ':(exclude)img/*'
```

### 10.1 Results, so nothing is re-measured

Constant density, 3ν, 1300 km, 60 energies, interleaved with a control that returned 1.00×:

| | |
|---|---|
| exponential alone, N=108 d=3 | **6.7×** faster than `eigh` (13.3× at d=2, N=1024) |
| constant-density / vacuum scans | **6–25×** by flavour count (2ν 17.3/24.7, 3ν 15.5/18.9, 4ν 7.2/7.4, 5ν 6.0/6.2) |
| PREM 60-energy scan, end to end | **2.11×** |
| vs NuOscProbExact, batched scan | **1.10 µs/energy against its 1.44** — Magnus wins |
| vs NuOscProbExact, single point | 33.8 µs against its 19.9 — we lose, and it is wrapper parameter resolution, not arithmetic |
| PREM 3ν vs NuOscProbExact | it is ~20× cheaper per call; Magnus reaches 3e-10 where its O(h²) discretisation stalls near 6e-5 |
| PREM 3+1, eV-scale splitting | **~1000× slower, with a `MagnusConvergenceWarning`**, and inherent — see §10.6 |

4ν/5ν gain less because the kernel covers d = 2 and 3 only; there is no practical closed form
for a 4×4 Hermitian eigenproblem and there never will be.

### 10.2 What to re-review, and why

Ranked. All of it was written *and* verified by one agent, which is the reason to look.

1. **`98024e7` — the least-reviewed code on the branch.** Nobody but its author has read it. It
   changes a cache on the hot path for every callable Hamiltonian, **deletes** a function
   (`oscprob._eval_mode_for`), and reroutes its only call site. §9 has the four findings it
   closes.
2. **`SEV_TOL` and its gate** (`expmkernels.py`, and the fallback in `magnus._expm_stack`). It is
   **conservative by design**: it also declines large-norm spectra that are *not* clustered and
   would have been fine. A reviewer who "tightens" it to remove those false positives reopens a
   7440× accuracy hole. Check the reasoning in the constant's docstring, not the number.
3. **The decline conditions in `_osc_prob_scan_separable_dispatch`.** Under-declining is a silent
   wrong answer; over-declining is merely slow. One gap (a NaN baseline, missed because
   `nan < L0` is False) was found *after* the author believed them complete.

Lower value, already swept hard by the first review and found clean: the Cayley–Hamilton algebra
(three independent derivations of `det X`, the root ordering, the divided differences and `Z²`),
the UPLO question (both backends read the numpy-index **lower** triangle; the upper-triangle
answer is O(1) away, so that test has teeth), `ordered_product` (matches `reduce` to n = 2049),
and the `_VACUUM_H_CACHE` copy discipline.

### 10.3 Refuted by measurement — do not retry

Beyond §4's two items, which still stand:

* **Compensating `m` and `n` in the cubic does not help.** Making both exact to 60 digits left
  the error at 4.79e-08 where the plain path gave 2.87e-08; making the *eigenvalues* exact
  reached `eigh` parity (1.20e-11). The loss is in the `(m, n) → λ` map, not in the invariants.
* **Newton polish on the characteristic cubic cannot fix it either.** At a double root `p′ = 0`
  as well, so the floor is `√(ε‖K‖³)` — the floor we already have. And at the tightest cells the
  trigonometric starting point is *further from the root than the roots are from each other*
  (3.56e-04 against a 1.00e-04 separation), so two starts can converge to the same root and
  collapse the interpolation nodes.
* **The clustering cannot be gated on.** The danger is a **band**, not a tail: at *exact*
  degeneracy the kernel is fine (0.3–1.8× of `eigh`, because at `u = +1` the coincident pair
  comes out of `cos(±2π/3)`, bit-identical by symmetry), the damage sits at intermediate
  separations, and `1/(1−u²)` is largest exactly where there is no problem. Calibration found no
  separating value. Gating on the *scale* works and is what ships.
* **`SEV_TOL` is robust across bases.** Re-checked in 8 random unitary bases, not the one it was
  calibrated in: minimum unsafe `m` = 1.100e5 against the gate at 1e4 — 11× margin, zero false
  negatives anywhere.

### 10.4 How to measure in this repo

The machine cannot be quieted below load ~1.2 (the desktop app itself). So use the house method
from `implementation_details.rst`: **interleave the alternatives round-robin, report minima, and
carry a workload the change cannot touch as a control.** A control that returns 1.00× is the
evidence the ratios are readable. Absolute microseconds are not transferable; ratios are.

Never assert "better than `eigh`" from random spectra alone — that is precisely the evidence that
produced the false claim this session had to retract. `tests/test_expm_backend.py` now crosses
separation with scale; extend that grid rather than writing a new sweep.

### 10.5 Traps paid for, each once

* **`earth.distance_traveled_inside_earth` returns KILOMETRES**; every `osc_prob` baseline is in
  natural units (`×gd.CONV_KM_TO_INV_EV`). Passing the raw value does not raise — it returns a
  converged, unitary, meaningless answer for a chord a few metres long, and the refinement ladder
  then agrees with *itself* at every tolerance. **A self-convergence study that reads exactly
  `0.000e+00` is degenerate, not converged.** This contaminated four PREM measurements and a
  figure that reached the CHANGELOG.
* **A flat absolute tolerance is wrong when the achievable floor scales with ‖K‖.** `exp(-iK)`
  carries a phase of size ‖K‖, so the floor is `ε‖K‖`: at scale 1e5 that is 3e-11, and asserting
  1e-12 asserts better-than-possible. Two of this session's own checks failed on this.
* **A ratio to `eigh` is the wrong yardstick.** 19.6× at 8.7e-14 is harmless; 7440× at 6.7e-08 is
  not. The author's own acceptance criterion ("ratio ≤ 2×") would have *rejected* the working fix.
* **`git checkout -- notebooks/` reverts the generator**, which lives beside its output. Restore
  artefacts with `git checkout -- 'notebooks/*.ipynb'`.
* **Stage before running the suite.** `test_tree_matches_git` compares `TREE` against `git
  ls-files`, so an untracked file is invisible locally and red in CI. `python
  tests/test_file_tree.py` does *not* run that assertion — only pytest does. Regenerate both doc
  trees with `python tests/test_file_tree.py --write`.
* **A `:func:`/`:data:` role pointing at a *private* name fails the docs gate** when it sits in a
  docstring autoapi renders (CI runs `-n -W`). Either export the name or use double backticks.
* **An incremental Sphinx build cannot re-find a warning whose file it has cached**, so
  `make html` on an existing `build/` is not the gate — CI builds from a clean tree. This
  branch carried a broken `:func:` role from `2debd51` to the eve of its first push, green on
  every incremental check in between. Always `make clean` before believing the docs gate.
* **`make_notebooks.py --only 25_`** rebuilds one notebook in ~10 s against ~30 min for all 26.
  A change under `src/` invalidates every notebook's cache, so CI re-executes all of them; it
  fails on execution errors, not on output diffs.
* **numba is in the `test` extra**, so CI exercises the kernel instead of skipping 100+ tests.
  That lets the resolver hold numpy back a minor release; accepted deliberately.

### 10.6 Still open

* **The next per-call bottleneck is `_expm_stack`'s own framing, not the exponential.** `eigh` on
  one 3×3 costs 3.6 µs; reaching it through `_expm_stack` costs 14.2 µs. The ~10 µs is the
  anti-Hermiticity probe and its temporaries, which do not shrink with the stack — measured at
  43–69% of the call at N=1 and skippable for the two constant-H callers that build `Ω` from a
  Hamiltonian this package constructed. Deliberately deferred: it reworks the most-used function
  in the core, and it does not close the single-point gap to NuOscProbExact, which is wrapper
  parameter resolution and wants a thin entry point instead.
* **`_samples_identical`'s docstring claims `np.array_equal` short-circuits.** It does not; NumPy
  materialises the comparison. Costs 17.6/7.4/4.8% of an order-4 `_magnus_gl` at n = 1/108/2048,
  on every refinement iteration, on smooth profiles where it can never fire.
* **`integration_method='nonsense'` and `min_n_slabs=-5` are accepted** by both the constant
  engine and the per-point path. Pre-existing and unchanged — but note the *ceiling* half of
  this is now closed: see §11 for the three bounds the constant engine was skipping and the
  mislabelled condition underneath them. These two are what is left, and they are accepted
  **consistently**, which is why they are a wart rather than a defect.
* **PREM 3+1 is ~1000× slower than the closed form and warns.** The cost is flat across
  tolerances, which is the diagnosis: the ladder runs to its slab ceiling rather than converging.
  An eV-scale `Δm²₄₁` over an 11 000 km chord needs a slab width far below what the ladder
  reaches. Two causes, one fixable: 4ν falls back to `eigh` (a small constant), and a truncated
  Magnus expansion needs narrower slabs as the phase grows (which is what the expansion *is*).
  Notebook 25 says so in those terms. Do not treat it as an engineering defect.
* §5's ladder question is untouched, and `docs/dev/adversarial_batteries/RUN_P4.md` is still unrun.

---

## 11. The second review — what it found, and what it cleared

Run 2026-08-10 from a fresh session at max effort, on the eight code commits above. The point of
running it twice is §10.2's: the first review's findings were fixed *and* verified by the agent
that wrote the code. This pass was told about §9 and §10.3 up front so it would not re-derive
them, and it did not.

**Four findings, none of which changes a computed probability.** All five fixes below were
confirmed by execution before being written, and each carries a test verified to fail when its
defect is put back.

1. **A cache hit dropped the `DensityUnitWarning` that the first review believed it had fixed.**
   `vcc_func_from_rho_func`'s constant memo skips the conversion both unit guards live inside.
   The earlier repair mirrored only the `density_matter_is_in_g_per_cm3=True` arm — the safer
   one — and left the arm that catches an *undeclared* g/cm³ density, which returns exactly the
   vacuum probability with no tell in the numbers. **Re-emitting from the cache site cannot
   work**: `warnings.warn`'s `stacklevel` attributes the call to a different frame, the frame is
   part of the interpreter's registry key, and the imitation therefore printed a *second*
   warning under the default filter where an uncached call printed one. That was the earlier
   repair's own side effect, unnoticed because it was only ever tested under
   `simplefilter('always')`. Densities that trip either guard are now **not cached at all**, so
   both fire from where they always did. Verified against `main`: 3/3 under `'always'`, 1/3
   under the default filter, matching in both directions.
2. **Three refinement ceilings the constant engine accepted and `osc_prob` rejects** —
   `max_n_slabs=0`, `rtol`/`atol` ≤ 0, `max_num_loops=0`. Answers were never wrong; the *error
   contract* depended on whether the caller's density happened to be constant.
   `_refinement_params_rejected` mirrors the ladder and declines, so the message stays
   `osc_prob`'s own.
3. **`EXPM_BACKEND` did not survive a process boundary.** loky re-imports magnus in each worker
   at `'auto'`, and the `oscprob` wrappers expose no `expm_backend` parameter, so the global was
   the caller's only control and `n_jobs != 1` ignored it — worst for the one use it is
   documented for, since a backend comparison in parallel compared `'auto'` with itself. Carried
   by value into the worker. **Any new parallel entry point must do the same.**
4. **The evaluation-mode cache's per-interval dict was unbounded.** Weak keys bound the outer
   map; a module-scope Hamiltonian never dies, so 1000 distinct baselines retained 1000 entries
   (~184 KB) for the life of the process. Bounded at 256, cleared wholesale, matching its two
   siblings.

**And one pre-existing defect the mirror in (2) uncovered.** `osc_prob`'s ceiling check tested
`max_n_slabs` while its message named `max_n_tpts_per_slab`: that parameter was never validated,
`max_n_slabs` was bounded at `> 2` against a message promising `> 1`, and `max_n_slabs=2` was
refused in the name of a parameter the caller had not passed. Both messages encode the same rule
— each ceiling clears its own floor, `min_n_slabs` 1 and `min_n_tpts_per_slab` 2 — so the
condition was the wrong half and now names `max_n_tpts_per_slab`.

### 11.1 Cleared — do not re-open these either

* **The Cayley–Hamilton algebra.** `det X`'s cross term, the divided differences, the root
  ordering and all six entries of `Z²` re-derived independently against the code. Correct.
* **The constant engine against an independent oracle.** Engine-vs-engine is worthless here —
  both in-package routes reach `_expm_stack` and agree to *exactly* 0.0. Against
  `scipy.linalg.expm` it matches to 3.4e-15 across ν/ν̄, three densities, per-point baselines,
  L0 ≠ 0, a single scalar point and vacuum.
* **The d = 2 kernel has no `SEV_TOL` gate** — `_ch2_core` returns a hardcoded `sev = 0.0` — and
  does not need one: 1.9e-11 against `eigh`'s 6.6e-11 at ‖K‖ = 1e5, and 1.9e-09 against 3.3e-09
  at 1e7, two orders of magnitude past its documented range.
* **`float(VCC_func)` on an array-valued density** is unreachable: `validate_input_battery`
  rejects array `rho_func` first, and that code is older than this branch.
* **`verbose=None`** raises `TypeError`, and already did on `main`.
* **`ordered_product` on an empty stack** raises `IndexError` where `reduce` raised `TypeError`.
  Both raise; no caller can produce an empty chain.

**And one the gate caught rather than the review.** A `:func:` role in
`implementation_details.rst` pointed at `magnus.magnus._expm_stack`, a private name autoapi does
not document, which fails `-n -W`. It was introduced by `2debd51` and survived to the eve of the
first push because **every check of that gate had been an incremental build**, and Sphinx does
not re-read a file it believes unchanged — so once the warning had been emitted and its file
cached, no later `make html` could reproduce it. CI builds from a clean tree and would have
failed on the first push. See §10.5.

**Still not independently verified, by either review:** `SEV_TOL`'s d = 3 calibration — the
separation-by-scale grid that puts the gate at 1e4 — and the PREM 3+1 diagnosis in §10.6. Both
are taken on §10.3's word. A third pass wanting something to do should start there.
