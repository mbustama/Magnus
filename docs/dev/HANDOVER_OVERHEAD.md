# Handover: reducing Magνs's per-call overhead

**Written:** 2026-08-09, at the close of a long session. **Read §§0–2 before touching anything.**

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
