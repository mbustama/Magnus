# Pending paper edits

Running list of items raised but not yet applied to `main.tex`. Deliberately untracked,
like `audit-report.md` and `HANDOVER-audit.md` (a tracked file needs a `TREE` entry in
`tests/test_file_tree.py` and a regenerated `docs/source/installation.rst`).

## 1. `N` -> `N_{\rm slabs}` where appropriate

Sec. 4.5 introduces the slab count as bare `$N$`; Sec. 4.3 now writes `$N_{\rm slabs}^{-2}$`.
Pick one. `N_{\rm slabs}` echoes the code's `n_slabs` argument and disambiguates against the
other two things `N` currently denotes.

Slab count -- rename these:

| line | usage |
|------|-------|
| 424  | `N^{-2}` (order-two convergence) |
| 459  | `$N$` slabs (the definition, Sec. 4.5 opening) |
| 472  | `N^{-p}` |
| 484  | `N^{-4}`, `N^{-6}` |
| 486  | `N^{-p}`, `$N$`, `N^{-2}` |
| 488  | `$N$` x3, `\sqrt{N}` |
| 517  | `N^{-p}` |
| 1336 | composing `$N$` constant-density operators |

NOT the slab count -- leave, or give their own subscript:

| line | usage |
|------|-------|
| 528  | `N = 108`, `N = 4096` -- stack size for `eigh` timing |
| 712  | cumulative scan -- number of baselines |
| 1364 | "a further factor of `$N$`" -- nuSQuIDS phase samples |

## 2. Appendix rename, `n` -> `k` for the Magnus order

Sec. 4.1 now uses `k` (matching `\equ{magnus}`'s `\sum_{k\geq1}\Omega_k` and line 366).
The appendix at line 1506 onward still has `\Omega_n`, `n-1`, `n-2` in four places.
Mechanical; `n` is the flavour count elsewhere, which is why it moved.

## 3. Sec. 4.7, the 4x4 closed form

"there being no closed form for a $4\times4$ or $5\times5$ Hermitian eigenproblem"
contradicts the Abel--Ruffini passage in Sec. 1, which grants that the quartic is soluble.
Insert "practical", or split the two cases: the quintic has no solution in radicals at all,
the quartic has one not worth evaluating at these sizes. `expmkernels.py:37` has it right.

## 4. Term counts collide

Paper: `\Omega_7` has 17 terms, `\Omega_{10}` has 129 -- right-nested chains in the lower
`\Omega_m`, matching `magnus.py:1046`. But `expansionterms.py:41` says 26 at order 6, 211 at
order 8, 1918 at order 10 -- fully expanded commutator words in `A` alone. Both correct,
different objects. One clause in either place saying which is being counted.

## 5. "the cost panel"

Fig. 2's caption introduces panels positionally (Top left / Top right / ...) then uses
"the cost panel" once, unglossed, at the end. Body says "the lower left panel". Use the
positional name in both.

## 6. Unitarity numbers: maxima or means?

Sec. 4.7's "4e-16 for a single 3x3 and 4e-15 for a stack of 4096" reproduces as a *maximum*
over the stack -- median and mean are flat at 7.5e-16 / 8.9e-16 from N=1 to N=40960, so
nothing accumulates and the growth is extreme-value sampling. Reworded on that basis.
Still to confirm: whether the probability-output figures (3e-12 -> 1.6e-11 across four
decades in the number of points) are also maxima. If they are means, that clause needs
different wording, since a rising mean would be a real effect.

## 7. Two thresholds in the order-two paragraph

The slab count "stops growing below 0.2 GeV" and the resonance leaves the profile below
0.27 GeV. If these are the same feature, name it once.

## 8. Fig. 2, lower left: drop "Ref." from the legend

`notebooks/make_notebooks.py:13809`

    CODES = [('dop853', 'Ref.: Runge-Kutta order 8 (DOP853)', INK, '-', None, None)]
                        ^^^^^^ drop

to `'Runge-Kutta order 8 (DOP853)'`. Not cosmetic: DOP853 is a *competitor* in this panel,
one of the seven timed configurations, while the referee is the mpmath midpoint slab product
Richardson-extrapolated three times. Labelling it "Ref." is what led the Sec. 4.3 draft to
describe the 50-digit reference as "a stand-alone Runge-Kutta order-8 solver, DOP853 as
implemented in scipy" -- which cannot be right, scipy being double precision.

The comment three lines above says the same thing ("the solver they are all measured
against"); worth rewording to "the solver they are all timed against" so competitor and
referee stay distinct.

Requires regenerating Fig. 2. Cheap relative to Fig. 13 but not free -- batch it with any
other figure edits rather than doing it alone.

## 9. Sec. 4.4: the order-six display is cited to the wrong paper

The three-commutator form the code implements --

    C1 = [a1, a2];  C2 = -(1/60)[a1, 2a3 + C1]
    Omega^(6) = a1 + a3/12 + (1/240)[-20a1 - a3 + C1, a2 + C2]

-- is Eq. (251) of the Blanes-Casas-Oteo-Ros review (`Blanes:2008xlr`, arXiv:0810.5488),
with the alphas its Eq. (257). The review attributes it to BCR, BIT 42 (2002) 262, and
states on p. 96 that three commutators is the *minimum* for sixth order.

`Blanes:2000bit` (BIT 40, 2000) gives a *different* order-six scheme, its Eq. (3.10), in
terms of B^(0), B^(1), B^(2) and needing four commutators. Same nodes, different
arrangement. So the display should cite `Blanes:2008xlr`, not `Blanes:2000bit`.

Add to the .bib (used for the minimality claim, and the order-8 scheme):

    @article{Blanes:2002opt,
          author         = "Blanes, S. and Casas, F. and Ros, J.",
          title          = "{High order optimized geometric integrators for linear
                            differential equations}",
          journal        = "BIT Numer. Math.",
          volume         = "42",
          pages          = "262--284",
          year           = "2002"
    }

Also: "the Gauss-Legendre schemes run out at three nodes" is true of this implementation,
not of the method -- BCR 2000 Sec. 4.1(iii) constructs an 8th-order four-node scheme
explicitly (Eqs. 3.21, 3.22, 4.3; ten commutators, six in the 2002 paper). Word it as a
property of Magnus.

## 10. Future direction: order-eight Gauss-Legendre collocation

Concrete and bounded. Four nodes at v1 = (1/2)sqrt((3+2 sqrt(6/5))/7),
v2 = (1/2)sqrt((3-2 sqrt(6/5))/7); B^(0..3) from BCR 2000 Eq. (4.3); Q1..Q7 from its
Eq. (3.21); assemble by Eq. (3.22). A fourth branch in `_magnus_gl`, a `_GL4_NODES`,
and `MAGNUS_EXP_ORDER_MAX_GL = 8`.

Verification is ready-made: the existing Simpson order-8 path is an independent oracle at
the same requested order, the mpmath reference bounds both, and h^8 convergence is a sharp
test. Payoff is the cost driver Sec. 4.3 already measured -- four evaluations of H per slab
against Simpson's sixty-to-five-hundred, which is what makes order ten cost over a minute
per probability against about 12 ms.

## 11. Order-8 collocation: what is verified, and the R^(4) derivation

Independently verified (Fable, reading BIT 40 (2000) at 300 dpi plus its text layer, and
re-implementing from scratch):

* The order-8 scheme of BIT Eqs. (3.21), (3.22), (4.3) as transcribed in item 10 is correct
  line by line. Local-error slope ~9 across six seed/t0 combinations, and the same with exact
  integrals instead of the 4-node quadrature, so order 8 does not depend on the quadrature.
* Q7 = -(1/42)[B0,[B0, Q3 - (1/3)Q4 + h Q5]]. BIT PRINTS -1/42; there is no typo in the
  coefficient. (An earlier note here claimed otherwise -- that was a misread of the stacked
  fraction at low resolution.) The one real defect in (3.21) is a stray comma, "Q3 - (1/3)Q4,
  +hQ5", which must be a plus.
* (3.21)-(3.22) reproduces the raw R-form (3.11)-(3.20) to 1.4e-16.
* Review Eq. (251) and BIT Eq. (3.10) are both order 6 but differ by O(h^7) on the same three
  nodes -- different constructions, confirming item 9's citation fix.

Alpha notation for order 8 is a DERIVATION, not a transcription: R^(4) appears in neither
paper. Derived exactly (Gauss-Jordan in rationals on T^(4)_ij = (1-(-1)^(i+j))/((i+j)2^(i+j))):

    R^(4) = [[  9/4,    0,  -15,     0],
             [    0,   75,    0,  -420],
             [  -15,    0,  180,     0],
             [    0, -420,    0,  2800]]

so alpha1 = (9/4)A0 - 15 A2, alpha2 = 75 A1 - 420 A3, alpha3 = -15 A0 + 180 A2,
alpha4 = -420 A1 + 2800 A3, with A^(i) = h B^(i). Its even sub-block reproduces R^(3), as it
must. NOT YET VERIFIED: the alpha-form Omega^(8) built on this. Rewriting (3.21) through
R^(4) needs its own convergence check before anything is printed.

Two traps for any hand conversion between the two papers:
* review A^(i) = h * B^(i) (BIT 3.4) -- normalizations differ by one power of h;
* BIT measures node offsets from the slab MIDPOINT, the review from the slab START.

---

## DONE 2026-09-05: re-ran the Earth column of Fig. 11

`notebooks/prem_chord_common.py`'s `vcc` evaluated the matter potential with a Python
list comprehension over every position: 6.1 us a point, against 0.063 us for the
vectorized form now in place. Same numbers, to 3.3e-16 -- only the timings differed.

**The handicap fell on one code only.** Magnus is driven through that helper
(`gen_prem_benchmarks.py` passes `prof['vcc']`); NuOscProbExact's Route A
(`append_npe_rtol_prem.py`) builds its potential from `earth.earth_slabs` and never
calls it. Share of Magnus's stored Earth time that was the helper, at 2 flavours:
38% at rtol 1e-3, 44% at 1e-4, 49% at 1e-6, 61% at 1e-8.

**What needs re-running**, on an idle machine:
* Magnus orders 4, 6 and 8, all four flavour cases, Earth chord only
  (`gen_prem_benchmarks.py`, after deleting the three Magnus series per case).
* Then notebook 28, then the figure.

**What does NOT:** NuOscProbExact (never touched the helper); the whole exponential
column (its `vcc` was already vectorized); `prem_chord_reference.json` (the helper
returned identical values, so the reference is sound).

**Paper numbers that move with it**, all in the two paragraphs after "The ceiling is
the flavor count" in `\subsection{A smooth profile...}`:
* "5.9 us per Magnus slab at order 4" and the "0.045 us" it is set against;
* the factor of 130, and the "366 overtakes 130" sentence;
* the ~600 us fixed refinement cost;
* the break-even costs per sample: 6.5 us, 1.7 us, 144 ns.
Unaffected: the 366x slab ratio and the sample counts (972 vs 34,782; 2606 vs 278,494),
which are geometry rather than timing.

Provisional arithmetic, not a re-run: corrected, Magnus at 1e-6 is about 515 us against
NuOscProbExact's 603 us, so it is already ahead there and the crossover sits nearer 1e-6
than 1e-8. Re-measure before printing any of it.

### And the mirror defect, in the exponential column (found 2026-09-05, verified here)

`gen_profile_benchmarks.npe_points` builds the Hamiltonian OUTSIDE the timed region:
`v = prof['vcc'](mid)` and the `H` broadcast happen before `call` is defined, so
`timed(call)` measures `probabilities_Nnu_slabs` alone. `magnus_points` times the whole
`osc_prob_matter_std_potential` call, sampling and H construction included.

Measured cost of the untimed part, as a percentage of what IS timed:

| | n=2048 | n=32768 |
|---|---|---|
| 2nu | 227% | 379% |
| 3nu | 68%  | 170% |

So NuOscProbExact's exponential-column times are understated by roughly 1.7x to 4.8x.

**Both columns of Fig. 11 therefore favour the closed form, for unrelated reasons:**
the exponential column does not time NPE's H construction; the Earth column charged
Magnus for a slow helper NPE never called. The Earth NPE route (`append_npe_rtol_prem`)
IS fair -- its `evaluate` builds H inside the timed call.

Fix before re-running: move the `v = ...` and `H = ...` lines inside `call` in
`npe_points`, so both codes are timed from the profile function onward. Then the
exponential column needs re-running too, not just the Earth one.


---

## Outcome of the Earth-column re-run (2026-09-05)

Magnus 1.24x to 2.89x faster than the contaminated numbers, median 1.72x. Errors
unchanged to round-off (max absolute shift 1.08e-14 in a probability, consistent with
re-associating a product over ~1200 slabs). NuOscProbExact untouched.

**The crossover moved about three decades on the Earth chord**, from ~1e-8 to between
1e-4 and 1e-6. At 2 flavours, order 4: NPE faster by 5.3x at rtol 1e-3 and 2.8x at 1e-4;
Magnus faster by 1.5x at 1e-6 and 7.8x at 1e-8.

Paper numbers updated in `\subsection{A smooth profile...}`:
* per-slab cost 5.9 -> **1.5 us** (Magnus), 0.045 us (closed form) unchanged;
* cost ratio 130x -> **33x**;
* fixed refinement cost 600 -> **220 us**;
* break-even per sample: 144 ns at 1e-6 -> **0.5 us at 1e-4, and Magnus already ahead
  from 1e-6 down**;
* and the paragraph no longer claims one crossover for both profiles: it is near 1e-8 on
  the exponential, near 1e-5 on the Earth chord.

STILL OPEN: `gen_profile_benchmarks.npe_points` builds H outside its timed region. That
series is stored but no longer plotted, so it does not affect the figure -- but anyone
reusing that function inherits the bias, which is how it was found.

---

## POSSIBLE NEXT STEP after the commutator: thread the batched scan's chunk loop

Measured 2026-09-05, `dev-overhead`. Not implemented; recorded so the reasoning
is not lost.

**What `n_jobs` does today, which is not what it looks like.** At
`src/magnus/oscprob.py:4529` the batched dispatch returns `NotImplemented`
whenever `n_jobs != 1`. So on the common scan -- many energies, one shared
baseline -- setting `n_jobs > 1` does not add parallelism to the batched path,
it *disables* it and falls through to the per-point engine fanned across joblib
workers. Verified by instrumenting both the dispatch and `Parallel`:

| case | n_jobs=1 | n_jobs=10 |
|---|---|---|
| 12 energies, one baseline | batched taken, no joblib | batched DECLINED, 1 joblib fan-out |
| 12 energies, own baselines | declined, serial per-point | declined, 1 joblib fan-out |

It is a trade, not an addition: array batching (all energies sharing one set of
potential samples) is exchanged for process parallelism. The docstring at
`oscprob.py:2777` already says the default "is usually fastest" and that
`n_jobs` "is not a pure performance knob", but no measurement stands behind it.

**Why the two could be combined.** The batched path already chunks over
energies inside each refinement level (`for i0 in range(0, len(active), chunk)`,
sized by `_tile_for_working_set`). Those chunks are independent within a level
and only rejoin for the convergence test. So the place to parallelise is the
INNER chunk loop, not the outer point loop where `n_jobs` is wired.

Use **threads, not processes**: `expmkernels._jit` is
`nb.njit(cache=True, fastmath=False, nogil=True)`, so the compiled kernels
release the GIL. Threads share the arrays; processes would pay to pickle a large
complex stack each way.

**Semantics would be preserved** if the parallelism stays strictly inside a
level and the convergence check stays serial across all chunks: the slab count
is still set by the hardest energy in the whole stack, so refinement decisions,
warnings and answers are unchanged. That is exactly what today's `n_jobs` does
NOT preserve, since it changes which engine runs.

**Caveats.** Memory binds -- `_tile_for_working_set` sizes a chunk to a working
set and N threads need N of them, so thread count and chunk size trade off. The
gain is bounded by the fraction of time in nogil kernels versus GIL-held numpy.
And the comment above `_jit` chose `parallel=False` precisely to avoid
oversubscription inside joblib workers, so any threading here must not nest
inside a fan-out.

**Do it after the kernels, not before.** Each kernel that lands raises the nogil
fraction and gives threading more to work with -- the composition is one
already, the commutator may become another.

### What it would take, assessed 2026-09-05 -- OPTIONAL, not committed to

Read-only analysis of the scan, the tiling and memory guards, the numba module,
the warning machinery and the test suite. Nothing forbids threading the chunk
loop; four things would have to change first, and the engine's own documentation
argues the ceiling is bandwidth rather than cores.

**Would have to change:**

1. **Divide the tile budget by the worker count** -- pass
   `max_entries = BATCH_WORKING_ENTRIES // n_threads` to `_tile_for_working_set`.
   Provably cannot move a result: `test_tiling_the_working_set_changes_no_number`
   establishes that tiling changes no number.
2. **Make `magnus._working_set_chunk` concurrency-aware.** THE UNPLEASANT ONE. It
   reads `_available_memory_bytes()` per call and claims `available/2`, so N
   threads collectively claim N x half: the guard under-counts by exactly the
   thread count. Bites on trapezoid/simpson at order >= 8 (~142 MB per thread at
   order 8, ~530 MB at order 10) and at 5 flavours with large slab counts, where
   the tile degenerates to one energy. That is the OOM class
   `docs/dev/BUG_IP_EXP_MEMORY.md` exists to prevent. The default `'gl'` branch
   has no working-set guard at all. The fix needs a signature change on
   `evolution_operators_from_samples` or a new module global -- and a global is
   the shared mutable state the rest of the design avoids.
3. **Emit slab-norm warnings from the serial join**, collecting per chunk.
   Otherwise "shown once per session" becomes a check-then-act race and
   `stacklevel=4` resolves into the pool's frames rather than the scan's.
4. Fix the stale tile comment -- DONE, `f664ef2`.

**Needs no work:** correctness and bit-identity (chunks copy in, write disjoint
slices of `P_new` out, no cross-chunk reduction, so scheduling cannot move a
bit); `EXPM_BACKEND`, only read on this path, its mutation confined to loky
worker processes; `VCC_func`, sampled serially before the loop, so user-supplied
callables never run in a worker thread; numba, whose first-call compilation
serialises behind its own lock.

**The argument against.** The `BATCH_WORKING_ENTRIES` docstring records that
these kernels are memory-bound and that the 1 MB tile was tuned for cache
residency N threads would share -- so change (1), the safe one, works against
the thing being optimised. And the available parallelism moves per refinement
level: early levels fit all energies in one tile and offer none at all, while
later levels have fewer active energies to spread.

Sequencing: behind the commutator, and behind measuring the n_jobs-versus-batching
trade, which may show the question is smaller than it looks.

## DONE 2026-09-05: the batched Jacobi eigensolver for 4nu and 5nu

**Landed in commit 334ab33, v1.0.9.** The gate this section called for was run
before implementing: `err_k <= 10*err_e` held in all 364 cells, worst ratio 5.49x
on the prototype and 6.39x on the kernel as shipped. Measured end to end
afterwards at 1.81-1.95x (4nu) and 1.46-1.59x (5nu), with two and three flavours
unmoved at exactly 0.0 -- the negative control the whole change rests on. Worst
probability shift 1.557e-12, below the 2.8e-12 the survey predicted. Full suite
1285 passed.

The record below is kept as written, because the reasoning it corrects is the
part worth preserving.

Surveyed 2026-09-05 by Fable; **not started, and not to be started without the
author's say-so**, because unlike every other optimisation of this session it
changes numbers.

**The premise that has expired.** `expmkernels.supports_dim` stops at d = 3, and
its docstring gives the reason: "A 4x4 or 5x5 Hermitian eigenproblem has no
practical closed form, so 4nu and 5nu stay on eigh." The first clause is a fair
engineering judgement -- NuOscProbExact *does* solve the quartic in closed form,
by Euler's resolvent cubic from the SU(4) invariants, but its own
`psi_roots_4nu` docstring records that those roots "carry only the accuracy that
I_2, I_3, I_4 carry", which is why that code also carries a double-double
arithmetic suite and a Newton root-polishing pass. Exact at four flavours is not
cheap: NuOscProbExact costs 1.70 us/slab/E at 3+1 against 0.130 at 3nu.

What has expired is the *conclusion*. Magnus does not need a closed form to win
at d = 4, because the closed form is not what it is losing to: `np.linalg.eigh`
is 66% of the whole d = 4 call and 63% at d = 5, and the ~2.3 us per 4x4 is
LAPACK `zheevd`'s fixed overhead rather than arithmetic.

**The proposal.** A numba `_jacobi_expm_core(K, out, lam) -> sev` in
`expmkernels.py`, same contract as `_ch2_core`/`_ch3_core`, dispatched from
`expm_herm_stack` for d in (4, 5); cyclic complex-Hermitian Jacobi, warm-started
from the previous matrix's eigenvectors (consecutive matrices are consecutive
slabs of one energy, so they arrive nearly diagonal), with modified Gram-Schmidt
on the warm-start basis at every step. The MGS is load-bearing, not tidiness:
without it non-unitarity compounds along the chain, 2.3e-11 after 13k matrices
against 3.9e-14 with it.

**Claimed gain: 1.77x at 3+1, 1.38x at 3+2**; stacked with the fusions, 2.09x at
3+1, which moves Magnus from 1.53x *slower* than NuOscProbExact to about 0.73x
-- 1.247 against 1.70 us/slab/E. Those numbers are Fable's and have not been
reproduced independently.

**Why it is the author's call and not a routine optimisation.** It is not
bit-identical. Worst full-matrix probability shift 2.8e-12, and the finest-grid
floors move (expo d=4 7.40e-12 -> 7.71e-12, d=5 2.44e-11 -> 2.56e-11), while
reference error is unchanged wherever discretization dominates. The argument for
accepting it is that the shipped d <= 3 backend's own documented guarantee
against `expm` is 5e-12 absolute, so the shift falls inside a tolerance the
package already accepts for a backend swap rather than opening a new class; and
that Jacobi is backward stable with no conditioning cliff, so unlike the
closed-form kernels it needs no `SEV_TOL` gate (return sev = 0.0, with the
30-sweep cap escalating to the existing eigh fallback in a corner never
observed).

**What the change can and cannot touch, checked 2026-09-05.** The often-quoted
worry -- that this would undermine the paper's claim to reach 2.9e-13 -- does not
apply, and an earlier draft of this note had it wrong. `docs/source/comparison.rst`
labels that table "Exponential profile, 3nu"; the reach claim is a *three-flavour*
result, and it rests on the closed-form d = 3 kernel, which this change does not
touch. `supports_dim` would go from {2, 3} to {2, 3, 4, 5}: two and three flavours
keep exactly the backend they have now.

The claims that *are* exposed are the four- and five-flavour ones, and they sit
far above the shift: "PREM 3+1, self-convergence depth" quotes the referee's own
floor at 4e-7, and "PREM 3+1, cost" is a cost claim this change would improve
rather than threaten. The paper's accuracy oracle is a DOP853 integration at
rtol = 1e-12, atol = 1e-14, so a 2.8e-12 probability shift sits *at* the
resolution of the instrument the paper validates with -- it could not be
distinguished by the figure that would have to show it.

**The real exposure is `tests/test_expm_backend.py`, and it is not readable off
the page.** That file holds 32 tight-tolerance assertions, several at 1e-14
against `scipy.linalg.expm` and at 1e-13 against the eigh backend, currently
exercised only at d <= 3. Report 01 recommends extending them to d = 4, 5, which
is right -- but it also measures Jacobi's operator-level dU against eigh at
1.7e-11 at norm 1e4, which is above both those tolerances. Whether that is a
problem depends entirely on the norms the existing tests use, since the report
also states that eigh has the *same* eps*||K|| scaling and that the two agree
cell by cell. Somebody has to measure it at the tolerances actually asserted
rather than reason about it. **That measurement is the gate on this decision**,
and it has not been done.

The battery to run is specified in
`docs/dev/overhead_survey/BRIEF_jacobi_tolerance_gate.md`, including the control
that makes it meaningful: eigh must be scored on the identical battery, or a
failure of an assertion that simply does not transfer to higher dimensions will
be misread as a failure of Jacobi.

The remaining argument against is the general one: a backend that is occasionally
worse than eigh would cost more than the speed is worth, and the only way to know
is the battery above.

**Tests that would have to change:** `test_dim_four_and_five_delegate_to_eigh`
and `test_public_kernel_entry_point_refuses_unsupported_dimensions` pin the
current routing; the ValueError boundary moves from ">= 4" to ">= 6". The
agreement, degeneracy, ascending-eigenvalue and no-worse-than-eigh tests should
be *extended* to d = 4, 5 rather than edited.

Full write-up: `docs/dev/overhead_survey/report-01-jacobi-expm-4nu-5nu.md`;
final prototype kernel in `docs/dev/overhead_survey/prototypes/jacobi_proto5.py`,
accuracy battery in `jacobi_acc.py`, end-to-end A/B in `jacobi_e2e.py`. Start
from that directory's `README.md`, and from report 08 before quoting any stacked
speedup.

## PARKED 2026-09-05: the 3+1 cost numbers in the paper are stale after Jacobi

**Parked at the author's instruction -- do not act on this without being asked.**
Recorded here so it is not rediscovered, and so that whoever next regenerates a
figure knows to fold it in.

Measured 2026-09-05, after the Jacobi 4x4/5x5 backend landed in the working tree.
**Nothing here is wrong yet -- it becomes wrong the moment that change is
committed and the figures are not re-run.**

**What moves.** Magnus at four flavours is 1.81-1.95x faster end to end
(marginal us/slab, PREM and the smooth profile, eigh arm vs Jacobi arm
interleaved in one process, control drift +3.6%); at five flavours 1.46-1.59x.
Two and three flavours do not move at all -- `max|dP| = 0.000e+00` and 0.98-1.00x
timing, which is the negative control on the whole change.

**Numbers that therefore need re-running or re-stating:**

1. `docs/source/comparison.rst`, "PREM, 3+1": "**Cost: NuOscProbExact, by about
   400x.** 56 000 us per probability against 127." At 1.9x the 56 000 becomes
   about 29 000 and the ratio about 230x. Both the sentence and the bolded
   headline change.
2. `resources/paper/main.tex` \tabl{summary}, the row "PREM $3+1$, cost &
   {\tt NuOscProbExact} by $\sim\!440\times$". Same correction, roughly 230x.
3. **Figure 11's $3+1$ Earth panel.** The Magnus cost curve shifts down by the
   same factor; the crossover with NuOscProbExact moves with it. This one is a
   re-run, not an edit, and re-running Fig. 11 is expensive.
4. Any five-flavour cost statement, by 1.5x. There is no competitor at five
   flavours, so this is self-referential rather than comparative, but the
   absolute microsecond figures still move.

**What does NOT move, and this is worth stating because it is the natural
worry.** The reach claim -- "the slab product floors at 2.5e-11 and Magnus
continues to 2.9e-13" -- is a *three-flavour* result (`comparison.rst` labels
that table "Exponential profile, 3nu"), and three flavours keep the closed-form
kernel untouched. The 3+1 *accuracy* claims also stand: the residual against the
referee is 4.5e-08, four orders above the 1.557e-12 worst probability shift this
change introduces.

**Sequencing.** The paper's exponential-backend passage (`main.tex` around line
648) describes a Cayley-Hamilton kernel "for $2\times2$ and $3\times3$" and stops
there. That is not falsified by this change, but it is now an incomplete
description of the code: four and five go to Jacobi. One or two sentences, in the
same place.

Decide whether to land the Jacobi backend *before* the next figure re-run, so
that the two are paid for once rather than twice.

## PARKED 2026-09-05: does methodology.rst's order-selection argument still hold?

**Parked -- not investigated, and not to be without being asked.** Noticed during
the docs audit that followed the Jacobi backend.

`docs/source/methodology.rst`, in the order-selection study, argues:

> Evaluation count turned out to be a poor proxy: the fixed per-slab overhead
> (array setup, the eigendecomposition for the matrix exponential, the slab
> product) outweighs the node count, so fewer slabs matters more than fewer
> evaluations.

**All three of the costs it names got cheaper on 2026-09-05**, and by different
factors: array setup and the slab product through the composition, commutator and
fused-Omega kernels; the eigendecomposition through the Jacobi backend at four and
five flavours only. The compound is 8.6-12x at two and three flavours against
1.9-2.9x at four and five, so the balance the argument rests on has not merely
shrunk, it has shifted *unevenly across flavour count* -- and the study's own
conclusion is a crossover ("order 2 wins on evaluations at 1e-4 but loses on wall
time"), which is exactly the kind of claim a shifted balance can move.

**Why it was left alone rather than edited.** The conclusion probably still holds:
node count did not get cheaper either, so the ratio may be roughly preserved. But
"probably" is not what that paragraph asserts, and settling it means re-running the
order-selection study rather than reasoning about it. Asserting either way from the
armchair is the failure this project has a name for.

If it is ever re-run, note that the wall-time arm is what produced the table above
it in that file, so the table moves with the argument.

## ACHIEVING MATURITY -- the five things left, assessed 2026-09-05

Written after a session that shipped v1.0.6 through v1.0.10. The engineering
discipline around this code is mature; these are the substantive gaps.

### 1. The Notebooks gate reports green on a broken notebook -- DIAGNOSED

**Not flaky. Consistently broken, and cached over.** `25_magnus_against_other_codes`
raises `KeyError: 'Magnus, order 6'` at `DIAL_STYLE[series['name']]`. The benchmark
file `notebooks/external_profile_benchmarks.json` carries five series names --
`Magnus`, `Magnus, order 6`, `Magnus, order 8`, `NuOscProbExact`,
`NuOscProbExact, rtol` -- and the notebook's `DIAL_STYLE` defines three. Orders 6
and 8 entered the data when they were added to Fig. 11; the notebook never learned
about them.

**Why it looks intermittent.** The workflow skips any notebook whose source
fingerprint matches a marker in `.nbcache`, and that cache persists between runs
through `actions/cache`. A run is green whenever notebook 25's *source* happened
not to change since it last passed -- which was before the data grew. So the gate
alternates with the contents of the cache rather than with the health of the code,
and a green Notebooks run currently means "not re-run", not "works".

**The fix is two dictionary entries** in the notebook's `DIAL_STYLE`, plus a
decision about the cache: a marker keyed on notebook source alone cannot see a
change in the *data* the notebook reads. Folding the benchmark JSON's hash into
the fingerprint would close that hole.

`29_magnus_pseudo_dirac` also fails on the same shard, from a different cause,
not yet examined.

### 2. The cumulative branch has never carried default traffic

Reachable in released code via `cumulative=True`, and `oscprob.py` records that
making it the default surfaced two latent defects within minutes -- a missing
scalar squeeze returning `(1, d, d)`, and a `convergence_info` keyword forwarded
to an engine that rejects it. Reports 05 and 06 propose speeding this path up;
the audit should come first. A faster kernel on an unaudited branch is the wrong
order of work.

### 3. A documented silent-wrong mode remains

`strategy='auto'` has a narrow-feature case that returns a confident wrong answer,
recorded as unfixable on a fixed grid. Known and documented is far better than
unknown, but a code that can be quietly wrong is not finished.

### 4. The sterile projector takes one scalar

`H_matt` factorises as `VCC(l) x P`, so the projector cannot follow the Earth's
layered `Y_e` the way the density does. Shipped with a warning and a deferred
per-layer fix.

### 5. The paper

Until it is out, the code's claims live in docstrings rather than anywhere
citable.

**Not on this list, having been checked:** v1.0.0 *is* released on GitHub
(2026-08-13). An earlier note here saying otherwise was stale.
