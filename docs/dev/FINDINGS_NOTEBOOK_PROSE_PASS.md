# Findings: the notebook prose pass

**Written:** 2026-09-22, on branch `parked-defects`, commits `ae9c845` … `ccdf4db`.

Every notebook's prose read against the cells beneath it. **111 markdown cells** rewritten
across 26 notebooks, **2 code cells** (both notebook 21), 31 files, +3201/−2976. Notebooks
19, 24, 25, 27 and 28 were never executed — their markdown was patched into the generator
and the `.ipynb` in lockstep — so every live timing survived.

---

## 0. The finding that reordered the work

**A notebook can match its generator perfectly and still be wrong.** The freshness guard
(`tests/test_notebooks_match_their_generator.py`) compares cell *sources* and never looks at
outputs. Nothing in CI can see prose describing numbers the package has since moved past.

Four notebooks were stale:

| notebook | what had moved underneath it | size of the move |
|---|---|---|
| 14 | a `src` fix had restored the `t_breakpoints` cure | 8.257e-01 → 4.307e-06, **51 680x** |
| 04 | `946769f` fixed the sign of West/South coordinates | Fermilab–Homestake 1207.2 → **1284.7 km** |
| 21, 22 | this branch's `SUN_RADIUS` correction, 6.947e5 → 6.957e5 km | ground truth moved by **0.056** |

`ac1811a`'s own commit message predicted this: *"A notebook can run clean while its prose
asserts something its new figures contradict."*

**Working rule adopted after notebook 14:** rebuild (or script-verify, for the seven that
time live) *before* writing prose, never after.

---

## 1. The defects, by kind

### Numbers that drifted from the cell that produces them

The commonest kind. A sample:

- **16**: summary table said `1e-14` where the cell prints `4.84e-14`, and `6e-16` where it
  prints `3.33e-16`.
- **26**: first release's 68 % width given as 0.024; the cell prints **0.035**. Narrowest
  release named as v4.0; the cell names **v4.1**.
- **24**: "pass an array of energies ~2.7x" against a printed **2.99x**; the palindrome
  "~1.8x" against a printed **1.68x** at a point and **2.48x** on a scan.
- **17**: "an order of magnitude smaller than the ordering effect" for a measured **factor
  of thirty** — understating the notebook's own headline result by three.
- **13**: the Sun–Earth path given as `~1e10` cycles; it is **~1e6** at 5 MeV.

### Internal contradictions

The most serious kind, because each one is visible without leaving the notebook.

- **23** attributed a 53x averaging suppression to notebook 13, which prints **0.84x** and
  whose entire headline is that averaging does *not* help there — as notebook 23's own
  closing paragraph said. The numbers are real and belong to
  `docs/dev/adversarial_batteries/p_avg_check.txt` (cubic row), the same ray with a
  different window. Re-attributed, and notebook 13's 0.84x / 1.47x now named in the closing.
- **29** calls a splitting "still coherent" in §5 and "NEITHER limit" in §6. Measured: at
  `dm2 = 2.6e-17` over 100 Mpc, `coherence_report` returns undecided for all three pairs, at
  2.03 rad. The *demonstration* is sound — `block_form` sums amplitudes with no phase
  factor, so the factor of two is algebra — but the claim of coherence was not.
- **27** said six of nine scenes animate and then listed seven; seven is also what the code
  writes, what `img/` holds, and what its own render table shows.
- **22** said "six engines, five families" directly above a print of **eight**. Both counts
  are right: the paper lists dispatch *rows*, where `constant` is the constant-density case
  of `separable` and `expm` is the ladder's first term. Now reconciled rather than picked.

### Measurements that had stopped measuring

- **24 §4** fixes the slab count "so the order is the only variable" and prints the **same
  error for all six orders**. See §3 below.

### Dead pointers

- **All 29 footers** linked to `implementation_details.html`, split into `engines.rst`,
  `performance.rst` and `diagnostics.rst` by `fee6ec0`. Repointed at the docs home, which
  survives the next split; `add_footers()` carries a comment saying why.
- **20**'s warning table listed **nine of fourteen** classes and cited the same dead file
  for false-alarm rates. Five rows added from each class's docstring; rates now quoted from
  `diagnostics.rst` (76 %, 59 %, 57 %, against 2 silent misses in 168 configurations).
- **16** and **20** both pointed at notebook 21 for warning material; notebook 21 contains
  **zero** mentions of any warning class.

### One error propagated between notebooks

**25 §12's phase table divides by 2π where it should divide by π** — one oscillation length
advances the sin-squared argument by π. Every entry was half: 2 352 → **4 705**, 940 981 →
**1 881 960**, 9 410 → **18 820**. Notebook 14 had inherited it as "some 940 000
oscillations". Both fixed; notebook 14's cell 9 prints **4 726** for the first row, which
the table now cites as an independent check.

---

## 2. Claims checked and found correct

About as much of the work as the fixes, and worth recording so nobody re-checks them:

- **`4 x 14 = 56` is every wrapper there is** — exactly 14 scenario stems exist at all four
  flavour counts; the two extras in the module (`matter_std`, `vacuum_std`) are the
  `oscprobstd` closed forms, outside the product set.
- **"Every wrapper takes `nubar=True` and does both"** — worth checking because **none of
  the 67 `osc_prob*` callables declares `nubar`**; it arrives through `**kwargs`, which is
  how a silently-ignored flag hides. Probed vacuum, constant-density, NSI and solar at 2, 3
  and 4 flavours: every one honours it. The single zero difference, 2ν vacuum, is physics.
- **"the sign is applied in exactly one place"** — `matter.py:680` is the *only* `nubar`
  sign-flip site in the package, and `VCC_func` takes no `nubar` at all.
- **"measured false about three quarters of the time"** — nearly deleted as unsourced;
  `diagnostics.rst:304` measures it at **76 %** (70 firings, 17 true, 53 false).
- **27's two size tables** match `img/anim_*.gif` on disk file by file, 14.56 MB over 7 files.
- **28's "the two routes agree to between 5 and 8e-5"** — the four measured rows give
  5.076e-05, 6.024e-05, 7.083e-05, 7.924e-05.
- Nine PREM boundaries; the exponential fit 2.4x too dense at the centre; `USE_PALINDROME`
  at `magnus.py:2772`; the partial-chord palindrome withdrawal at `magnus.py:2888`.

**Two numbers I was about to "correct" and should not have:** the three-quarters figure
above, and notebook 25's `9.9e-10`, which is the 128 000-slab sweep point rather than the
printed best. Both times, reading the source before editing caught it.

**One correction I made and had to undo.** I removed the "18 shock configurations, improved
7, worsened 11" tally from notebook 14 as unreproducible. It is stated in four places —
`oscprob.py:5745` (inside the warning the package emits), `recipes.rst:309` (which *links to
notebook 14* and calls it "that measurement"), notebook 18, and notebook 14. Deleting it
from the one place the other three point at broke all three. Restored, with the worked case
kept alongside. **Grep the repo for a claim before deleting it from one file.**

---

## 3. Open, and why

### Notebook 24 §4 — broken, needs two code fixes

`vcc_prem_at` builds its DOP853 referee with the **default electron fraction of 0.5**, while
`osc_prob_3nu_earth` takes Y_e from PREM layer by layer — 0.466 through core and mantle,
0.555 in the crust. Proven by rebuilding the referee both ways:

| referee | order 2 | order 4 | order 6 |
|---|---|---|---|
| flat Y_e = 0.5 (as committed) | 1.137e-01 | 1.137e-01 | 1.137e-01 |
| PREM per-layer Y_e | 6.764e-05 | 6.733e-05 | 6.733e-05 |

The first row reproduces the committed output exactly.

**Correcting Y_e is not sufficient.** The referee still floors at 6.7e-05 and cannot
separate the orders, because an adaptive DOP853 crosses the ~16 PREM density jumps without
being told they are there. It needs to integrate piecewise *between* the boundaries.

**Retracted — see §6.** I wrote here and in the notebook that an independent measurement
gives order 2 → 4 = 5591x. It does not: that figure came from a Magnus reference judging
Magnus, and an independent referee puts order 4's error 250 times higher. Both the notebook
claim and this one are withdrawn.

Cost: notebook 24 times in seven cells through its `best_of` helper. `fig/expansion_order.pdf`
is a docs figure, **not** a paper figure — the paper draws from `resources/paper/figs/`.

### Notebook 25 §11 — not broken, two regenerations behind

`external_shock_benchmarks.json` was rewritten by `ccdb6b2` (09-06) and `50e6e56` (09-11);
the notebook's outputs date from `7188a4d` (09-05). The JSON as of `ccdb6b2^` reproduces
cell 60's printed numbers exactly, which settles the provenance.

| 70 km front | notebook prints | file says today |
|---|---|---|
| Magnus best | 9.600e-10 at 576 410 µs | **1.925e-11 at 178 426 µs** |
| NuOscProbExact best | 5.600e-06 at 2 461 µs | 5.060e-06 at 2 320 µs |

A rebuild would move the result **in Magnus's favour**: "22 times more accurate for three
times the cost" becomes about **27 times for 1.3 times**, and the tail 9.9e-10 → 1.9e-11.

**But a rebuild is worse than the staleness.** Notebook 25 times Magnus live in five cells
and compares it against competitor numbers frozen in `external_*.json` — cell 21 prints
"frozen dataset … our chord at that costhz". Re-timing Magnus today while the competitors
stay frozen from August destroys the pairing, which is the contamination
`magnus-per-call-overhead` records. §11's own numbers already come from a file where both
codes were measured together under an interleaved control (0.996).

`external_profile_benchmarks.json` also moved after the run, but immaterially
(2.899e-13 → 2.863e-13 on the quoted row). The other six data files predate the run.

### Not covered

- **Notebook 01** never had the pass; only its footer was touched.
- **Notebooks 10 and 12** parked by the author. Notebook 12's deferred caveat (the hybrid
  route's self-certification blind spot, with a pointer to notebook 22) is markdown-only and
  can land through `prose_only.py` with no rebuild and no re-timing.

---

## 4. The instrument

`scratchpad/prose_only.py` patches the generator and the `.ipynb` markdown in lockstep
**without executing anything**, exploiting the fact that the freshness guard compares
sources. Two things it must get right: call `mk.add_footers()` before reading
`mk.books[name]`, because `build()` does and otherwise the built notebook is one cell short
of disk; and write with `json.dumps(..., indent=1, ensure_ascii=False)`, which round-trips
these files byte-for-byte. It refuses any edit that changes the cell count or touches a code
cell.

This is what made notebooks 19, 24, 25, 27 and 28 reachable at all.


---

## 5. Second pass, same day: sweeps rather than a second read

Commit `18eb86e`.  Pass 1 read notebook by notebook; pass 2 built a check for each defect
class pass 1 had turned up and ran it over all twenty-nine at once.

### Found

- **Notebook 01 had never been read** -- only its footer was ever touched, and it is 74
  markdown cells and 3900 words.  Six typos, two broken sentences, an unbalanced
  parenthesis, a cross-reference to `10_magnus_matrix_exponential.ipynb` (it is **11**; 10
  is the averaged probability), and "the NuFit 6.0 global fit" where cell 2 loads
  **NuFIT 6.1**.
- **It taught two settings that do nothing.**  `n_jobs = 10` in four cells, with prose
  claiming it pays off "as here, with a thousand baselines" -- those baselines are a Python
  loop calling `osc_prob` once each, and `oscprob.py:6359` gates parallelism on
  `n_points > 1`.  Inert.  Had the scan been batched it would have been *worse* than inert:
  `oscprob.py:4925` makes the batched engines decline whenever `n_jobs != 1`.  Alongside
  it, `magnus_exp_order = 3` "for speed", where 3 and 4 share the two-node Gauss--Legendre
  scheme (notebook 24: 9.67 ms against 10.01 ms).  **Both removed, and all ten figures came
  back byte-identical** -- which is the proof rather than the argument.
- **Two stale notebooks**, both now current: 01 in its version banner only, and **27's
  section 7** (0.1817 -> 0.2120 and two others) from the same `SUN_RADIUS` correction that
  moved 21 and 22.  No prose quoted those three.  Rebuilding 27 regenerates no GIFs:
  `RENDER = False`, every GIF write is behind it, there is no `savefig` outside that guard,
  and `img/` stays untouched.
- **A miss from pass 1:** notebook 19's footer still said "nine warnings" after notebook
  20's table became fourteen.  The blurb lives in `READING_ORDER` and renders in the
  *previous* notebook -- the same trap as notebook 22's blurb.  The stale count had reached
  `tutorials.rst` and `installation.rst` too; the latter's tree is generated from
  `tests/test_file_tree.py`, so the fix went to the source and the tree was regenerated.
  `tutorials.rst` was also still quoting notebook 21's pre-solar-radius `2.5e-2`.
- Eight literal `Magνs` in the two parked notebooks, one "colours", one micro sign.

### The engine count: two sixes, different sets

Worth naming rather than smoothing over, and it predates this work:

| source | the six it counts | omits |
|---|---|---|
| paper, Fig. 1b and `tab:engines` | average, hybrid, ip_exp, separable, cumulative, magnus | `constant`, `expm` |
| `docs/source/engines.rst` | constant, hybrid, ip_exp, separable, cumulative, magnus | `average`, `expm` |
| `oscprob.ENGINE_FAMILIES` | all eight | -- |

`engines.rst` also states that `expm` "is not an engine but is used as an oracle", which
pass 1's notebook 22 wording contradicted.  Notebook 22 now names both sixes and says which
omits what.  **Reconciling the paper with the docs is the author's call.**

### Sweeps that found nothing, so they need not be repeated

Broken notebook links (0); repo paths named in prose (1 hit, a false positive); per-notebook
number check against each notebook's own outputs and data files (17 candidates, all sourced
in code comments or set during pass 1); near-duplicate sentences across notebooks (5, all
deliberate restatements carrying no numbers); unbalanced inline math (0); reST roles left in
markdown (0, after pass 1 removed the one `:doc:`); non-ASCII outside the expected set (the
`ö` of Schrodinger, and the one micro sign).

### A process error worth recording

Removing those two keywords, I ran `str.replace` over the **whole generator** rather than
notebook 01's block: 24 sites across several notebooks instead of 8.  `git diff --stat`
caught it before any rebuild.  **Scope every replacement to the notebook block it belongs
to** -- the generator is one file holding twenty-nine documents.


---

## 6. Deep audit of the numbers this work introduced, 2026-09-22

Prompted by a justified challenge.  I had put "a factor of **5600**" into notebook 24's
committed prose, derived from a self-refereed comparison -- in the same paragraph where I
withheld the 4 -> 6 step *because* it was self-refereed.

### Method

The 151 markdown lines added across `ae9c845..18eb86e` that carry a number, sorted by
**provenance**, since provenance predicts this failure:

- transcribed from a cell's own printed output;
- read from a sourced file in the repository;
- computed from table data or arithmetic, independently checkable;
- **measured by comparing the package against a reference I built** -- the only class that
  can fail the way 5600x failed.

Exactly two claims fell in the last class: notebook 16's convergence statement, and
notebook 24's factor.

### Results

| claim | verdict |
|---|---|
| nb12 2.9e-02 blind spot | stated verbatim in `adiabatic.py` |
| nb13 X range, 37 % / 76 %, fit 1.6x and 2.4x | recomputed from the BS05 table |
| nb14 70 m front | arithmetic |
| nb14 resonance 30 323 / 28 058 km | re-derived; one crossing, as claimed |
| nb15 collapse above 6 GeV, 0.03 vs 0.02, peaks 2.8 / 2.1 / 2.4 | all reproduced |
| nb16 "fifth decimal ... sixth" | **audited against an independent referee -- sound** |
| nb17 0.015, 0.49 / 0.42, thirty times | reproduced (29.9) |
| nb19 ratio tends to twice the contrast | reproduced; surface term heads for 0.5 |
| nb24 Y_e "0.466 core and mantle, 0.555 crust" | **WRONG -- corrected** |
| nb24 "a factor of 5600" | **WRONG -- removed** |
| nb25 phase table 4 705 / 1 881 960 / 18 820 | two derivations agree |
| nb26 0.039 between 1.1 and 1.2 | from the printed table |

Notebook 16 is worth singling out: its self-movement under a tighter tolerance is 3.39e-05
and its true error against an independent referee is 3.40e-05.  There the ladder's own
estimate is honest, and the claim survives -- now checked in a way the notebook itself
cannot check.

### The two failures

**The 5600x.**  Order 4 at 600 slabs measures 1.597e-10 against a Magnus reference at 4000
slabs and **4.018e-08** against an independent one.  The two referees agree about order 2
(8.93e-07 and 8.72e-07) and differ 250-fold about order 4, because the Magnus reference
shares order 4's error and it cancels in the difference.  Removed, and nothing replaces it:
the honest statement is that the cell cannot price the order, which the notebook now says.

**The electron fractions.**  PREM has **four** Y_e zones, not two: core 0.4656 (r <= 3480),
mantle 0.4957 (to 6346.6), crust 0.4952 (to 6368), ocean 0.5551 above.  I wrote "0.466
through the core and mantle, 0.555 in the crust" -- 0.466 is the core alone, and 0.555 is
the *ocean*.  Corrected.  **The same wrong pair is in commit `c4d6b9a`'s message**, which
cannot be rewritten; this is the correction of record.

### What fixing notebook 24 actually requires

Established by measurement during the audit, and more than §3 first proposed:

1. per-layer `electron_fraction`, **and**
2. per-layer `ratio_number_neutrons_to_protons` = (1-Y_e)/Y_e -- without it a residual of
   6.7e-05 survives and the orders still do not separate; and
3. a **Richardson-extrapolated** referee that reports its own uncertainty.  A plain midpoint
   product self-converges to only ~1.2e-07 at 16 000 slabs, which cannot resolve orders 3
   and up, sitting at 4.0e-08.

With (1) and (2) alone the notebook could say 2 -> 4 is worth at least about tenfold, and no
more.  Only (3) prices 4 -> 6.

### Reading for next time

The rule I stated and then broke: **a reference built from the code under test cannot
measure that code's error.**  It measures the difference between their discretizations,
which is smaller -- here by 250x -- and looks like accuracy.  The tell was in the data and I
did not take it: the two referees agreed about order 2 and disagreed about order 4, which is
the signature of a shared error cancelling.
