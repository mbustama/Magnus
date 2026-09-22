# Handover: notebook 24's expansion-order measurement

**Written:** 2026-09-22, on branch `parked-defects`, after the notebook prose pass
(`ae9c845`…`754a722`).  **Parked at the author's request.**  Nothing below has been started.

**The one-line version.** Section 4 of notebook 24 claims to measure what
`magnus_exp_order` buys, and measures nothing: its referee describes a different Earth from
the code it referees, so all six orders report the same error.  Fixing it takes three
changes, not one, and the third is the one that is easy to skip.

---

## 0. Read this before touching the numbers

I put a wrong number into this notebook during the prose pass and it was caught on review.
The failure is worth understanding first, because the obvious repair reproduces it.

I wrote that order 2 → 4 is worth **a factor of 5600**, measured against "a Magnus reference
at 4000 slabs".  That reference is *the code under test at a finer discretization*.  It
shares the error it is meant to measure, the shared part cancels, and what survives is the
difference between two discretizations — not accuracy.  Measured properly:

| referee | order 2 error | order 4 error | implied step |
|---|---|---|---|
| Magnus at 4000 slabs (**wrong**) | 8.93e-07 | 1.597e-10 | 5591x |
| independent, composition-matched | 8.72e-07 | 4.018e-08 | ~22x, and see §3 |

The two agree about order 2 and differ 250-fold about order 4.  **That disagreement is the
tell**, and it was in my data before I wrote the sentence.  The claim is retracted in
`754a722`; `FINDINGS_NOTEBOOK_PROSE_PASS.md` §6 is the record.

**Rule for this work: a reference built from Magnus cannot measure Magnus's error.**

---

## 1. What is wrong, with the proof

`24_magnus_performance.ipynb` cell 13 fixes the slab count at 600 "so the order is the only
variable" and scores six orders against `dop853_earth`.  Every row prints **1.137e-01**.

The referee builds its potential with

```python
matter.vcc_func_from_rho_func(rho, density_matter_is_in_g_per_cm3=True)
```

— the **default** electron fraction of 0.5 — while `osc_prob_3nu_earth` takes `Y_e` from
PREM layer by layer.  The two integrate different Earths, and that gap swamps the
truncation.

Rebuilt both ways, at `COSTHZ = -0.9`, 12 energies from 1 to 10 GeV, `n_slabs = 600`:

| referee | order 2 | order 4 | order 6 |
|---|---|---|---|
| flat `Y_e = 0.5` (as committed) | 1.137e-01 | 1.137e-01 | 1.137e-01 |
| PREM per-layer `Y_e` | 6.764e-05 | 6.733e-05 | 6.733e-05 |

The first row reproduces the committed output exactly, which settles the cause.  **The
second row shows the composition fix alone is not enough** — the orders still do not
separate.

### PREM's electron fractions, since I got these wrong once

Four zones, from `earth.electron_fraction_func_prem`:

| zone | radius | `Y_e` |
|---|---|---|
| core | r ≤ 3480 km | 0.4656 |
| mantle | 3480 – 6346.6 | 0.4957 |
| crust | 6346.6 – 6368 | 0.4952 |
| ocean | r > 6368 | 0.5551 |

None is 0.5, which is the point.  (I had written "0.466 through the core and mantle, 0.555
in the crust": 0.466 is the core alone and 0.555 is the ocean.  Corrected in `754a722`.)

---

## 2. Why the composition fix alone leaves 6.7e-05

Because `oscprob` does **not** treat `V_CC` as linear in `Y_e`.  It derives the average
nucleon mass from the composition, layer by layer, via `r = (1 - Y_e)/Y_e`.  Notebook 28's
cell 93 says so in its own output, and `magnus-ye-matching-has-a-floor` records the same
thing from the other direction.

A referee must therefore pass **both**:

```python
matter.vcc_func_from_rho_func(
    rho,
    electron_fraction=ye,
    ratio_number_neutrons_to_protons=earth.neutron_to_proton_ratio_from_electron_fraction(ye),
    density_matter_is_in_g_per_cm3=True)
```

With that, at `n_slabs = 600`:

| order | error vs the independent referee |
|---|---|
| 1, 2 | 8.715e-07 |
| 3, 4, 5, 6 | 4.018e-08 |

The orders separate.  The pairs still collapse, as they must — `gl` gives three distinct
schemes wearing six names.

---

## 3. Why that is still not enough to publish a number

The referee above is a midpoint slab product.  Its **own** convergence, measured:

| refinement | movement |
|---|---|
| 2000 → 8000 slabs | 2.40e-06 |
| 8000 → 32000 slabs | 1.50e-07 |

So at 16 000 slabs its own uncertainty is around **1.2e-07** — *larger* than the 4.018e-08
it reports for orders 3 and up.  **Those rows are below its resolution.**  All that is
established is that 2 → 4 is worth at least about tenfold; 4 → 6 is unresolved.

What is wanted is the construction notebook 25 already uses: `prem_referee` +
`refereed()` in its cell 11 — a slab product on the *continuous* profile, **Richardson
extrapolated**, returning `(value, uncertainty)` so the floor is visible rather than
assumed.  Notebook 25 reports "its own residual discretization error: 4.34e-07" and then
says Magnus agrees to within it.  That is the shape to copy.

**Do not quote a step the referee cannot resolve.**  If the extrapolated referee lands at,
say, 1e-09, then 4 → 6 becomes quotable; if it lands at 1e-07, say so and quote only 2 → 4.

---

## 4. The timing side: an idle machine, and why the notebook cannot tell

Notebook 24 times in **six** cells (2, 4, 6, 10, 13, 17) through
`best_of(call, repeats=3)`, which keeps the fastest of three.

**It carries no interleaved control**, unlike notebooks 25 and 28.  A contaminated run
therefore cannot announce itself.

Measured — three consecutive runs of the committed notebook:

| quantity | three runs | spread |
|---|---|---|
| array vs loop | 3.48, 3.51, 3.36 | 4.5 % |
| palindrome, expensive scan | 2.66, 2.60, 2.63 | 2.3 % |
| palindrome, expensive point | 1.82, 1.78, 1.90 | 6.7 % |
| **palindrome, plain PREM** | 0.82, 0.79, 1.01 | **28 %** |
| **rtol 1e-5 relative cost** | 0.73, 1.07, 0.70 | **53 %** |

Quantities measured over hundreds of milliseconds are stable to a few percent.  Quantities
measured over single-digit milliseconds move by a third to a half **between consecutive runs
with nothing heavy alongside**.  Those two are exactly the rows the committed notebook
reports as **1.10x** and **2.14x**, and none of the three runs reproduces either.

Two conclusions:

1. **A rebuild needs an idle machine.**  The window is short — the notebook executes in
   **32 s** — so this is easy to arrange.
2. **`best_of(3)` is not enough resolution for a 2 ms operation.**  Either raise the repeat
   count for those rows, or report them as "no measurable difference", which is what the
   notebook's own prose already says in words ("About 1.00x") while its table says 1.10x,
   and what notebooks 18 and 19 now say after the prose pass.

---

## 5. The work, in order

1. **Replace `dop853_earth` in cell 13** with a segmented, Richardson-extrapolated referee
   returning `(value, uncertainty)`.  Model it on notebook 25 cell 11's `prem_referee` /
   `refereed`.  Pass `electron_fraction` **and** `ratio_number_neutrons_to_protons` per
   layer.  Segment at `earth.prem_layer_edges_along_chord(COSTHZ)` and subdivide *within*
   each segment, so no slab straddles a boundary.
2. **Validate the referee's own floor** before believing anything it says.  This is the step
   I skipped last time and it is how the 5600x happened.
3. **Add an interleaved control**, as notebook 25 does, and print it.
4. **Rewrite cell 14** — it currently explains why the cell cannot resolve the orders — and
   **cell 19**, which says "The PREM chord above cannot price it".
5. **Decide what to do with the two millisecond-scale rows** (§4.2).
6. **Rebuild on an idle machine**, then re-verify every timing-derived number in the prose:
   `~3x`, `1.7x`/`2.5x`, `1.68x`/`2.48x`, `~1.0x`, `costs ~2x`, `1.87x`, `1.86x`.

### Verification before committing

- `tests/test_notebooks_match_their_generator.py`, `tools/lint_notebook_cells.py`.
- The pairs must still collapse exactly (1=2, 3=4, 5=6): that is a property of `gl`, and if
  it breaks, the change broke something.
- The referee's stated uncertainty must be **below** every error it reports, or the
  corresponding rows are not measurements.
- Check `git status` for `fig/` churn and for figures that differ only by CreationDate.

---

## 6. Blast radius

- **The paper is untouched.**  `main.tex` sets `\graphicspath{{figs/}}`; all 39 of its
  figures resolve inside `resources/paper/figs/`, written by notebook 28 alone.  Notebook 24
  writes to `../fig/`.
- **Those three figures are not even tracked.**  `fig/expansion_order.pdf`,
  `expansion_unitarity.pdf` and `dispatch_vs_tolerance.pdf` are gitignored — only 2 of the
  72 files in `fig/` are tracked — and none is referenced from `docs/`, the README or the
  paper.
- **No paper data is re-measured.**  `MAGNUS_PAPER_RETIME` stays unset; `paper_figure_cache.json`
  is not opened; notebook 28 is not rebuilt.
- **What does change:** notebook 24's own committed outputs, including its six live timings.

## 7. Estimate

**1½ – 2 hours of working time**, plus about **30 seconds of machine time per rebuild**.
The machine is not the constraint.  The referee validation (§5.2) and the control (§5.3) are.

An earlier estimate of 45 minutes was for a two-part fix that does not work; see §2 and §3.

---

## 8. Still open elsewhere, for context

**Notebook 25 section 11** is two regenerations behind `external_shock_benchmarks.json`
(`ccdb6b2`, then `50e6e56`; the notebook's outputs date from `7188a4d`).  A rebuild would
improve the result — "22 times more accurate for three times the cost" becomes about 27
times for 1.3 times — but it re-times five cells whose numbers are paired against competitor
measurements frozen in August, which would make the comparison mixed-date.  **Recommendation
on record: leave it.**  Refresh it only by re-running the whole benchmark, both codes
interleaved.

**The engine count.**  The paper and `engines.rst` each say "six engines" and mean different
sets; `ENGINE_FAMILIES` registers eight.  Notebook 22 now names both.  Reconciling the paper
with the docs is the author's call.
