# Per-call overhead survey, 2026-09-05

Where Mag(nu)s's remaining time goes, and what could be done about it, surveyed
under three standing constraints: **no accuracy sacrificed**, **no change to what
any routine takes or returns**, and surgical edits rather than rewrites.

Measurements are marginal us/slab/E fitted over two slab counts on the PREM chord
and the smooth exponential profile, arms alternated, first timing discarded, with
an interleaved control workload reported alongside so machine drift shows up
instead of hiding in a ratio. Where a report quotes a ratio, the control ratio is
quoted with it.

## Status, updated 2026-09-05 after the session that ran it

Gains below are **measured end to end by the repository owner**, not the survey's
own figures, which were mostly stage-level and so mostly larger.  Where the two
differ the end-to-end number is the one a caller experiences.

### Landed

| Report | Change | Measured gain | Bit-identical | Commit |
|---|---|---|---|---|
| 02 | Fuse the order-4 GL `Omega` into one pass | 1.25x median, 1.52x at 2nu | yes, exactly 0.0 | `683a66c`, v1.0.8 |
| 02+ | The same for orders 6 and 8 (the report's "order-6 sibling", extended to 8) | 1.29x and 1.75x median | yes, exactly 0.0 | `220c1fb`, v1.0.8 |
| 07 | Cache the pass-through keyword names | 1.19-1.22x on a small call | yes | `683a66c`, v1.0.8 |
| 01 | Batched Jacobi eigensolver at 4x4 and 5x5 | 1.81-1.95x (4nu), 1.46-1.59x (5nu) | **no** -- 1.557e-12 worst shift | `334ab33`, v1.0.9 |

Report 01 was gated before it was implemented: `err_k <= 10*err_e` against `eigh`
across 364 cells, worst ratio 6.39x as shipped.  Two and three flavours are
untouched by it, measured at exactly 0.0 -- the control the change rests on.

Before the survey, on the same branch: the composition-loop kernel
(`_ordered_product_batched`, v1.0.6) and the commutator kernel
(`_commutator_batched`, v1.0.7), which orders 4, 6 and 8 all still use in their
NumPy fallbacks.

**Compound over all of it**, baseline `a55b8a4` to `334ab33`, measured across
separate checkouts rather than by multiplying the rows above (see report 08 for
why multiplying is wrong): **11-12x at 2nu, 8.6-11x at 3nu, 2.8x at 4nu, 1.9-2.1x
at 5nu** on the default order-4 path.

### Declined after measurement

| Item | Why |
|---|---|
| Order-2 `Omega` fusion | Prototyped and measured: 1.4-1.8x marginal, but 0.82-1.04x *per call* at two flavours.  One ufunc has no boundaries to fuse, and the dispatch's fixed cost exceeds the streaming saving on a branch costing 0.01-0.05 us/slab. |

### Surveyed, not started

| Report | Change | Claimed gain | Bit-identical |
|---|---|---|---|
| 03 | Fuse `_expm_stack`'s anti-Hermiticity framing into one compiled pass | ~4.7x on the framing stage (12% of a d=3 call, 19% at d=2) | yes in measurement |
| 04 | Compile the interaction-picture engine's slab fold | ~2.1x on that engine | not assessed |
| 05 | Compile the cumulative baseline scan's running product | 2-4x on the scan at large slab counts | not assessed |
| 06 | A compiled cumulative integral for the trapezoid/simpson path | scipy's cumulative-Simpson is 57% of the simpson-path call | not assessed |

Report 03 is the closest thing to a next step: it is the only untouched item whose
report claims bit-identity, and it sits on the same hot path the landed work has
been shortening.  Reports 04-06 are each confined to one engine that the default
path does not use.

### Ruled out, with reasons (report 09)

- The refinement ladder's ~4x is **structural** and near its theoretical floor.
- The compiled kernels are near their floor at d = 2-3; parity with
  NuOscProbExact there is not reachable by more of the same.
- The false `_samples_identical` docstring claim -- fixed, in `683a66c`.

## Read these two first

- **09 -- negative results.** Three well-supported "do not bother"s, recorded so
  the next session does not re-derive them.
- **08 -- why the ratios did not multiply.** 3.49 x 2.1 != 3.66 because speedups
  measured against a shared endpoint do not compose; *time saved* adds. A 2x2
  factorial shows the two kernels' savings additive to within noise. Worth
  reading before quoting any stacked figure.

## Also here

- `BRIEF_fuse_gl2_gl6_gl8.md` -- the specification for extending report 02's
  fusion to the order-2, order-6 and order-8 branches, including the acceptance
  bar (bit-identity, not agreement to 1e-14) and the array-widths trap that bit
  once already.
- `prototypes/` -- working scripts behind the numbers: `jacobi_proto5.py` is
  report 01's final kernel, `stages.py` and `prof*.py` the stage splits,
  `*_e2e.py` the end-to-end A/B harnesses, `factorial.py` report 08's design.
  The `om_d4.npy` / `om_d5.npy` operator caches these read were **not** copied
  (8.7 MB, regenerable by `grab_om.py`).

Provenance: produced by subagents in one session and copied out of a scratchpad
that did not survive it. The gains in the second table are *claimed*, not
independently reproduced -- treat them as leads with a stated method, not as
results.
