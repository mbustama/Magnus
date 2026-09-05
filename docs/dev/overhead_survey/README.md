# Per-call overhead survey, 2026-09-05

Where Mag(nu)s's remaining time goes, and what could be done about it, surveyed
under three standing constraints: **no accuracy sacrificed**, **no change to what
any routine takes or returns**, and surgical edits rather than rewrites.

Measurements are marginal us/slab/E fitted over two slab counts on the PREM chord
and the smooth exponential profile, arms alternated, first timing discarded, with
an interleaved control workload reported alongside so machine drift shows up
instead of hiding in a ratio. Where a report quotes a ratio, the control ratio is
quoted with it.

## Landed on `dev-overhead`

| Report | Change | Gain | Bit-identical |
|---|---|---|---|
| 02 | Fuse the order-4 GL `Omega` into one kernel, constant-sample check included | 2.0-5.2x on the `_magnus_gl` stage | yes, exactly 0.0 |
| 07 | `lru_cache` on `_passthrough_kwarg_names` (59 us of `inspect.signature` per call) | 11% of a single-point call's ~541 us fixed overhead | yes |

Earlier the same day, on the same branch: the composition-loop kernel
(`_ordered_product_batched`) and the commutator kernel (`_commutator_batched`),
which orders 4, 6 and 8 all use.

## Surveyed, not implemented

| Report | Change | Claimed gain | Bit-identical |
|---|---|---|---|
| 01 | Batched Jacobi eigensolver for the 4x4 and 5x5 exponential | 1.77x at 3+1, 1.38x at 3+2 | **no** -- worst probability shift 2.8e-12 |
| 03 | Fuse `_expm_stack`'s anti-Hermiticity framing into one compiled pass | ~4.7x on the framing stage | yes in measurement |
| 04 | Compile the interaction-picture engine's slab fold | ~2.1x on that engine | not assessed |
| 05 | Compile the cumulative baseline scan's running product | 2-4x on the scan at large slab counts | not assessed |
| 06 | A compiled cumulative integral for the trapezoid/simpson path | scipy's cumulative-Simpson is 57% of the simpson-path call | not assessed |

Report 01 is the one that needs a decision rather than a measurement: it is the
only item here that changes numbers. Its status is tracked in
`resources/paper/pending-edits.md`.

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
