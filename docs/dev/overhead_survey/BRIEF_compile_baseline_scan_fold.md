# Brief: compile the baseline scan's running product (report 05)

Design and evidence: `report-05-cumulative-scan-compiled-fold.md` in this
directory.

## What this is, and what it is not

`_osc_prob_cumulative_scan` in `src/magnus/oscprob.py` serves the public
:func:`osc_prob_energy_baseline` -- every baseline in `L_out` from a single
traversal, because each requested answer is a prefix of the next. **This is not
the cumulative-quadrature engine** (trapezoid/Simpson, report 06), despite the
shared word. It is a used public API with its own tests, not an under-exercised
branch.

The Python running-product loop costs ~6 us/slab against ~1.25 us/slab for
building the operators it consumes. Replace it with a compiled fold.

## The pattern is established -- follow it

This is the fourth instance of the same fix: v1.0.6 (`_ordered_product_batched`),
report 04 (`_ordered_product_into`, currently uncommitted in the tree), and the
Gauss-Legendre Omega kernels. Read `docs/dev/overhead_survey/README.md` first for
what landed and at what gain, and **read report 08 before quoting any stacked
speedup** -- ratios against a shared endpoint do not multiply.

Check whether `_ordered_product_into` (added by report 04) already does what this
loop needs, or can with a small extension. A second kernel that duplicates it
would be worse than reusing one, and this loop differs mainly in that it must
*snapshot* the running product at requested baselines rather than only at the end.

## The invariant most at risk

`tests/test_oscprob.py::test_energy_batched_scan_matches_per_point` asserts the
batched scan reproduces the per-point path. Changing how the running product folds
is exactly what could break it. Note it currently runs at three flavours only and
asserts 1e-12 rather than bit-identity -- so **passing it is necessary, not
sufficient**. Check the invariant yourself at 2, 4 and 5 flavours too.

## Measure what a caller sees

Every stage figure in this survey has been far larger than what reaches a caller:
2.0-5.2x on a stage became 1.25x end to end, 4.7x became 1.11x. Report both, and
state the difference. Report 04's gain also turned out to depend strongly on slab
count -- ~1.0x at 512 slabs, 2.3x at 32768 -- so **sweep the slab count** and say
where the gain appears rather than quoting one number.

**Verify which engine answers.** Constructing a workload that reaches the intended
code path took four attempts on report 04: `separable`, `magnus` and `hybrid` each
took the request first. Use `oscprob._engine_probe(info=...)` and assert on the
engine actually used, in both arms.

## Bit-identity

Aim for it, and **establish up front whether it is achievable, before
implementing**. If the fold routes through BLAS in the old path, it is probably
not -- MKL's `zgemm` uses FMA and a `fastmath=False` kernel cannot reproduce its
ordering, which is what happened in report 04. If so, say so at the start of your
report, quantify the shift, and score the new fold against a high-precision
reference so we know which side is closer to exact. Preserve the fold's
association order regardless.

## Constraints

- **Surgical edits.** Do not rewrite files or reflow untouched code.
- **Update or add docstrings and in-line comments** the edit affects; new kernels
  get a real docstring in the style of `_gl4_omega_core` in `magnus.py`.
- **No signature changes**, no change to what any routine takes or returns.
- Do NOT touch `CHANGELOG.md` or `pyproject.toml`. Do NOT commit.
- Run `python3 -m pytest tests/test_oscprob.py tests/test_engines.py -q`. Not the
  full suite -- it takes 19 minutes and the owner runs it.
- Measure on an idle machine, arms interleaved, control drift reported.

## Report

What changed; whether `_ordered_product_into` was reusable; the fold's gain; the
gain a caller sees, swept over slab count; bit-identity established up front with
its evidence; and anything you could not make exact.
