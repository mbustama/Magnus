# Audit: the cumulative branch, 2026-09-05

`src/magnus/oscprob.py` carried a standing request: two latent defects had been
found within minutes of making the cumulative branch the default, and "two found
in minutes says the branch needs its own audit before it carries the default
traffic". This is that audit.

**Finding: nothing further.** The branch behaves correctly on both classes it has
previously failed in, its refusals are deliberate, and it is no longer thinly
tested. The sentence that asked for the audit had itself gone stale.

## The two known defect classes, re-probed

Both defects were interface failures rather than physics ones, so both were probed
directly rather than through accuracy.

**Shape.** The first defect returned `(1, d, d)` where `(d, d)` was due. Scalar
energy and scalar baseline, at d = 2, 3, 4, 5, on `cumulative=False` and
`cumulative=True`: all eight return `(d, d)`. Clean.

**Keyword forwarding.** The second raised `TypeError` on a `convergence_info`
keyword the branch rejected. Nine keywords -- `n_jobs`, `integration_method`,
`growth_factor_n_slabs`, `max_num_loops`, `min_n_slabs`, `max_n_slabs`,
`min_n_tpts_per_slab`, `max_n_tpts_per_slab`, `magnus_exp_order` -- all reach the
branch and return a probability. Clean.

## Three refusals, all deliberate

The audit's first pass reported twenty discrepancies. **Every one was the audit's
own construction**, and saying so is the point of writing this down: a list of red
lines looks like a finding.

1. **Differing energies with `cumulative=True`** raises: "cumulative=True scans
   baselines at one energy; the given energies differ. Baselines nest and energies
   do not, so there is nothing to reuse across energies." Correct, and
   `cumulative=True` is specified to raise rather than fall back.
2. **A supplied `t_slab_edges`** raises: "cumulative=True builds its own slab grid
   and cannot also honor t_slab_edges." Correct.
3. **Mismatched energy and baseline lengths** raises from the shared input check,
   on both settings. Correct; the audit passed three energies against four
   baselines.

The remaining rows compared **different grids**: at a fixed `n_slabs`, the
per-point path gives each baseline that many slabs over its own length while the
cumulative scan walks the whole profile once and snapshots. A difference there is
discretization, not disagreement. The suite's own test avoids this by handing the
per-point arm the scan's grid, and that comparison was made during the v1.0.12
work: 4.9e-15 to 2.2e-14 at d = 2 through 5, against a 1e-12 bar.

## Coverage

Not thin. 24 `cumulative=True` call sites across five test files, twenty tests
named for the branch, fourteen touching `convergence_info`. All 25 selected tests
pass. They assert routing (that it engages on a real scan and stands aside
otherwise), ordering, agreement with the per-point path on an identical grid,
every requested baseline landing on a slab edge, agreement with `expm` for a
constant Hamiltonian, probe strictness, and that the probe does not warn about
answers it discards.

## What remains open, and is not this

The exposure recorded in the same comment -- that the hybrid path declines to the
cumulative scan at a threshold where accuracy steps -- is a *routing* question,
not a correctness one, and is tracked at `FINDINGS_ROBUSTNESS_PROGRAMME.md` §13
with a stated surgical fix. This audit does not touch it.

GitHub issue #52, which proposes compiling this branch's quadrature, was filed
conditional on this audit. The condition is met.
