# Plan: rebuild Figure 12 in Figure 11's shape

**Not started. This is a scoping document.** Figure 12 currently has two panels;
the target has ten. The work is data generation, not plotting.

## Where the two figures stand

|  | Fig. 11 (`fig:smooth_reach`) | Fig. 12 (`fig:shock_plane`) today |
|---|---|---|
| panels | 5 rows x 2 cols | 1 row x 2 cols |
| columns | exponential profile, PREM chord | shock width 1e-06, width 1e-03 |
| rows | profiles, then 2nu, 3nu, 3+1, 3+2 | one flavour count |
| swept | 12 energies at one baseline | 61 distances at one energy (15 MeV) |
| dial | rtol, both codes | n_slabs, both codes |
| Magnus series | orders 4, 6, 8 | order 4 |
| reference | DOP853, rescored against mpmath | frozen `solve_ivp` |

The shock data is keyed by `width` and carries no `flavours` field at all;
3+1 and NSI live in separate files (`external_shock_4nu.json`,
`external_shock_nsi.json`) with their own three- and four-point sweeps.

## Three decisions needed before any compute

1. **What is swept along the x-axis?** Figure 11's cost axis comes from sweeping
   *energies*; Figure 12's from sweeping *distances* along the ray. "The same
   format" admits either. This choice determines the reference, the cost, and
   whether the existing shock reference can be reused at all. Everything below
   assumes the 61-target distance sweep is kept, because that is what the shock
   physics is about; say if not.
2. **Does the thin front get `t_breakpoints`?** Figure 11's PREM column declares
   the sixteen layer crossings, and the caption says so. The shock analogue is the
   front itself. This is not a detail: a front of width 1e-06 is *exactly* the
   narrow-feature case that `HiddenFeatureWarning` exists for, where every engine
   is wrong together and refinement never helps because it never puts a point
   inside the feature. Without declared edges the Magnus column is measuring the
   blind spot rather than the method.
3. **The 3+2 row is Magnus-only**, as in Figure 11 -- NuOscProbExact has no
   five-flavour route. Confirm that is wanted rather than dropping the row.

## Phase 1 -- references (the long pole)

Eight references are needed, four flavour counts x two widths; two exist, at
three flavours. Mirror Figure 11's method exactly, or the two figures are not
comparable: a tight `DOP853` integration at `rtol = 1e-12`, `atol = 1e-14`, its
own convergence verified per configuration by tightening to `1e-13`, then
rescored against an mpmath reference (`gen_mp_reference.py`,
`rescore_against_mp_reference.py`). Shock analogues of both do not exist.

`make_shock_reference.py` freezes the current three-flavour `solve_ivp` ground
truth and is the place to start, but it is a different construction from Figure
11's and cannot simply be reused.

**Unknown cost, and the main risk to the schedule.** A DOP853 solution across a
shock front is stiff where the front is sharp.

## Phase 2 -- Magnus, orders 4, 6 and 8 by tolerance

120 points: 3 orders x 5 tolerances x 4 flavour counts x 2 widths.

**Cost is the problem.** At three flavours the tightest current point runs
554,565 us per probability; over 61 targets that is 34 s per call, and `timed`
takes the best of five, so **2.8 minutes for one point**. Four and five flavours
cost several times more. A conservative reading puts Phase 2 in the several-hour
range, and it should be run overnight rather than interactively.

**Risk: the ladder may ceiling instead of converging.** On a sharp front the
refinement can run to `max_n_slabs` rather than meet the tolerance -- which is
what Figure 11's 3+1 Earth panel shows, where `1e-3` and `1e-4` return the same
number and the points sit on top of each other. If that happens on the thin-front
column, a tolerance sweep says *less* than the slab sweep it replaced. Test one
column at two tolerances before committing to all 120 points.

## Phase 3 -- NuOscProbExact by tolerance

30 points, and cheap: its tightest current point is 2.4 ms per probability
against Magnus's 554. Route it through `slabs._n_for_tolerance`, as
`append_npe_rtol_prem.py` does -- that helper is geometry-agnostic, taking an
`evaluate(n)` callable and refining against self-convergence, so it works on a
shock ray without modification.

## Phase 4 -- the figure

Little work. Figure 11's cell already loops
`for col_key, axes, bench in zip(('exp', 'earth'), cols, (BENCH, PREM_BENCH))`,
so the two columns become the two widths by passing different benchmark files.
`SERIES_STYLE`, the axis handling, the unreachable-point filter and the
label-offset machinery all carry over. A new `RTOL_LABEL_OFFSETS` key is needed
per column, and the labels want dragging with `docs/dev/nudge_fig11_labels.py`.

## Phase 5 -- caption and prose

The current caption's argument -- that the width of the front decides, not the
flavour content -- survives and is strengthened by having four flavour rows to
show it on. Table 2's "Supernova shock" row and Section 8's discussion both
quote from this figure and would move with it.

## Two constraints that apply throughout

**Both codes must be timed in one session.** The machine has been measured
drifting 12-20% on overhead-bound work against a commensurability probe that
reported 2.5%, so Magnus and NuOscProbExact numbers from different sessions
cannot share a time axis. See the Figure 11 note above.

**Verify which engine answers, in both arms**, with
`oscprob._engine_probe(info=...)`. A shock profile is exactly the input where the
dispatcher's choice is not obvious.
