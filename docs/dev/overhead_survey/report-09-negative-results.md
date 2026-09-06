# Negative results: where the remaining headroom is small or structural

Three well-supported "do not bother"s, so the next session does not re-derive them.

## 1. The refinement ladder's ~4x is structural, and near its theoretical floor

Measured today at ~4x the cost of the same accuracy at a fixed slab count. The
arithmetic: a ladder with growth factor g that certifies at level n_f has paid
n_f * (1 + 1/g + 1/g^2 + ...) ~= n_f * g/(g-1); and certification means n_f *agreed
with* level n_f/g, so the fixed count that would have sufficed is n_f/g. The overhead
is therefore **g^2/(g-1)**, which is minimized at g=2 where it equals exactly **4**,
and is 4.5 at the default g=1.5. The measured ~4x sits at the structural floor already;
no tuning of g can go below 4, and the phase seed already truncates the early rungs.

Under the constraint "same convergence decisions, same slab counts" there is nothing to
reclaim: the ladder *is* the decisions, GL nodes do not nest between levels (so no
sample reuse), and slab operators at different widths share nothing. The honest lever is
the one the docs already state -- a good n_slabs floor / t_breakpoints skips rungs --
plus making each rung cheaper (reports 01-03 multiply through the ladder unchanged).
What would reopen this: a decision to allow extrapolation between rungs or one-sided
certification, which is a *behaviour* change and out of scope by constraint 2.

## 2. The compiled kernels are near their floor at d=2-3; parity with NuOscProbExact
   needs the full multiply-out, and only there

After reports 02-03, the d=2 marginal decomposes as: ch2 exponential 0.022, fused GL4
0.014, ordered product 0.012, At build 0.018, framing 0.005, remainder ~0.013
us/slab/E, total 0.084 vs NuOscProbExact's 0.049. The three compiled kernels
(0.048 together) already run at 27-62 ns/matrix; there is no order-of-magnitude left in
any of them (the Jacobi report's 20-30% rotation-update note is the largest known
crumb). The remaining structural gap vs NPE is honest: order 4 samples the potential at
2 nodes per slab and composes an extra commutator; NPE samples once at the midpoint and
is order 2. Closing to parity at d=2/3 would need the single fused
build->Omega->exponentiate->compose kernel, i.e., the multiply-out across (order,
method, flavour) the author flagged -- judged **not worth it**: the measured remaining
overhead it could recover (At build + remainder, ~0.03 us/slab/E at d=2) buys ~1.3x
while freezing the engine's structure into one kernel per configuration. The exception
worth a second look is the narrow Variant B of report 02 (At built inside the GL4
kernel), which takes ~half of that 1.3x for one extra scan-only kernel.

## 3. The `_samples_identical` docstring claim is false, but fix it inside report 02

"array_equal short-circuits on the first differing element" -- numpy's array_equal does
not short-circuit (it reduces a full == array), and the measured 77 us/chunk at d=3 on
a smooth profile is full-comparison cost. The fix is already inside report 02's kernel
(an actual early-exit scan); if that report is declined, a standalone two-line pre-check
(compare a handful of elements before falling through to array_equal) still recovers
most of the 5-9% -- and either way the docstring sentence should go, since it currently
documents an optimization numpy never performed.

## Also examined, no action proposed

- **V evaluation and Vmat build** in the separable scan: 0.0003 and 0.003 us/slab/E --
  noise.
- **Probability extraction** (|U|^2 transpose): 0.0007-0.002 us/slab/E -- noise.
- **Threading**: deliberately excluded; 408f53c already records the cost/benefit
  without committing.
- **The At broadcast build on its own** (0.018-0.062 us/slab/E by flavour): only worth
  touching via report 02's Variant B; a standalone out=-buffer micro-fix measured (in
  kind) as too small to carry a change.
