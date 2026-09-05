# Compile the interaction-picture engine's slab fold

**The Python `@` loop is 56% of the IP engine's pass; a compiled accumulator-in fold is
worth ~2.1x on the engine end-to-end** (estimated from the measured stage split; the fold
itself shrinks from 427 ms to an estimated 30-60 ms per pass).

## The cost, and how it was measured

`_osc_prob_ip_exp_core` (src/magnus/oscprob.py:4628; the loop at ~4856) folds each
block's slab operators into a cross-block accumulator one Python iteration at a time:

```
for k in range(U_slab.shape[1] - 1, -1, -1):
    acc = U_slab[:, k] if acc is None else acc @ U_slab[:, k]
```

Stage split of one fixed-slab pass (2-level, exponential profile, 64 energies, 32768
slabs, e_chunk=64, blk=32, 1024 blocks; `scratchpad/survey/ip_stages.py`, stages sum to
749 of a measured 757 ms):

| stage                         | per pass  | share |
|-------------------------------|-----------|-------|
| fold (python @ loop)          | 426.6 ms  | 56%   |
| arg + exp + I                 | 136.4 ms  | 18%   |
| `_expm_stack`                 |  78.3 ms  | 10%   |
| U_free diagonal exp           |  60.5 ms  |  8%   |
| Omega build + max             |  45.4 ms  |  6%   |
| V0 eval                       |   2.2 ms  |  0.3% |

Each fold iteration is a (64, 2, 2) x (64, 2, 2) batched matmul costing ~13 us through
the gufunc machinery -- ~30x its arithmetic. The engine's per-slab marginal is 0.36
us/slab/E, which is *the same as the full 3nu GL scan* for what is a 2x2 elementwise
method; that is all dispatch.

## The change

A numba kernel in `magnus.py` beside `_ordered_product_batched`:

```
_ordered_product_into(acc, U) -> None
    # acc: (nE, d, d), updated in place; U: (nE, nb, d, d)
    # for k in nb-1 .. 0:  acc <- acc @ U[:, k]
    # accumulator on the LEFT, k descending -- the same association the
    # Python loop uses, element sums accumulated the obvious triple-loop way.
```

The `acc is None` first block seeds `acc = U_slab[:, -1].copy()` and folds the rest.
Everything else in the engine -- tiling, block order, elementwise builds, `_expm_stack`
-- is untouched. The cross-block accumulator semantics are preserved exactly because the
accumulator is an *argument*, not recomputed per block: this is why the kernel must take
`acc` in rather than reduce each block and multiply afterwards, which would change the
parenthesization the tiling comment explicitly forbids.

## Correctness

- Association identical; per-element rounding may move at the 1e-14 class because
  numba's triple-loop dot need not match numpy's small-matmul kernel bit for bit --
  the same situation as `_ordered_product_batched` (bc393b8), which was accepted at a
  measured worst shift of 1.28e-14. Quantify the same way (mpmath references) before
  landing; expect the same class.
- **The exact-equality tiling test is the gate**: tests/test_oscprob.py:1829 pins a
  tiled run against an untiled one with `np.array_equal`. Since both runs would use the
  same kernel, tiled-vs-untiled stays exactly equal -- the test verifies the tiling, not
  the matmul backend -- but run it first; if it fails, the accumulator plumbing is wrong.
- Convergence decisions unchanged: `max_omega`, the trust threshold, the certification
  bound, and the agreement logic never touch the fold.
- No signature changes, no module globals.

## Secondary opportunity in the same engine (measure before bothering)

The three elementwise stages (arg+exp+I, U_free exp, Omega build; 32% together) allocate
~10 temporaries per block and call complex `np.exp` twice per element class. A single
elementwise numba kernel per block could roughly halve them (~16% of the pass), with the
caveat that complex exp dominates the arithmetic and numba's is no faster than numpy's.
Worth doing only if the fold lands first and the engine still matters on the profile;
the fold is the story.

## Tests that bear on it

- tests/test_oscprob.py:1829 (tiled == untiled, exact) -- the load-bearing pin.
- The IP dispatch/certification tests in test_engines.py (trust threshold, give-up
  exits) exercise decision paths that must not move.
- The solar-average benchmark (notebook 25, comparison.rst) is the user-visible number
  this engine feeds; re-measure it after landing (~0.7 s today).

## Risks, and what would falsify it

- Modest. The fold touches the accumulator that the memory-bug fix (BUG_IP_EXP_MEMORY)
  was built around; keep `live_arrays=8` bookkeeping intact -- the kernel adds no
  same-shape temporary, so the budget is unchanged.
- Falsifier: a probability shift beyond the ~1e-14 class on the mpmath-referenced
  profiles, or any change in which calls certify vs give up on the certification test
  battery. Either kills the "no behavioural change" claim.
