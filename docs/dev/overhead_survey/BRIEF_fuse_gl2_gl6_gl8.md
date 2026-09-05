# Brief: fuse the order-2, order-6 and order-8 Gauss-Legendre Omega branches

## What already exists (read this first, it is the pattern)

`src/magnus/magnus.py`:

- `_commutator_batched_core` / `_commutator_batched`  (line ~310, ~336) -- compiled
  commutator, already used by orders 4, 6 and 8.
- `_gl4_omega_core` / `_gl4_omega_kernel`  (line ~385, ~436) -- the fused order-4
  kernel, landed today. **This is the template.** Read its docstring: it states the
  three things it fuses and why the result is bit-identical.
- `_magnus_gl`  (line ~1224) -- the dispatcher holding all four order branches. The
  order-4 dispatch block is at ~1261; copy its shape (guard, flatten, widths
  broadcast, reshape back).
- `_samples_identical`  (line ~1205) -- the constant-A test each branch runs before
  doing any work.

## The job

Fuse the remaining three branches, one kernel each, following the order-4 pattern.

### Order 6 (`if order <= 6:`)
Currently: two `_samples_identical` passes, then three commutators with ~19 array
temporaries around them (`a1`, `a2`, `a3`, `C1`, `C2` and the argument builds).

### Order 8 (final branch)
Currently: three `_samples_identical` passes, then `S1/S2/R1/R2`, `B0..B3`,
`a1..a4`, and six chained commutators `C1..C6` -- roughly 65 array temporaries.
`_GL4_W1/_GL4_W2/_GL4_V1/_GL4_V2` are module constants; inline or pass them, but the
arithmetic must associate exactly as the NumPy expression does.

### Order 2 (`if order <= 2:`)
Currently one ufunc call: `return h*An[..., 0, :, :]`. There is nothing to fuse --
one operation has no boundaries to remove. The only conceivable gain is that
`An[..., 0, :, :]` is a stride-2 view, so NumPy reads it with a buffered strided
loop where a kernel could read it directly.
**So: write it, measure it, and report the number. Keep it only if it actually
wins. A kernel that loses to NumPy is worse than no kernel.** Report either way --
the decision is the author's, not yours.

## Why this is a better target than order 4 was

NumPy evaluates each operation across *all* slabs before starting the next, so every
intermediate is a full `(nB, d, d)` array written to memory and read back. A fused
kernel finishes one slab completely before moving to the next, so the entire chain of
intermediates lives in a few `d x d` scratch buffers in L1. At order 8 that is ~65
streamed arrays collapsing to a handful. The saving should be proportionally larger
than the 1.4x order 4 gave.

numba detail: hoist the `d x d` scratch buffers **outside** the slab loop and reuse
them. Allocating inside the loop is a per-slab malloc and will eat the gain.

## Hard constraints

1. **Bit-identical output.** Not "agrees to 1e-14" -- exactly 0.0 difference against
   the current code, on every shape and flavour. This is the whole basis on which the
   order-4 kernel was accepted. To get it:
   - accumulate each commutator in the same interleaved order as
     `_commutator_batched_core` (`s += X[i,k]*Y[k,j] - Y[i,k]*X[k,j]`);
   - form scalar factors in the same association as the expression they replace;
   - preserve each `C_i` as a stored intermediate where NumPy stores one -- do not
     algebraically re-fold the chain, even where the algebra is valid.
2. **No change to what any routine takes or returns.** Signatures, argument names,
   return shapes and dtypes all stay. The optimisation is invisible to callers.
3. **No accuracy sacrificed** anywhere.
4. **Preserve the constant-A fast paths exactly**, including that order 6 tests
   `A1==A2 and A2==A3` and order 8 tests all three pairs. Fuse the scan with early
   exit on the first differing element, as the order-4 kernel does.
5. **Guards must match the order-4 dispatch**: kernel present, `dtype == complex128`,
   correct node count on `An.shape[-3]` (3 for order 6, 4 for order 8). Everything
   else falls back to the NumPy expression, which stays in place.
6. **Surgical edits only.** Do not rewrite files or reflow untouched code. Touch the
   lines that must change and no others.
7. **Update every docstring and in-line comment the edit affects**, and write real
   docstrings for the new kernels in the style of `_gl4_omega_core`'s -- what it
   fuses, why it is bit-identical, and that the pure-Python form exists only as
   compilation input.
8. Note in the docstrings that numba and non-numba paths may differ at the 1e-14
   level, consistent with what landed for the commutator kernel.

## Verification you must do and report

- Bit-identity (expect exactly 0.000e+00) against the kernel-disabled path across:
  several batch shapes, d = 2..5, scalar **and array** widths, the constant-A path,
  and the edge cases `(0, m, d, d)` and `(1, m, 1, 1)`.
  *Array widths matter*: `h` already carries two trailing singleton axes, so passing
  `h` where the kernel wants one scalar per slab silently works for scalars and
  fails for arrays. That bug was made once today; do not repeat it.
- Guard behaviour: complex64, float64, F-order and reversed-stride inputs must all
  fall back and still agree at 0.000e+00.
- Confirm the kernel is actually reached on the real `osc_prob` path at each order,
  and that the *other* orders' behaviour is unchanged.
- `python3 -m pytest tests/test_magnus_expansion.py tests/test_engines.py -q`.
- Speed: marginal us/slab, kernel off vs on, interleaved, for each order.

## Report

Say what you changed, the measured gain per order, the bit-identity evidence, and --
explicitly -- anything you were unable to make bit-identical and why. If order 2 does
not win, say so plainly and leave it unimplemented rather than shipping it.
