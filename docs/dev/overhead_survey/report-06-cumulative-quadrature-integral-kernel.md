# A compiled cumulative integral for the trapezoid/simpson path

**scipy's cumulative-Simpson machinery is 57% of the simpson-path call; the cumulative
integrals plus the integrand temporaries are ~60% of the trapezoid one. A complex
cumulative-integral kernel is worth an estimated ~2x on the simpson path and ~1.5x on
trapezoid.** The author's 7.5% commutator figure is confirmed from the other side: the
commutator is 23% under cProfile here, and the cost is indeed elsewhere.

## The cost, and how it was measured

cProfile of the PREM chord scan at d=3, order 4, n_tpts_per_slab=9, 1088 slabs, 12
energies (`scratchpad/survey/prof_quad.py`):

- **simpson**: 234 ms/call (GL at the same setting: 5.5 ms). Of it:
  `sp.integrate.cumulative_simpson` 134 ms (57%), `_commutator_batched` 55 ms (23%),
  `cumsum` inside scipy 19 ms, `ascontiguousarray` copies 11 ms.
- **trapezoid**: 147 ms/call. `_cumulative_integral` 54 ms (its own numpy expression:
  pairwise means, cumsum, concatenate-with-zeros), commutator 36 ms,
  `_magnus_terms_quadrature`'s own integrand builds (`-0.5*C3 + F1*(D12 + D21)` and
  friends, full-(...,m,d,d)-stack temporaries) 32 ms.

Two structural taxes in `_cumulative_integral` (src/magnus/magnus.py:722):

1. scipy's `cumulative_simpson` "silently discards the imaginary part", so the routine
   calls it **twice** (real, imag) and then pays `re + 1j*im` -- three extra full-array
   passes on top of scipy's own several temporaries per call.
2. The trapezoid branch materializes the pairwise means, the cumsum, a zeros block, and
   a concatenate -- five passes for what is one running sum.

## The change

A numba kernel in `magnus.py`:

```
_cumulative_integral_kernel(y, ds, out, method_id) -> None
    # y: (nB, m, d, d) complex contiguous; out: same shape;
    # cumulative trapezoid or Simpson along axis -3, complex arithmetic inline,
    # out[..., 0, :, :] = 0, reproducing scipy's composite formulas:
    #   equal-interval Simpson sub-integrals h/12*(5f0+8f1-f2) and
    #   h/12*(-f0+8f1+5f2), summed in ascending index order (= np.cumsum's order).
```

`_cumulative_integral` dispatches to it for complex128 contiguous stacks and keeps the
existing branches (including the m<3 and no-`cumulative_simpson` fallbacks) verbatim for
everything else. `_full_integral`'s simpson branch (43 ms total under simpson -- small)
can stay on scipy.

Expected: the integral machinery drops from 134 ms to an estimated 10-20 ms (one pass,
~10 flops/sample over ~1M complex samples/call), taking the simpson call from ~234 to
~110-120 ms. Estimated, not prototyped -- this path was profiled but no kernel was
built, unlike reports 01-03.

## Correctness

- The values must reproduce scipy's *composite formulas*, not just "a" Simpson rule.
  Reproducing them bit-for-bit is plausible (the sub-integral expressions and a
  sequential cumsum are simple enough to mirror exactly) but must be **verified against
  scipy on random complex stacks at exact equality first**; if any element differs, the
  change lands under the mpmath-quantified-shift protocol instead, like the commutator
  kernel did.
- One behavioural trap to preserve, not fix: the m<3/old-scipy fallback silently uses
  the trapezoid rule for the *cumulative* integrals while `_full_integral` stays
  Simpson, so results are scipy-version-dependent today. The kernel must not change
  which branch any input takes.
- No signature changes, no module globals.

## Tests that bear on it

- The cumulative-path order/convergence tests in test_magnus_expansion.py (the
  order-table docstring rates: cumulative order 4 -> global rate 6, etc.) are the
  behavioural pins; they tolerate rounding-level shifts by construction.
- A new exact-equality test kernel-vs-scipy on random complex stacks (m odd and even,
  m=3 boundary) is the load-bearing addition.

## Risks, and what would falsify it

- Medium implementation cost for a **non-default path**: 'gl' is the default and the
  recommended method; trapezoid/simpson matter for orders above 8 and for the
  paper's cumulative-order comparisons. If those workloads are rare, the ~2x may not be
  worth the new kernel's maintenance -- the author should weigh usage, not speed.
- Falsifier: if the kernel cannot match scipy exactly AND the mpmath-quantified shift
  comes out above the accepted class (unlikely -- the arithmetic is fixed-order sums),
  or if a re-profile after reports 02-03 land shows the commutator share ballooning
  instead (it will not: it is compiled already).
