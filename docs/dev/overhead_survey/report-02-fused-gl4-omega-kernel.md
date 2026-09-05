# Fuse the order-4 GL Omega into one kernel (constant-sample check included)

**Bit-identical** (|dP| exactly 0.0 end-to-end at every slab count and flavour tested),
and worth 2.0-5.2x on the `_magnus_gl` stage. Together with report 03 it delivers the
measured end-to-end **1.42x at 2nu and 1.41x at 3nu** (their separate shares below).

## The cost, and how it was measured

At order 4 (the default), `_magnus_gl` (src/magnus/magnus.py:1167) runs four steps whose
cost is mostly numpy overhead, not arithmetic (stage timings, PREM chord, order-4 GL,
1104 slabs x 12 energies, us/slab/E; `scratchpad/survey/stages.py`):

| stage                                   | d=2    | d=3    | d=4    | d=5    |
|-----------------------------------------|--------|--------|--------|--------|
| `_samples_identical(A1, A2)`            | 0.0138 | 0.0231 | 0.0352 | 0.0534 |
| `_commutator_batched(A2, A1)`           | 0.0288 | 0.0616 | 0.1274 | 0.2341 |
| linear combination -> Omega             | 0.0296 | 0.0491 | 0.0807 | 0.1272 |
| whole `_magnus_gl`                      | 0.0733 | 0.1331 | 0.2417 | 0.4095 |
| (whole call, for scale)                 | 0.152  | 0.409  | 2.55   | 3.51   |

Three separate inefficiencies add up:

1. **`_samples_identical`'s docstring is wrong.** It claims "`array_equal` short-circuits
   on the first differing element". numpy's `array_equal` does not short-circuit: it
   builds the full `==` array and reduces it, which is why the check costs 77 us per
   chunk at d=3 on a smooth profile where the very first element already differs. On
   PREM this pure-overhead check is 6% of the whole call at d=3, 9% at d=2.
2. `A1 = An[..., 0, :, :]` and `A2 = An[..., 1, :, :]` are **strided views**, so
   `_commutator_batched` pays `ascontiguousarray` copies of both, and the linear
   combination gathers strided memory through four temporaries.
3. The linear combination `0.5*h*(A1+A2) + (sqrt(3)/12)*h*h*C` allocates ~4 full-stack
   temporaries for ~6 flops per element.

## The change

One numba kernel in `magnus.py`, general in d (loop bounds from the shape):

```
_gl4_omega_from_At(At, h, out) -> int
    # At: (nB, 2, d, d) contiguous;  h: (nB,) float;  out: (nB, d, d)
    # 1. early-exit equality scan over At[:,0] vs At[:,1]; on the first
    #    difference, jump to step 2.  If ALL equal: out = h*A1, return 1
    #    (the constant-A fast path, exactly as today).
    # 2. out[b] = (0.5*h[b])*(A1+A2) + ((sqrt(3)/12)*h[b]*h[b])*(A2@A1 - A1@A2)
    #    with the commutator accumulated interleaved, exactly as
    #    _commutator_batched_core does.
```

Dispatched from `_magnus_gl`'s `order <= 4` branch when the input is a complex128
C-contiguous stack with m=2 (the GL path always is; anything else falls through to the
existing expression, as `_commutator_batched` already does). The leading axes flatten to
nB and h broadcasts to it -- ~15 lines of dispatch, mirroring `_commutator_batched`'s.
The kernel's return value states which path it took, replacing `_samples_identical` on
this branch (the helper itself stays for the order-6/8 branches).

Working prototype: `scratchpad/survey/fusion_proto2.py`; A/B harness
`fusion_speed2.py`; end-to-end `fused_e2e.py`.

Measured on the real memory layout (strided-view baseline, i.e., what `_magnus_gl`
actually pays): 798 -> 152 us (5.2x) at d=2, 443 -> 186 us (2.4x) at d=3, 267 -> 136 us
(2.0x) at d=4, 438 -> 251 us (1.7x) at d=5, per chunk.

## Correctness

- **Bit-identical, by construction and by measurement.** The commutator accumulates in
  the same interleaved order as `_commutator_batched_core`; the scalar factors are
  formed in the same association (`(0.5*h)`, `((c*h)*h)`); the sum order matches. dOm
  measured exactly 0.0 on every configuration, and |dP| exactly 0.0 end-to-end on the
  PREM chord at 544/1088/2176 slabs, d=2-5. No mpmath re-quantification needed.
- The constant-sample decision is the same exact-equality test, evaluated with early
  exit -- the same *result* array_equal produces, cheaper. A castle-wall profile still
  takes the h*A1 path bit-identically.
- No signature changes (private helper gains an internal fast path), no module globals.
  Same warnings (none live here), same exceptions.

## Optional extensions, to judge separately

- **Variant B -- build At inside the kernel.** A second entry point taking
  (HE_c, V, mA, widths) would fold the scan's At broadcast build (0.018-0.062 us/slab/E,
  report 09 has the table) into the same pass and halve the kernel's memory traffic.
  It binds the separable structure H = H_E + V*h_matt into magnus.py, which is the
  multiply-out the author flagged: one more kernel, scan-only, still general in d.
  Worth ~1.15x further at d=3. Defensible, but it is the first step down the road where
  every engine wants its own fused kernel; decide deliberately.
- **Order-6 sibling.** The `order <= 6` branch (3 commutators, a1/a2/a3 builds) fuses
  the same way for a similar share. More code (~40 lines); order 4 is the default and
  the paper's workhorse, so do 4 first and let 6 wait for a demonstrated need.

## Tests that bear on it

- `tests/test_magnus_expansion.py` (34 tests): the per-slab constancy tests (castle-wall
  profiles) pin the h*A1 fast path; the GL-order convergence tests pin the arithmetic.
  Both should pass unchanged -- bit-identity makes this the rare change where "no test
  moves" is the expected outcome, and any test that moves falsifies the claim.
- `tests/test_engines.py` engine-equivalence checks cross-check the scan output.

## Risks, and what would falsify it

- Low. The one real risk is drift between the kernel and the numpy fallback if the
  order-4 expression is ever edited in one place and not the other -- same standing risk
  `_commutator_batched` already carries, mitigated the same way (a test comparing the
  two paths at exact equality on random stacks).
- Falsifier: any nonzero |dP| in the `fused_e2e.py` A/B, or a nonzero dOm on random
  inputs, means the accumulation order was not reproduced and the change must not land
  as "bit-identical" (it would then need the mpmath quantification instead).
