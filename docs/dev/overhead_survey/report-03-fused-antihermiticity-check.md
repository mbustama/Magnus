# Fuse `_expm_stack`'s anti-Hermiticity framing into one compiled pass

**~4.7x on the framing stage, bit-identical end-to-end in measurement** (|dP| exactly
0.0); the framing is 12% of the whole call at d=3 and 19% at d=2, and at d=2 it costs
more than the exponential it guards (0.029 vs 0.022 us/slab/E).

## The cost, and how it was measured

`_expm_stack` (src/magnus/magnus.py:1608) opens with

```
K = 1j*Om
Kh = np.conj(np.swapaxes(K, -1, -2))
scale = np.max(np.abs(K))
... np.max(np.abs(K - Kh)) <= 1.e-12*scale
```

-- five full-stack temporaries and ~7 memory passes to answer "is this stack
anti-Hermitian", on every chunk, before any exponential is computed. Stage cost
(us/slab/E): 0.0292 (d=2), 0.0506 (d=3), 0.0895 (d=4), 0.1470 (d=5). HANDOVER_OVERHEAD
S9.3 already named this framing as "the next bottleneck, not the exponential" for
single-matrix calls; it is a per-slab cost in the scan too.

## The change

A numba helper in `magnus.py`:

```
_antiherm_scale_dev(Om) -> (scale, dev)
    # one pass over the (nB, d, d) stack, no temporaries:
    #   scale = sqrt(max over elements of re^2+im^2)          == max|K|
    #   dev   = sqrt(max over elements of |Om + conj(Om^T)|^2) == max|K - K^H|
    # (|i z| = |z| exactly, so working on Om equals working on K)
```

`_expm_stack` then computes `scale, dev = _antiherm_scale_dev2(flat_Om)`, keeps the
`scale == 0.0` and `dev <= 1e-12*scale` branches exactly as they are, and only builds
`K = 1j*Om` (one pass) after the branch is taken. Non-complex128, non-contiguous, or
numba-less input falls through to the existing expression untouched.

Important detail found by prototyping: computing `abs()` per element in numba (hypot) is
*slower* than numpy's vectorized passes -- the first prototype measured 0.4-0.7x. The
win comes from comparing **squared magnitudes** and taking one sqrt at the end
(max is order-preserved): 4.6-4.9x measured at d=2-5
(`scratchpad/survey/fusion_proto2.py` `_antiherm_scale_dev2`, harness
`fusion_speed2.py`).

## Correctness

- The argmax element is identical (monotonicity); the returned *value* can differ from
  numpy's hypot-based one by <=1 ulp on that single element, because
  sqrt(re^2+im^2) != hypot(re, im) in the last bit for some inputs. Measured difference
  on the benchmark workloads: exactly 0.0, and |dP| exactly 0.0 end-to-end. The branch
  decision could flip only for an input whose deviation-to-scale ratio sits within one
  ulp of 1e-12 exactly -- four decades from where any real input lands (anti-Hermitian
  inputs measure ~1e-16, non-Hermitian ones >> 1e-12). If even that is unwanted, the
  kernel can compute hypot per element for exact value-identity at roughly half the
  speedup; the squared form is the recommendation.
- `scale` and the eigenvalues still feed `_warn_slab_norm` identically; the
  `scale == 0.0` identity branch is unchanged (sqrt(0) == 0 exactly).
- No signature changes, no module globals. The scipy fallback for genuinely
  non-anti-Hermitian input is reached under exactly the same condition.

## Tests that bear on it

- `tests/test_expm_backend.py`: `test_non_anti_hermitian_input_still_reaches_scipy`,
  `test_slab_norm_warning_still_fires_on_the_kernel_path`,
  `test_real_valued_hermitian_input_is_accepted`, `test_non_contiguous_input_is_handled`
  -- all pin exactly the behaviors this touches, and all should pass unchanged.

## Risks, and what would falsify it

- Lowest of the three kernel proposals. Falsifier: a nonzero ds/dd on any workload in
  `fusion_speed2.py`, or any |dP| != 0 in `fused_e2e.py`.
- Note `_expm_stack` is also called by the IP engine per block (report 04) and by every
  other engine, so this gain is shared; the IP engine's framing share (78 ms of a 757 ms
  pass) shrinks with the same patch.
