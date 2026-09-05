# Brief: fuse `_expm_stack`'s anti-Hermiticity framing (report 03)

Design and evidence: `report-03-fused-antihermiticity-check.md` in this directory.
Prototype: `prototypes/fusion_proto2.py` (`_antiherm_scale_dev2`), harness
`prototypes/fusion_speed2.py`.

## The job

`_expm_stack` (`src/magnus/magnus.py`, now line 2016 -- the report's 1608 is stale,
the file has grown) opens by asking "is this stack anti-Hermitian" in five
full-stack temporaries and about seven memory passes, before any exponential is
computed. Replace that with one compiled pass returning `(scale, dev)`, keeping
both branches and building `K = 1j*Om` only after the branch is taken.

## The one decision that is not the report's to make

The report recommends comparing **squared magnitudes** and taking one `sqrt` at
the end, because computing `abs()` per element in numba is *slower* than NumPy's
vectorized passes -- its first prototype measured 0.4-0.7x, and the squared form
measured 4.6-4.9x. Take that recommendation, but be exact about what it costs:

`sqrt(re^2 + im^2)` and `hypot(re, im)` **differ in the last bit for some inputs**.
The argmax element is identical by monotonicity, so the *decision* is stable, but
the returned `scale` can differ from NumPy's by <= 1 ulp -- and `scale` feeds
`_warn_slab_norm`, so it is not purely internal.

**Every other kernel in this series has been bit-identical by construction, not by
measurement, and that is the standard the repository owner has applied all along.**
This one cannot be, in the recommended form. So:

1. Implement the squared form.
2. **Also measure the `hypot`-per-element variant**, which the report says is
   exactly value-identical at roughly half the speedup.
3. Report both numbers side by side, and state plainly what each buys and costs.

Do not decide between them. Implement the squared form, measure both, and let the
owner choose. If the hypot variant turns out to cost less than the report's "half"
estimate, say so -- that would change the decision.

## Re-measure the shares; the report's are stale

Report 03 says the framing is 12% of a d=3 call and 19% at d=2, and quotes stage
costs of 0.0292 / 0.0506 / 0.0895 / 0.1470 us/slab/E at d = 2/3/4/5. **Those
predate the Jacobi backend**, which landed today and made the d=4 and d=5
exponentials 1.8-1.9x and 1.5-1.6x cheaper end to end. The framing's *share* at
four and five flavours is therefore larger now than the report says, and its
absolute per-slab cost may have moved too. Measure the shares yourself at all four
flavour counts rather than quoting the report.

## Constraints

- **Surgical edits.** Do not rewrite files or reflow untouched code; touch only
  what must change.
- **Update or add every docstring and in-line comment the edit affects**, and give
  the new helper a real docstring in the style of `_gl4_omega_core`'s -- what it
  fuses, what it returns, why the squared form is used, and the 1-ulp caveat.
- **No signature changes**, no change to what any routine takes or returns.
- Preserve both branches exactly: `scale == 0.0` (identity) and
  `dev <= 1e-12*scale`. `sqrt(0) == 0` exactly, so the identity branch is safe.
- Non-complex128, non-contiguous and numba-less input must fall through to the
  existing expression untouched.
- Do **not** touch `CHANGELOG.md` or `pyproject.toml`; versioning is the owner's.
- Do **not** commit.

## Verification

- End-to-end `|dP|` against the unfused path, both profiles, d = 2..5: report the
  worst value, and say whether it is exactly 0.0 or merely small.
- The four tests the report names in `tests/test_expm_backend.py` pin exactly this
  behaviour and should pass unchanged:
  `test_non_anti_hermitian_input_still_reaches_scipy`,
  `test_slab_norm_warning_still_fires_on_the_kernel_path`,
  `test_real_valued_hermitian_input_is_accepted`,
  `test_non_contiguous_input_is_handled`.
- Run `python3 -m pytest tests/test_expm_backend.py tests/test_engines.py tests/test_magnus_expansion.py -q`.
  Do **not** run the full suite; it takes 25 minutes and the owner will run it.
- Speed: marginal us/slab, fused vs unfused, arms interleaved, control workload
  reported. The machine is otherwise idle; keep it that way.

## Report

What you changed; the measured gain at each flavour count; the framing's share
*as re-measured*, not as the report has it; the two variants' numbers side by
side; the worst `|dP|`; and explicitly anything you could not make exact.
