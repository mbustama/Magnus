# Brief: the Jacobi 4nu/5nu backend -- measure the gate, then implement

Two phases, and **phase 2 happens only if phase 1 passes**. The gate is a single
ratio. If it fails, stop and report; do not implement, and do not widen a
tolerance to make it pass.

Design and evidence: `report-01-jacobi-expm-4nu-5nu.md` in this directory, with
the final prototype kernel at `prototypes/jacobi_proto5.py`
(`_jacobi_expm_warm_mgs2`).

---

# Phase 1 -- the gate

## What is already settled, so you do not re-derive it

Checked by reading the assertions rather than running them.

The suite's tight tolerances are **not** flat across norms. `_herm` is called
with `scale` in {1.0, 1e-8, 1e-13}, so the fixed-tolerance assertions run near
norm 1, where any backward-stable method sits at ~eps. The large-norm case has
its own test, `test_kernel_is_no_worse_than_eigh_at_any_norm`, which sweeps
`scale` over 1e-150 to 1e4 and sets

    tol = max(1e-15, 1e-14*scale)

At norm 1e4 that is 1e-10, against the 1.7e-11 report 01 measures there --
inside, with ~6x of margin. So the absolute tolerances are **not** the binding
constraint, and the 1.7e-11 figure was being read against the wrong yardstick.

## The binding constraint, which is the gate

The second assertion in that same test:

    assert err_k <= 10.0*max(err_e, 1e-17)

The kernel must stay within **10x of eigh at every norm**. Report 01 claims
parity of scaling but never reports this ratio. It is the one number that decides
adoption.

## What to measure

`err_k / err_e` at **d = 4 and d = 5**, across the same norm sweep the test uses
(`scale` in 1e-150, 1e-8, 1.0, 10.0, 100.0, 1e3, 1e4), on stacks built the way
that test builds them:

    K = _herm((32, d, d), rng); K *= scale/np.max(np.abs(K))

with `err = max|U - scipy_expm(K)|` for each backend. Build the stacks directly
with `_herm`; do **not** regenerate the `om_d4/om_d5.npy` caches, which the
prototypes read but which were not copied into the repository.

Also cover the two spectrum structures the suite cares about at d <= 3, since
they are where a method cliffs rather than degrades: clustered eigenvalues
(splittings down to 1e-14 and exactly 0.0) and exact degeneracies, each crossed
with the norm sweep. That crossing is not decoration -- a one-axis-at-a-time
sweep over this exact code missed a 7440x accuracy hole once, and the comment
block above `test_two_condition_grid` in the test file records it.

## The control that makes the measurement mean anything

**Score eigh on every point you score Jacobi on.** The assertion is *relative*:
a Jacobi error without its eigh partner is uninterpretable. Report both, and
report margins rather than verdicts -- "1.2x of eigh" and "40x of eigh" call for
opposite decisions, and "fails" hides which one you saw.

## Phase 1 output

A table of scale x dimension x spectrum-structure, giving `err_k`, `err_e` and
the ratio. Then one plain sentence: does `err_k <= 10*err_e` hold everywhere, and
if not, where and by how much.

---

# Phase 2 -- implement, ONLY if the ratio holds everywhere

Follow report 01's design. Integration points, from that report:

- `expmkernels._jacobi_expm_core(K, out, lam) -> sev`, same contract as
  `_ch2_core` / `_ch3_core`; one kernel with loop bounds from the shape, not
  specialized per dimension.
- `expmkernels.supports_dim`: return True for 4 and 5, **and rewrite its
  docstring** -- it currently says "A 4x4 or 5x5 Hermitian eigenproblem has no
  practical closed form, so 4nu and 5nu stay on eigh", and it calls itself "the
  one place that decision is made". The premise is fine; the conclusion is what
  changes, and the new docstring should say why: eigh was the cost, not the
  missing closed form.
- `expm_herm_stack`: route d = 4, 5 to the new kernel. d = 2, 3 untouched.
- `_expm_stack`: no edit needed -- the numba branch already fires when
  `supports_dim` says yes, and the `sev > SEV_TOL` hook already provides the eigh
  fallback.
- The `ValueError` for unsupported d moves from ">= 4" to ">= 6".

Load-bearing details, not optional:

- **Modified Gram-Schmidt on the warm-start basis at every step.** Without it
  non-unitarity compounds along the chain: 2.3e-11 after 13k matrices against
  3.9e-14 with it.
- **`fastmath` stays off.**
  `test_compiled_kernel_matches_the_same_source_uncompiled` requires the compiled
  and uncompiled forms to agree *bitwise*, and the stability argument depends on
  cancellations the compiler must not reassociate.
- Return `sev = 0.0`; if any matrix hits the 30-sweep cap, return a `sev` above
  `SEV_TOL` so the existing hook recomputes with eigh.
- Eigenvalues ascending, to honour the documented `expm_herm_stack` contract.

## Tests

- `test_dim_four_and_five_delegate_to_eigh` and
  `test_public_kernel_entry_point_refuses_unsupported_dimensions` pin the current
  routing and **must be updated** -- they test internal dispatch, not public
  contract.
- `test_kernel_is_no_worse_than_eigh_at_any_norm` hardcodes a `(32, 3, 3)` stack;
  extend it to d = 4 and 5. This is a real edit, not a parametrization.
- Extend the agreement, exact-degeneracy, ascending-eigenvalue and two-condition
  tests to d = 4, 5 as well.
- Run `python3 -m pytest tests/test_expm_backend.py tests/test_engines.py -q`.
  Do **not** run the full suite; it takes 25 minutes and the repository owner
  will run it.

## Standing constraints

- **Surgical edits.** Do not rewrite files or reflow untouched code.
- **Update every docstring and in-line comment your edit affects.**
- **No signature changes**, and no change to what any routine takes or returns.
- Do **not** touch `CHANGELOG.md` or `pyproject.toml`; versioning is the
  repository owner's call.
- Do **not** commit. Leave the work in the tree.

## Report

Phase 1's table and verdict. Then, if you implemented: what you changed, the
accuracy you measured against eigh and scipy at d = 4, 5, which tests you edited
and why, and explicitly anything you could not make pass. If you widened any
tolerance, say so in the first line of the report -- that is a finding, not a fix.
