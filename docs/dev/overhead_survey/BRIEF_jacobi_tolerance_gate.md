# Brief: does the Jacobi backend pass the tolerances the suite actually asserts?

**This is a measurement, not an implementation.** Do not change `src/`. The
question is whether the batched Jacobi eigensolver of report 01 could be adopted
at d = 4 and 5 without weakening any assertion the test suite already makes. The
answer decides whether the option is taken up at all, so it must be measured
rather than argued.

## Why this is the gate

Report 01 proposes extending `expmkernels.supports_dim` from {2, 3} to
{2, 3, 4, 5}. Two and three flavours keep the closed-form kernels they have, so
the paper's three-flavour claims -- including the 2.9e-13 reach -- cannot move.
The four- and five-flavour claims sit far above the reported shift.

What is *not* settled is the test suite. `tests/test_expm_backend.py` holds 32
tight-tolerance assertions, several at 1e-14 against `scipy.linalg.expm` and at
1e-13 against the eigh backend, and today they run only at d <= 3. Report 01
recommends extending them to d = 4, 5 -- and separately measures Jacobi's
operator-level `dU` against eigh at **1.7e-11 at norm 1e4**, which is above both
tolerances. That may be entirely benign, because the same report states eigh has
the same `eps*||K||` scaling; or it may not. Nobody has run the extended
assertions.

## What to do

1. Read `tests/test_expm_backend.py` and extract, for each tight-tolerance
   assertion, the **actual configuration it uses**: matrix norms, spectrum
   structure (clustered, exactly degenerate, generic), batch shapes, and the
   tolerance asserted. Do not guess these; the norms are the whole question.
2. Build the d = 4 and d = 5 analogue of each such assertion.
3. Score **three** backends on that battery, not one:
   - the Jacobi prototype (`prototypes/jacobi_proto5.py`, kernel
     `_jacobi_expm_warm_mgs2`),
   - `np.linalg.eigh`, which is what d = 4, 5 use today,
   - `scipy.linalg.expm` as the reference where the assertion uses it.
4. Report, per assertion and per dimension: **pass or fail, and the margin** --
   the measured value against the asserted tolerance, for Jacobi *and* for eigh.

## The control that makes this worth running

**eigh must be scored on the identical battery.** If an extended assertion fails
for eigh too, the assertion does not transfer to higher dimensions and the
finding is about the test, not about Jacobi. Reporting Jacobi's failures without
eigh's would produce exactly the wrong conclusion, and this codebase has been
burned by a one-armed comparison before. Every number you report for Jacobi must
have its eigh counterpart beside it.

Likewise, report margins rather than verdicts. "Fails 1e-14, measured 1.2e-14"
and "fails 1e-14, measured 3e-11" call for opposite decisions.

## Practical notes

- The prototypes read `om_d4.npy` / `om_d5.npy`, which were **not** copied into
  the repository (8.7 MB). `prototypes/grab_om.py` regenerates them.
- Do not run the full suite; it takes 25 minutes.
- Measure on an idle machine. Check first that no other job is running.

## Report

A table of assertion x dimension x backend with margins, then a plain statement
of one of: (a) Jacobi passes everything eigh passes, (b) it fails assertions eigh
also fails, so the assertions do not transfer, or (c) it fails assertions eigh
passes -- and if so, exactly which, at what margin, and whether the failure is a
conditioning corner or a systematic gap. Do not recommend adoption or rejection;
report the measurement and let the repository's author decide.
