# Handover — allow a batched Hamiltonian in NuOscProbExact's tolerance route

**For:** a session working in `~/Research/NuOscProb/NuOscProbExact` (at `1a6a704` when this
was written).
**From:** the Magνs session, which needs the change but must not make it.

**Do the audit in Part 1 and report back. Change nothing until that report is reviewed.**

---

## 1. What is wanted, in one sentence

`probabilities_{2,3,4}nu_profile(...)` currently accepts a `hamiltonian_of` callable that
returns exactly one Hamiltonian per position, shape `(n, d, d)`. The request is to let it
also accept a leading batch axis, `(..., n, d, d)`, so that a whole array of energies can
be refined to a tolerance in one call instead of one call per energy.

## 2. Why, and who is asking

The Magνs paper's Figure 11 compares Magνs against NuOscProbExact on a smooth exponential
profile: deviation from a converged reference against wall-clock time per probability, at
2, 3, 4 and 5 flavours. Today the two codes are dialled by *different knobs* — a requested
tolerance for Magνs, a slab count for NuOscProbExact — which makes the curves awkward to
read against each other. Putting NuOscProbExact on the same `rtol` sweep would let every
point on both curves answer the same question: "you asked for 1e-8; here is what it cost
and what you got."

The obstacle is timing fairness. That benchmark times **twelve energies in one batched
call**, and its own source notes that timing them one at a time "would flatter Mag(nu)s
roughly fivefold". `probabilities_Nnu_slabs` takes the batch and is what the current curve
uses. `probabilities_Nnu_profile` does not, so a naive switch to the tolerance route would
put a fivefold handicap on NuOscProbExact's time axis for reasons that have nothing to do
with tolerances. That would be a misleading figure, and it is not going in the paper.

## 3. The precise site

`src/slabs.py`, inside `_probabilities_profile` (defined at **line 1016**), in the nested
`evaluate`:

```python
        h = np.asarray(hamiltonian_of(midpoints), dtype=complex)
        if h.shape != (n, n_flavors, n_flavors):
            raise ValueError(
                '%s: hamiltonian_of returned shape %s for %d positions; it '
                'must return one %dx%d Hamiltonian per position, of shape '
                '(%d, %d, %d)' % (caller, (h.shape,), n, n_flavors,
                                  n_flavors, n, n_flavors, n_flavors))
        return np.asarray(routine(h, np.diff(edges)), dtype=float)
```

The proposed relaxation is to make the comparison a suffix match:

```python
        if h.shape[-3:] != (n, n_flavors, n_flavors):
```

with the message reworded from "must be of shape (n, d, d)" to "must end in (n, d, d)".

Its three public callers are `probabilities_2nu_profile` (line 1177),
`probabilities_3nu_profile` (1284) and `probabilities_4nu_profile` (1364).

## 4. What the Magνs session already verified — confirm, do not re-derive

Each of these was checked by reading the source in this repo. Please confirm them
independently; if any is wrong, the whole proposal changes.

1. **`_n_for_tolerance` (line 384) is shape-agnostic.** Its acceptance tests are
   `np.all(4.0*gap/3.0 <= a + r*np.abs(p_coarse))` and
   `np.all(gap/3.0 <= a + r*np.abs(p_fine))`; its diagnostic is `float(np.max(gap/3.0))`.
   Nothing assumes a scalar or a single energy. No change should be needed there.

2. **`probabilities_Nnu_slabs` already accepts `(n_energies, n, d, d)`.** The Magνs
   benchmark calls it that way today and indexes the result as `P[..., col]`, so batched
   input returns shape `(n_energies, d**2)` and single input returns `(d**2,)`. The last
   axis is the flattened probability index in both cases.

3. **The pattern already exists in this package.** `src/earth.py` line 1881 passes array
   energies through `slabs._n_for_tolerance` and then re-scalarises with
   `if np.ndim(energy) == 0:`. So batched tolerance refinement is established behaviour
   here; `_probabilities_profile` is the outlier, not the innovation.

4. **The guard is validation, not an algorithmic constraint.** It was written to catch a
   caller returning the wrong shape. It rejects a batch axis as collateral damage.

## 5. PART 1 — THE AUDIT. Report before changing anything.

Please establish, and report, the following. Where the answer is "nothing to change", say
so explicitly; a short list of confirmed non-issues is as useful as a list of issues.

### 5.1 Would it break anything?

- **Callers inside this repo.** Every call site of `probabilities_{2,3,4}nu_profile` and of
  `_probabilities_profile`. Does any of them rely on the exact-shape rejection — for
  instance a test asserting that a wrong shape raises?
- **The test suite.** Search for tests that assert the `ValueError` from this guard, by
  message text or by shape. Those will need updating, and their intent should be preserved:
  a test that `(d, d, n)` is rejected is still valid; a test that *any* non-`(n,d,d)` shape
  is rejected is not.
- **Return-shape assumptions downstream.** With a batch axis, `evaluate(n)` returns
  `(n_energies, d**2)` and `_probabilities_profile` returns that shape too. Does anything
  in the package or its tests assume the profile routines return a flat tuple of `d**2`
  floats? Note `earth.py` line ~1886 has precedent for handling both.
- **`return_n_slabs=True`.** Confirm the second return value stays a scalar `n` under
  batching, which is the intended semantics (see 5.3).

### 5.2 What weakens, and is that acceptable?

A suffix match still rejects `(n, d)`, `(d, d, n)` and `(n, d, d+1)`. It newly *accepts*
`(k, n, d, d)` from a caller who meant `(n, d, d)` and got their axes wrong — a mistake
that is currently caught. Please judge whether that trade is acceptable, and whether the
reworded error message and docstring make the batch axis discoverable enough that the
mistake is unlikely. If you would rather not weaken it, an explicit `batched=False`
keyword is the alternative; say so and the Magνs side will adapt.

### 5.3 One behavioural point that must be documented

Refinement becomes **all-entries-at-once**: because the tests are `np.all`, the slab count
is set by the hardest entry in the stack, and every entry is returned at that count. A
user batching twelve energies therefore gets a different (never coarser) slab count than
calling twelve times. That is the correct semantics for the Magνs comparison — Magνs also
refines an energy array until the whole array satisfies the request — but it is a real
behavioural statement and should be in the docstring rather than discovered.

### 5.4 Where documentation would need to change

Please locate and report, without editing:

- The `hamiltonian_of` parameter description in all three public `..._profile` docstrings
  (they say "returns the Hamiltonian at each, as an array of shape `(len(positions), d, d)`").
- Their `Returns` sections, which would gain the batched shape.
- Any `Examples` blocks in those docstrings — note whether they are executed in the docs
  build, since an executed example is a test.
- `docs/source/*.rst` — `methodology.rst` and `functions.rst` are the likely ones; check
  all of them for statements about the profile routines' input or output shape.
- `notebooks/*.ipynb` — there are at least eight. Check whether any uses
  `probabilities_*_profile`, and whether any prose describes the one-Hamiltonian-per-
  position contract. **Do not edit `.ipynb` files directly if this repo generates them from
  a script**, as the Magνs repo does; check for a generator first.
- `README.md`, `CHANGELOG.md`, `RELEASE_NOTES.md`.

### 5.5 Versioning

`probabilities_3nu_profile` carries `.. versionadded:: 1.12.0`. This is a backward-
compatible API extension, so it wants a minor bump and a changelog entry. Please report
the current version and where it is defined — a `__version__` grep came back empty from
outside the repo, so it may live in `pyproject.toml` or `globaldefs.py`.

## 6. PART 2 — the change, once the audit is reviewed

Only after the report is agreed:

1. The one-line guard relaxation and the reworded message.
2. Docstring updates in the three public routines: the batch axis in `hamiltonian_of`, the
   batched return shape, and the all-entries-at-once refinement note from 5.3.
3. Whatever documentation and notebook changes the audit identified.
4. Changelog entry and version bump.

## 7. Acceptance tests

The change is right when:

- **Equivalence.** For an energy array `E`, a batched `hamiltonian_of` returning
  `(len(E), n, d, d)` gives probabilities equal to looping `probabilities_Nnu_profile` over
  the energies **at a fixed `n_slabs` with `rtol=None, atol=None`** — that comparison is
  exact and should agree to round-off. Under a tolerance the two need not match, because
  batched refinement uses one slab count for all (5.3); assert instead that the batched
  answer is at least as accurate as each per-energy answer.
- **Unchanged behaviour.** An unbatched `(n, d, d)` caller returns bit-identical results to
  the current release. Worth a regression test against captured values —
  `tests/bit_capture.py` and `tests/bit_compare.py` suggest that machinery exists.
- **The guard still guards.** `(n, d)`, `(d, d, n)` and `(n+1, d, d)` still raise.
- **Timing.** One batched call over twelve energies is meaningfully faster than twelve
  calls. That is the entire point; if it is not, report back before proceeding, because the
  premise is then wrong.

## 8. What the Magνs side will do with it

Add a `NuOscProbExact, rtol` series to `notebooks/external_profile_benchmarks.json`,
swept over the same tolerances as Magνs (1e-3, 1e-4, 1e-6, 1e-8, 1e-10), **keeping** the
existing `n_slabs` series in the file, and plot the tolerance one in Figure 11.

One thing to be aware of: `N_SLABS_MAX = 1024` (`src/slabs.py` line 324) is the default
`n_max`, and it is far too small for this use — the existing `n_slabs` curve needs 32,768
slabs to reach 2.6e-11, and `_n_for_tolerance` *raises* on exhaustion rather than returning
its best effort. Magνs will pass a large explicit `n_max`, so no change is required here,
but if you think the default is low for tolerance-driven use generally, that is worth
raising separately.
