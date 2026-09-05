# Cache `_passthrough_kwarg_names` (59 us of `inspect.signature` per call)

**Fixed 59 us per public-entry call, removable with a one-line lru_cache -- 11% of the
~541 us fixed overhead of a single-point call.** Irrelevant to the marginal us/slab
metric; real for the single-point exposure the project tracks.

## The cost, and how it was measured

`_passthrough_kwarg_names` (src/magnus/oscprob.py:3751) recomputes
`inspect.signature` over three functions (osc_prob, osc_prob_energy_baseline,
magnus_expansion_multislab) on **every** call to any wrapper that validates kwargs.
Measured directly (timed, control interleaved): 59.1 us/call, essentially all of
`_check_passthrough_kwargs`'s 58.9 us. A single-energy PREM call at n_slabs=16 costs
551 us total, of which ~541 us is fixed overhead (fit vs n_slabs=544), so this one
helper is ~11% of the whole fixed cost.

## The change

`@functools.lru_cache(maxsize=1)` on `_passthrough_kwarg_names` (it takes no arguments
and returns a frozenset), or an explicit lazy module-level memo. The derived-not-listed
design principle is preserved: the set is still read off the signatures, just once per
process instead of once per call.

## Correctness, and the one honest cost

- The signatures are module-level `def`s; they cannot change at runtime except by
  monkeypatching, so the cached value is the same value.
- **This is a module global** (the cache), which constraint 4 says to flag plainly:
  a test that monkeypatches `osc_prob`'s signature (none currently does -- verified by
  the absence of signature patching in tests/) and then relies on the error message
  listing the patched keyword would see stale names. `lru_cache` at least exposes
  `.cache_clear()` for such a test to call.
- Same errors with the same messages (the difflib suggestion path is untouched), same
  everything else.

## What would falsify it

- If profiling shows the remaining ~480 us of fixed overhead has larger single items
  (it does -- dispatch probes and Hamiltonian setup; see HANDOVER_SINGLE_POINT_EXPOSURE
  for that program), this is still worth taking because it is one line, but it is not
  the lever that changes the single-point story on its own.
