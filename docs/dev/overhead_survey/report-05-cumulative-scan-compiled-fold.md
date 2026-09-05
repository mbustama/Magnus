# Compile the cumulative baseline scan's running product

**The Python running-product loop costs ~6 us/slab against ~1.25 us/slab for building
the operators it consumes -- about 2-4x on the scan at large slab counts.**

## The cost, and how it was measured

`_osc_prob_cumulative_scan` (src/magnus/oscprob.py:5577; the loop at ~5650) walks the
profile once and snapshots the running product at every requested baseline:

```
for k in range(stop - start):
    running = U[k] @ running
    while (pos < len(order)) and (out_idx[order[pos]] == start + k + 1):
        P[j] = np.transpose(running.real**2 + running.imag**2); ...
```

Measured on a 3nu vacuum scan of 400 baselines (cProfile, `scratchpad/survey/
prof_cum.py`): the function's own body is 5 of 7 ms while the two
`compute_evolution_operator_multiple_slabs` calls that build every operator cost ~1 ms
-- the fold is a single (d,d)@(d,d) numpy matmul (~2-3 us dispatch) plus a transpose per
snapshot, per slab, in Python. Per slab: ~6 us fold vs ~1.25 us construction. The
absolute numbers grow linearly with `n_acc`: a fine accuracy grid (1e4-1e5 slabs, which
is what tight rtol produces on long baselines) pays 60-600 ms per (energy) call in pure
dispatch.

## The change

A numba kernel taking the block's operators, the in/out accumulator, and the snapshot
indices, so the fold and the snapshot writes happen compiled:

```
_running_product_snapshots(U, running, snap_k, P_out, p0) -> None
    # U: (nb, d, d); running: (d, d) in/out; snap_k: sorted local snapshot
    # indices; P_out: (n_out, d, d) float, filled from p0
    # for k in 0..nb-1:
    #     running <- U[k] @ running          (same association: new on the LEFT)
    #     while next snapshot == k: P_out[j] = transpose(|running|^2)
```

The Python caller keeps building `snap_k` from `out_idx` exactly as it does now (the
argsort bookkeeping is cheap); only the per-slab loop moves. Note the association here
is `U[k] @ running` -- accumulator on the RIGHT -- the mirror of the IP engine's fold;
the kernel must preserve each engine's own order, which is a reason to implement them as
two small kernels (or one kernel with a side flag set at the call sites, not exposed).

## Correctness

- Association preserved; rounding may move at the 1e-14 class exactly as for
  `_ordered_product_batched` -- quantify against the references before landing.
- Snapshot placement is pure integer bookkeeping; a test comparing snapshot positions
  and count against the current implementation on randomized L_out grids (including
  baselines coinciding with edges and with L0) pins it. The `out_idx == 0` identity
  snapshots stay in Python -- they precede the loop.
- No signature changes, no module globals; the memory-bounding block structure is
  untouched (the kernel processes one already-materialized block).

## Tests that bear on it

- test_oscprob.py's cumulative-scan tests (grep `cumulative` in tests/test_oscprob.py,
  test_engines.py) -- per-point vs cumulative agreement is the behavioural pin.
- DECISION_OSCPROB_CUMULATIVE.md records why the scan exists; nothing there depends on
  fold speed.

## Risks, and what would falsify it

- Small, but so is the payoff unless baseline scans at fine grids are a real workload;
  the notebooks that plot P(L) curves are (each figure's L-scan pays this loop).
  If profiling a real figure regeneration shows the scan under a few percent of the
  figure's cost, drop this one -- it is the least of the three fold proposals.
- Falsifier: probability shift beyond the accepted class, or any snapshot landing on a
  different slab index than today for any L_out (that would be a bug, not rounding).
