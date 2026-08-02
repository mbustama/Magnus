# Bug: a batched solar call can allocate hundreds of gigabytes

**Found:** 2026-08-02, while measuring the energy axis for `DECISION_OSCPROB_CUMULATIVE.md`.
**Status:** not fixed. Pre-existing — reproduced unchanged at `155e01e`, before the
`n_slabs` floor landed.
**Severity:** takes the machine down, not just the process. It OOM-killed the desktop
application it was running under three times before it was tracked down.

---

## What happens

```python
import numpy as np
import magnus.globaldefs as gd, magnus.oscprob as oscprob

En = np.linspace(5.0, 6.0, 8)*gd.UNIT_MEV          # eight solar-neutrino energies
L = gd.SUN_RADIUS*gd.CONV_KM_TO_INV_EV
oscprob.osc_prob_2nu_sun(En, L, 0.0, gd.S12_NO_BF_NUFIT_6_0, gd.D21_NO_BF_NUFIT_6_0)
```

```
MemoryError: Unable to allocate 977. MiB for an array with shape (8, 2000000, 2, 2)
             and data type complex128
```

Measured peak RSS, same call, varying only the number of energies:

| energies | peak RSS | wall time |
|---|---|---|
| 1 | 1.56 GB | 9.6 s |
| 2 | 2.84 GB | 15.8 s |
| 4 | 5.34 GB | 32.5 s |
| 8 | `MemoryError` under a 6 GB cap | — |

Linear in the number of energies, at roughly **1.3 GB and 8 s per energy**. Nothing about
this call is pathological: `osc_prob_2nu_sun` is public and documented, an array of energies
is its advertised batched form, and 0.1–100 MeV over a solar radius is the documented use
case. Notebook 03 cell 104 scans **1000** energies over exactly that range — through raw
`osc_prob` calls, which is the only reason the notebooks have never hit this.

## Why

`_osc_prob_ip_exp_core` (`oscprob.py`, the `while True` at ~2836) allocates, *per refinement
level*, several simultaneous arrays of shape `(nE, n_slabs, d, d)` complex128 — `arg`, `I`,
`Omega_t`, `U_slab`, plus whatever `magnus._expm_stack` needs — while

```python
n_slabs_cap = IP_EXP_N_SLABS_CAP      # 2_000_000
growth = 2.0
```

doubles `n_slabs` toward two million. One such array at the cap, for a single 2-level
energy, is `2e6 * 4 * 16 = 128 MB`; there are several live at once, and every one of them is
multiplied by `nE`.

The reasoning recorded in the source is explicitly about *time*, and it is correct as far as
it goes:

> each slab here costs one 2x2 eigendecomposition, so even the full ceiling completes in a
> couple of seconds

What it does not consider is memory, or that batching over energies multiplies the working
set by `nE`. The slab budget is also deliberately decoupled from the caller's
`max_n_slabs`, so a caller who tries to bound the cost cannot.

A second observation, worth separating from the memory bug: in this configuration the ladder
runs all the way to the cap without converging, taking ~10 s per energy. The per-point path
answers the same query in **72 ms**. So on this configuration the "fast path" is also about
130× *slower* than the general one.

## Suggested directions (not yet evaluated)

- Bound the working set rather than the slab count: process the energy batch in chunks sized
  so that `nE_chunk * n_slabs * d^2` stays under a budget, or stream slabs the way
  `DECISION_OSCPROB_CUMULATIVE.md` §6.1 describes for the cumulative scan.
- Reconsider `IP_EXP_N_SLABS_CAP = 2_000_000` in light of `nE`: a cap that is defensible for
  one energy is not defensible for a thousand.
- Investigate why the ladder does not converge here at all, since that is what drives it to
  the cap; the trust threshold is tied to `sqrt((atol+rtol)/2)` and may simply be
  unreachable for this profile at these energies.
- Whatever the fix, a regression test should pin peak memory for a batched solar call, not
  only its accuracy — this failure is invisible to every existing test because the suite
  never batches enough energies to make the allocation large.
