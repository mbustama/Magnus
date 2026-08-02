# Bug: a batched solar call can allocate hundreds of gigabytes

**Found:** 2026-08-02, while measuring the energy axis for `DECISION_OSCPROB_CUMULATIVE.md`.
**Status: FIXED** on branch `fix-ip-exp-memory`; see §"The fix" below.
Pre-existing — reproduced unchanged at `155e01e`, before the `n_slabs` floor landed.
**Severity:** took the machine down, not just the process. It OOM-killed the desktop
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

## The fix

**A. Tile the working set over (energy, slab).** The per-level temporaries are now built one
tile at a time against a fixed entry budget (`BATCH_WORKING_ENTRIES`, via
`_tile_for_working_set`), so peak memory is a property of the library rather than of the
call. Measured, at 65 536 slabs:

| energies | 4 | 16 | 64 | 256 |
|---|---|---|---|---|
| peak | 79.3 MiB | 78.7 MiB | 78.6 MiB | 78.6 MiB |

Flat across a 64-fold increase in the energy count. The call that used to die — eight solar
energies over a full solar radius — now completes.

The tiling is **exact, not approximate**. Within a tile the arithmetic is elementwise, so
slicing moves no value; and the slab product is folded in the same descending order with the
accumulator on the left, so the parenthesis nesting — the only thing that could shift a
floating-point result — is unchanged. Verified across a 2000-fold range of budgets, with the
peak tracking the budget (2.3 MiB at a 16 384-entry budget, 627 MiB at 33 554 432) and the
result bit-identical at every one of them:

| budget (entries) | 16 384 | 131 072 | 1 048 576 | 4 194 304 | 33 554 432 |
|---|---|---|---|---|---|
| peak | 2.3 MiB | 4.5 MiB | 21.5 MiB | 80.1 MiB | 626.8 MiB |
| identical to untiled | yes | yes | yes | yes | yes |

Time is flat across that range (12–16 s, no trend), because the cost is dominated by the
per-slab Python fold rather than by array size. There is no speed/memory trade to make here.

Two traps met while verifying this, both worth carrying forward:

- **A module constant consumed as a default argument cannot be monkeypatched.**
  `_tile_for_working_set` first took `max_entries=BATCH_WORKING_ENTRIES`, which binds at
  import; the test that varied the budget was therefore comparing two identical runs and
  passed for free. It now defaults to `None` and reads the constant at call time. The budget
  sweep above is what exposed it — the peak was suspiciously identical at every budget.
- **Check that a passing test can fail.** With the correct block order the tiled product is
  exactly equal; with the blocks walked the wrong way it is off by O(1) — so `array_equal` is
  a genuine discriminator here, not an assertion that happens to hold.

One subtlety worth recording, because the first attempt got it wrong: the budget has to
bound the engine's **whole** working set, not one array of it. Eight temporaries of the same
shape are live at the peak (the argument, the `exp` temporary, the slab integral, `Omega`,
the matrix exponential's eigenvectors and workspace, the slab operators, and the
accumulator's operand), so a budget applied per array overshoots eightfold — 625 MiB where
79 MiB was intended. `_tile_for_working_set` now takes `live_arrays` explicitly.

**B. Refuse up front when certification is provably impossible.** Certification requires
`max|Omega_t|` below the trust threshold, and that maximum is bounded below by the diagonal
entries, which have a closed form here (`Delta_jj = 0`, so the slab integral collapses to
`l_scale*(1 - exp(-w/l_scale))`). If that bound still exceeds the threshold at the slab
ceiling, no reachable slab count can certify, and the twenty-one doublings that follow are
guaranteed waste. Two evaluations of `VCC_func` and no allocation detect it.

It is a bound, never an estimate: it can only say *impossible*, so it cannot abandon a case
that would have worked. Where it bites, at 10 MeV over a solar radius:

| requested tolerance | 1e-2 | 1e-3 | **1e-4** | **1e-6** |
|---|---|---|---|---|
| slabs needed | 455 558 | 1 440 611 | **4 555 621** | **45 556 251** |
| vs cap 2 000 000 | climbs | climbs | **refuses at once** | **refuses at once** |

Note honestly that it does **not** fire at the default tolerance below one solar radius:
`n_min = 1 440 611` sits just under the ceiling, and — because `max|Mt_jj|` is the matter
mixing angle, which stays near 1 — that number is the same at every energy from 0.1 to
100 MeV. Fix B covers tighter tolerances and longer baselines; Fix A is what makes the
default case safe.

**C. Refuse a result that cannot fit, before allocating it.** Tiling bounds the engines;
nothing shrinks the answer. `osc_prob_energy_baseline` now checks the requested
`N*d*d` against the operating system's free-memory figure and raises a `MemoryError` naming
the size, rather than letting an overcommitting kernel turn it into a machine-wide kill. The
check costs one multiply below a 64 MiB floor, and never blocks on a platform where free
memory cannot be read.

## Still open

- **The ladder still climbs to ~1.44M slabs at the default tolerance**, taking ~10 s per
  energy against the per-point path's 72 ms. That is now merely slow rather than dangerous,
  but it means the fast path is ~130x *slower* than the general one on the case it exists
  for. Whether `_osc_prob_ip_exp_dispatch`'s gate should be narrowed — or the method given a
  second interaction picture that absorbs the diagonal matter phase exactly, which would
  destroy its closed-form slab integral and is therefore research — is a separate question,
  and should be settled by measuring `ip_exp` against `hybrid` and per-point across the
  solar range now that doing so is safe.
- The suite still has no memory coverage outside these tests. `_osc_prob_scan_separable` was
  measured safe (100 MiB at 10 energies to 330 MiB at 4000) because it already tiled; it now
  shares the same helper and constant, so the knob is in one place.
