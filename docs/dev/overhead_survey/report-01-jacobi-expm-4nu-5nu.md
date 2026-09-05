# Batched Jacobi eigensolver for the 4x4 and 5x5 exponential

**Measured end-to-end gain: 1.77x at 3+1, 1.38x at 3+2** (marginal us/slab/E on the PREM
chord, order-4 GL, fixed slabs, control ratio 1.00). Stacked with the bit-identical fusions
of reports 02-03, the combined stack reaches 2.09x at 3+1 -- which takes Mag(nu)s from 1.53x
*slower* than NuOscProbExact at 3+1 to **0.73x, i.e., faster** (1.247 vs 1.70 us/slab/E).

## The cost, and how it was measured

`expmkernels.supports_dim` stops at d=3 ("no practical closed form for a 4x4 or larger
Hermitian eigenproblem, so 4nu and 5nu stay on eigh"). That sentence is true of *closed
forms* and yet the conclusion no longer holds: on the current code, `np.linalg.eigh` is
**66% of the whole d=4 call and 63% at d=5** (cProfile, PREM chord, 1088 slabs, 12
energies; confirmed by stage timing: the eigh-plus-reconstruction path costs 2.02 of 2.55
us/slab/E at d=4 and 2.75 of 3.51 at d=5). The per-matrix cost is ~2.3 us for a 4x4 --
LAPACK `zheevd`'s fixed overhead, not arithmetic. This is the answer to "higher flavour
counts, where the arithmetic may now dominate": it is not the arithmetic, it is eigh.

Stage timings from `scratchpad/survey/stages.py`; end-to-end A/B from
`scratchpad/survey/jacobi_e2e.py` / `all_e2e.py` (alternated arms, first timing discarded,
marginal fitted over 1104->2192 actual slabs, interleaved control 1.00).

## The change

A numba kernel `_jacobi_expm_core(K, out, lam) -> sev` in `expmkernels.py`, same contract
as `_ch2_core`/`_ch3_core`, dispatched from `expm_herm_stack` for d in (4, 5) (one kernel,
loop bounds from the shape -- it is not specialized per dimension). Per matrix:

1. **Cyclic complex-Hermitian Jacobi**: for each pivot (p, q), the standard unitary
   rotation with phase `e^{i phi} = A[p,q]/|A[p,q]|` and angle from
   `tau = (A[q,q]-A[p,p])/(2|A[p,q]|)`, `t = sign(tau)/(|tau|+sqrt(1+tau^2))`;
   rows/columns p, q of A and columns p, q of the accumulated V updated in O(d).
   Sweeps until the off-diagonal Frobenius norm is <= 1e-16 of the matrix norm
   (5.0 sweeps cold at d=4 on the real stacks, 6.0 at d=5), hard cap 30.
2. **Warm start**: initialize V from the previous matrix's converged eigenvectors and
   rotate `A <- V0^H K V0` (two small matmuls). Consecutive matrices are consecutive
   slabs of the same energy, so A arrives nearly diagonal: sweeps drop to 3.05 (d=4) and
   4.75 (d=5), and per-matrix cost from 1409 to 919 ns (d=4) and 3122 to 1717 ns (d=5),
   vs eigh's 2288/2939 ns. At an energy boundary in the flattened stack the warm start
   is merely less effective for one matrix; correctness never depends on it.
3. **Modified Gram-Schmidt on V0 at every step.** Without it, warm-start non-unitarity
   compounds linearly along the chain: measured |dU| vs eigh of 2.3e-11 after 13k
   matrices, vs 3.9e-14 with MGS. This is load-bearing, not optional.
4. **Reconstruction from the MGS'd V** (`U = sum_j e^{-i lam_j} v_j v_j^H`, permutation-
   invariant, exp hoisted out of the inner loop): per-slab unitarity then matches eigh
   (composed-product unitarity 1.3e-12 vs eigh's 1.0e-12 at 4368 slabs; before this
   refinement it was 3-4x worse). Eigenvalues insertion-sorted ascending to honor the
   documented `expm_herm_stack` contract.

Working prototypes, already validated: `scratchpad/survey/jacobi_proto5.py` (the final
kernel, `_jacobi_expm_warm_mgs2`), with `jacobi_proto2.py` the cold-start variant.

Integration points:
- `expmkernels.supports_dim`: return True for 4 and 5, and rewrite its docstring -- it is
  "the one place that decision is made", which is exactly why this is a small change.
- `expm_herm_stack`: route d=4,5 to the new kernel; d=2,3 untouched.
- `_expm_stack` needs no edit at all: the numba branch already runs whenever
  `supports_dim` says yes, the `sev > SEV_TOL` hook already provides the eigh fallback,
  and the eigh path remains reachable via `expm_backend='eigh'` exactly as documented.
- The `ValueError` for unsupported d moves from ">= 4" to ">= 6".

## Correctness

- **No conditioning cliff, so no gate.** The closed-form d<=3 kernel needed `SEV_TOL`
  because arccos loses sqrt(eps)*||K|| digits on clustered spectra at large norm. Jacobi
  is backward stable with no such corner: measured against eigh AND scipy.expm at norm
  1e4 with clustered (1e-9 split) and exactly degenerate spectra, dU is 1.7e-11 across
  the board -- eps*||K||, the same scaling eigh itself has, and identical between the
  clustered and unclustered cells (large-norm battery in `scratchpad/survey`). Return
  sev = 0.0; if any matrix hits the 30-sweep cap without converging, return a sev above
  `SEV_TOL` instead so the existing hook recomputes the stack with eigh -- the same
  safety architecture d<=3 already has, exercised only in a corner never observed
  (Jacobi on Hermitian matrices converges unconditionally; the cap is a backstop).
- **Scored against the independent mpmath references** (prem_chord_reference.json and
  mp_reference_profile.json, numu->numu, both profiles, d=4 and 5, 1088 and 4352 slabs):
  reference error unchanged where discretization dominates (5.7890e-07 -> 5.7890e-07);
  at the finest grid's floor, expo d=4 moves 7.40e-12 -> 7.71e-12 and expo d=5
  2.44e-11 -> 2.56e-11. **Worst full-matrix probability shift 2.8e-12.** That is larger
  than the 6.7e-14 of the commutator kernel, and the right yardstick is the shipped d<=3
  backend itself: its documented guarantee vs expm is **5e-12 absolute**, so this shift
  is in the family the package has already accepted for a backend replacement, not a new
  class. Unitarity of returned probabilities matches the eigh arm at every setting.
- Warnings: the kernel returns eigenvalues, so `_warn_slab_norm(max|lam|)` sees the same
  quantity computed by a different (equally accurate) algorithm -- same conditions up to
  the rounding class above. Exceptions, dispatch, convergence decisions: untouched.
- No new module globals, no shared mutable state, no signature changes anywhere. The
  warm start and MGS live entirely inside one kernel invocation's locals.

## Tests that bear on it

- `tests/test_expm_backend.py` (36 tests): `test_dim_four_and_five_delegate_to_eigh` and
  `test_public_kernel_entry_point_refuses_unsupported_dimensions` **pin the current
  routing and must be updated** -- they test internal dispatch, not public contract.
  `test_backends_agree_on_random_hermitian_stacks`, `test_exact_degeneracies`,
  `test_eigenvalues_are_returned_ascending`, `test_kernel_is_no_worse_than_eigh_at_any_norm`
  and the calibration battery should be *extended* to d=4,5 (the large-norm battery above
  is the template; a `docs/dev/calibrate_sev_tol.py`-style sweep would be the thorough
  version).
- `test_probabilities_agree_between_backends` exercises exactly this A/B at the
  probability level.

## Risks, and what would falsify it

- The 2.8e-12 shift is the decision point: if the author holds the line at ~1e-13, use
  the cold-start variant (no cross-matrix state at all, dU vs eigh 2.1e-14/5.8e-14) --
  but it only wins 1.61x on the exponential stage at d=4 and ties eigh at d=5, so most
  of the d=5 gain is the warm start. Falsifier: run `survey/jacobi_acc.py` with the cold
  kernel; if its shift still exceeds the bar, the whole idea is off.
- Speedup depends on the slowly-varying-slab structure (warm start). On a profile with
  genuinely uncorrelated consecutive slabs the sweeps revert toward 5-6 and the d=5 gain
  shrinks to ~1.0x; d=4 keeps ~1.3-1.6x. The PREM chord and the exponential profile both
  have the structure; a fuzz/adversarial profile may not, but correctness is unaffected.
- numba compile time: one more kernel (~0.5-1 s first call), same class as the existing
  ones, amortized by `cache=True`.
- Headroom left inside the kernel itself: the rotation updates traverse full rows and
  columns where the Hermitian-aware update needs about half, worth an estimated further
  20-30% -- worth doing during implementation, not worth re-benchmarking the design for.
