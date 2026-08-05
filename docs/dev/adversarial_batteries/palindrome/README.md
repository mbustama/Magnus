# Palindromic-profile battery

Reproduces every number in `../../PLAN_PALINDROMIC_PROFILES.md` §3d, which corrects two claims
made in §0/§3b/§3c of that plan and in §§2.3/2.7/2.8 of
`../../HANDOVER_PALINDROME_AND_SESSION_LOG.md`.

Run from this directory (`common.py` inserts `src/` on the path itself):

```bash
for f in check_*.py; do echo "== $f"; python "$f"; done
```

| script | what it establishes |
|---|---|
| `check_1_grid_and_Uj.py` | the Earth chord's grid is palindromic to 4.30e-15 **relative** (rounding, not physics); and `U_j != U_{n-1-j}` at Magnus order >= 2, which is why NuOscProbExact's composer does not port |
| `check_2_signrule.py` | the sign-rule residual falls **4.0x per doubling** of `n_tpts_per_slab` — O(h^2), the trapezoid's own rate — and sits 4-5 orders below the quadrature error the same grid already commits |
| `check_3_generic_and_transpose.py` | the same O(h^2) convergence on a **generic random complex** `A(t)` at orders 2-8, so it is not an artefact of the PREM Hamiltonian's structure |
| `check_4_transpose.py` | the transpose identity `U = F^T F` needs `H^T = H`: exact at `dCP = 0`, fails at 7.4e-01 at `dCP = 3.70`. Uses **matched slab widths** — an earlier version compared a 64-slab full evolution against a 64-slab half and mistook the discretisation mismatch (1.3e-06) for the identity failing |
| `check_5_gl_prototype.py` | the handover's `gl_mirror` prototype is exact to 1e-15 at even slab counts and returns **uninitialised memory for the middle slab** at odd ones |
| `check_6_external.py` | confirms the sign rule against `solve_ivp` rather than against magnus internals: sign-flipped terms converge to the *reversed* propagator through order 6 exactly as unflipped terms converge to the forward one. Were the rule false at k >= 3 this would stall at order 3 |
| `check_7_gate.py` | a widths-only gate passes a **monotonic** profile and the mirror is then wrong by 3.34e-01; and the sampled `A` is never bitwise palindromic (~5e-16), so an `array_equal` gate on it never fires |

## Reading these

`check_6` is the one to trust most: it is the only check that does not use the code under test as
its own reference. `check_2` and `check_3` are internal-consistency checks and would agree with
each other even if `_magnus_terms_quadrature` were wrong in a way that affected forward and
reversed samples alike.

`check_7` once carried a third section timing three gate routes. It was removed rather than
fixed: the comparison routes called `_expm_stack(warn_wide=False)` while the shipped entry point
uses `warn_wide=True`, which runs an SVD per slab, so the "speed-up" was largely the absence of
that check. **No timing claim in §3d rests on this battery** — the surviving results are
exact-arithmetic facts, which is why they need no timing harness. Any future timing work should
follow the traps in the handover §4 (read minima not medians; sub-millisecond timings are noise).
