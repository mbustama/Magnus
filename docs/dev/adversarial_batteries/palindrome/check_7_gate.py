"""What can the gate actually test, and what does each choice cost?

The mirror exploits  A(mirror slab samples) == reverse(A(forward slab samples)).
A gate must establish that BEFORE the work is skipped.  Three candidates:

  W  widths only          -- cheap, computable from edges alone
  S  the sampled A        -- exact, but requires evaluating A everywhere
  D  caller declares it   -- free, but moves responsibility to the producer
"""
import time
import numpy as np
import common
import magnus.magnus as mg
import magnus.globaldefs as gd

L, A, H, vcc = common.chord_setup()
proj = np.zeros((3, 3)); proj[0][0] = 1.0
hv = H(0.0) - vcc(0.0)*np.diag([1.0, 0, 0])
V0 = 4.0e-13


def A_sym(l):
    l = np.asarray(l, dtype=float)
    v = V0*(1.0 + 0.5*np.cos(2.0*np.pi*l/L))
    return -1j*(hv + v[..., None, None]*proj)


def A_mono(l):                      # monotonic: solar-like, NOT symmetric
    l = np.asarray(l, dtype=float)
    v = V0*np.exp(-3.0*l/L)
    return -1j*(hv + v[..., None, None]*proj)


print("=" * 78)
print("1.  Does 'widths are palindromic' imply the mirror is valid?")
print("=" * 78)
n = 64
edges, w = common.symmetrised_edges(L, n)
print(f"  uniform grid, widths symmetrised: array_equal(w, w[::-1]) = "
      f"{np.array_equal(w, w[::-1])}")
for name, Af in (("symmetric profile", A_sym), ("MONOTONIC profile", A_mono)):
    At, wid = common.sample_A(Af, edges, None, 'gl', 4)
    Us = mg.magnus_expansion_multislab(Af, edges, order=4,
                                       integration_method='gl',
                                       validate_input=False)
    # what the mirror would produce: evaluate first half, mirror the rest
    m = n//2
    An = At[:m]
    Om_f = mg._magnus_gl(An, wid[:m], 4)
    Om_b = mg._magnus_gl(An[:, ::-1], wid[:m], 4)
    Om = np.empty((n, 3, 3), dtype=complex)
    Om[:m] = Om_f; Om[n-m:] = Om_b[::-1]
    Up = mg._expm_stack(Om)
    print(f"    {name:20s}  widths pass gate, max|U_mirror - U_shipped| = "
          f"{np.max(np.abs(Up - Us)):.3e}")

print()
print("=" * 78)
print("2.  Is the sampled A bitwise palindromic on a symmetrised grid?")
print("=" * 78)
for order in (2, 4, 6):
    At, wid = common.sample_A(A_sym, edges, None, 'gl', order)
    mirror = At[::-1, ::-1]
    print(f"  gl order {order}: array_equal(At, At[::-1,::-1]) = "
          f"{str(np.array_equal(At, mirror)):5s}   "
          f"max|diff| rel = {np.max(np.abs(At-mirror))/np.max(np.abs(At)):.2e}")


# ---------------------------------------------------------------------------
# A third section timed three gate routes here.  It was REMOVED rather than
# fixed, because it did not measure what it claimed: the comparison routes
# called _expm_stack() with the default warn_wide=False while the shipped
# magnus_expansion_multislab() passes warn_wide=True, which runs an SVD per
# slab for the convergence check.  The "speed-up" was mostly the absence of
# that check.  A replacement must call the shipped entry point on both sides,
# or pass warn_wide identically.
#
# The two sections above need no timing harness: both are exact-arithmetic
# facts about what a gate can test, and they are what constrains the design.
# See PLAN_PALINDROMIC_PROFILES.md section 3d(ii).
# ---------------------------------------------------------------------------
