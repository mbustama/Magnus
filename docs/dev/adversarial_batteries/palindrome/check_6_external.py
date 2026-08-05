"""Independent confirmation against scipy's solve_ivp, not magnus internals.

Claim under test:  Omega_k[A_reversed] = (-1)^(k+1) Omega_k[A]  in the continuum.

Build Omega_k from A on a fine grid, sign-flip the even ones, exponentiate, and
compare against the propagator that solve_ivp produces for the REVERSED problem.
If the identity is true, the mismatch is pure Magnus truncation: it must fall as
the Magnus order rises.  If the identity is false at k>=3, raising the order
cannot help and the mismatch will stall.
"""
import numpy as np
from scipy.integrate import solve_ivp
import magnus.magnus as mg

rng = np.random.default_rng(7)
d = 3
c0, c1, c2 = (rng.normal(size=(3, d, d)) + 1j*rng.normal(size=(3, d, d)))
a, b = 0.0, 0.9


def A_fwd(t):
    t = np.asarray(t, dtype=float)[..., None, None]
    return 0.6*(c0 + c1*np.sin(3.0*t) + c2*np.cos(2.0*t)*t)


def A_rev(s):
    return A_fwd(a + b - np.asarray(s, dtype=float))


def exact_U(Afun):
    def rhs(t, y):
        return (Afun(t) @ y.reshape(d, d)).ravel()
    sol = solve_ivp(rhs, (a, b), np.eye(d, dtype=complex).ravel(),
                    rtol=1e-13, atol=1e-14, method='DOP853', dense_output=False)
    return sol.y[:, -1].reshape(d, d)


U_rev_exact = exact_U(A_rev)
U_fwd_exact = exact_U(A_fwd)
print("  reference propagators from solve_ivp (rtol=1e-13)")
print(f"    max|U_fwd| = {np.max(np.abs(U_fwd_exact)):.4f}, "
      f"max|U_rev| = {np.max(np.abs(U_rev_exact)):.4f}")
print(f"    U_fwd != U_rev by {np.max(np.abs(U_fwd_exact - U_rev_exact)):.3e} "
      "(so this is a real, non-trivial test)")
print()

N_TPTS = 20001                       # fine enough that quadrature error is negligible
s = np.linspace(0.0, 1.0, N_TPTS)
w = b - a
Bt = w*A_fwd(a + w*s)[None, ...]     # one slab, shape (1, m, d, d)

print("=" * 76)
print("  Sign-flipped forward terms vs the REVERSED problem's exact propagator")
print("=" * 76)
print(f"  {'order':>6} | {'max|expm(signrule) - U_rev_exact|':>34} | falls with order?")
print("  " + "-" * 72)
prev = None
for order in (1, 2, 3, 4, 5, 6):
    T = mg._magnus_terms_quadrature(Bt, order, 'trapezoid')     # (order,1,d,d)
    sg = np.array([(-1.0)**(k + 1) for k in range(1, order + 1)])
    Om = (sg[:, None, None, None]*T).sum(axis=0)[0]
    e = np.max(np.abs(mg._expm_stack(Om[None, ...])[0] - U_rev_exact))
    tag = "" if prev is None else f"{prev/e:8.1f}x better"
    print(f"  {order:>6} | {e:>34.4e} | {tag}")
    prev = e

print()
print("  Control: the SAME terms without sign flips, vs the FORWARD propagator")
print("  " + "-" * 72)
prev = None
for order in (1, 2, 3, 4, 5, 6):
    T = mg._magnus_terms_quadrature(Bt, order, 'trapezoid')
    Om = T.sum(axis=0)[0]
    e = np.max(np.abs(mg._expm_stack(Om[None, ...])[0] - U_fwd_exact))
    tag = "" if prev is None else f"{prev/e:8.1f}x better"
    print(f"  {order:>6} | {e:>34.4e} | {tag}")
    prev = e
