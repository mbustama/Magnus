"""Is the sign rule Omega_k -> (-1)^(k+1) Omega_k algebraically false for k>=3,
or is the 3.7e-08 residual a convergent quadrature artifact?

Discriminator: refine n_tpts_per_slab.  An algebraic failure sits still.
A quadrature artifact goes to zero.

The mirror slab's samples are built by REVERSING the forward slab's samples,
so the two slabs are exact mirrors by construction -- no profile symmetry or
grid symmetry is involved, and nothing but the identity itself is under test.
"""
import numpy as np
import common
import magnus.magnus as mg

L, A, H, vcc = common.chord_setup()

N_SLABS = 64
edges, w = common.symmetrised_edges(L, N_SLABS)


def residuals(order, method, n_tpts):
    At, widths = common.sample_A(A, edges, n_tpts, method, order)
    Bt = widths[:, None, None, None]*At
    T_f = mg._magnus_terms_quadrature(Bt, order, method)          # (order, n, d, d)
    T_m = mg._magnus_terms_quadrature(Bt[:, ::-1], order, method)  # true mirror
    signs = np.array([(-1.0)**(k + 1) for k in range(1, order + 1)])
    T_sr = signs[:, None, None, None]*T_f
    per_term = [np.max(np.abs(T_m[k] - T_sr[k])) for k in range(order)]
    scale = [np.max(np.abs(T_m[k])) for k in range(order)]
    om_err = np.max(np.abs(T_m.sum(axis=0) - T_sr.sum(axis=0)))
    om_scale = np.max(np.abs(T_m.sum(axis=0)))
    return per_term, scale, om_err, om_scale


print("=" * 78)
print("A.  Sign-rule residual on Omega, vs n_tpts_per_slab  (64 slabs, PREM chord)")
print("=" * 78)
for method in ('trapezoid', 'simpson'):
    for order in (2, 4, 6):
        print(f"\n  {method}, order {order}")
        print(f"    {'n_tpts':>7} | {'max|Om_mirror - signrule|':>26} | {'relative':>10}")
        print("    " + "-" * 52)
        prev = None
        for n_tpts in (11, 21, 41, 81, 161, 321, 641):
            _, _, e, s = residuals(order, method, n_tpts)
            ratio = "" if prev is None else f"  ({prev/e:5.1f}x down)" if e > 0 else ""
            print(f"    {n_tpts:>7} | {e:>26.4e} | {e/s:>10.2e}{ratio}")
            prev = e

print()
print("=" * 78)
print("B.  Per-term breakdown (trapezoid, order 6) -- which k actually breaks?")
print("=" * 78)
for n_tpts in (41, 161, 641):
    pt, sc, _, _ = residuals(6, 'trapezoid', n_tpts)
    print(f"\n  n_tpts = {n_tpts}")
    print(f"    {'k':>3} | {'max|Om_k(mirror) - (-1)^(k+1) Om_k|':>36} | {'relative':>10}")
    print("    " + "-" * 56)
    for k in range(6):
        print(f"    {k+1:>3} | {pt[k]:>36.4e} | {pt[k]/sc[k]:>10.2e}")

print()
print("=" * 78)
print("C.  Is the residual below the scheme's OWN truncation error?")
print("    (compare against a 4001-point trapezoid reference for the same order)")
print("=" * 78)
for order in (4, 6):
    print(f"\n  trapezoid, order {order}")
    print(f"    {'n_tpts':>7} | {'sign-rule residual':>19} | {'quadrature error':>19} | ratio")
    print("    " + "-" * 66)
    At_ref, widths = common.sample_A(A, edges, 4001, 'trapezoid', order)
    Bt_ref = widths[:, None, None, None]*At_ref
    Om_ref = mg._magnus_terms_quadrature(Bt_ref, order, 'trapezoid').sum(axis=0)
    for n_tpts in (41, 161, 641):
        _, _, e, _ = residuals(order, 'trapezoid', n_tpts)
        At, widths = common.sample_A(A, edges, n_tpts, 'trapezoid', order)
        Bt = widths[:, None, None, None]*At
        Om = mg._magnus_terms_quadrature(Bt, order, 'trapezoid').sum(axis=0)
        q = np.max(np.abs(Om - Om_ref))
        print(f"    {n_tpts:>7} | {e:>19.4e} | {q:>19.4e} | {e/q:8.3f}")
