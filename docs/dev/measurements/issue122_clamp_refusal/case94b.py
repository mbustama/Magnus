"""The #94 test's case (3nu, order 4, 'gl', cap 115, rtol = atol = 1e-6): the true error of what it returns."""
import warnings, numpy as np
import magnus.oscprob as op, magnus.globaldefs as gd
rho = lambda l: 3e3*np.exp(-np.asarray(l, dtype=float)/(100.0*gd.UNIT_KM))
osc = gd.load_nufit_params('NuFIT 6.1')
E = np.linspace(0.1, 0.4, 40)*gd.UNIT_GEV
def run(**k):
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter('always')
        P = np.asarray(op.osc_prob_matter_std_potential(3, rho, E, 250.0*gd.UNIT_KM, osc, L0=0.0,
                       density_matter_is_in_g_per_cm3=True, strategy='magnus', **k), float)
    return P, sorted({x.category.__name__ for x in w})
P, w = run(rtol=1e-6, atol=1e-6, max_n_slabs=115)
R4, w4 = run(rtol=1e-12, atol=1e-14, magnus_exp_order=4)
R6, w6 = run(rtol=1e-12, atol=1e-14, magnus_exp_order=6)
print('refs: order 4 vs 6 differ by %.1e' % np.abs(R4 - R6).max())
tol = 1e-6 + 1e-6*np.abs(R6)
def score(X):
    err = np.abs(X - R6)
    bad = np.any(err > tol, axis=tuple(range(1, err.ndim))).sum()
    return 'max |P - ref| = %.2e, max err/tol = %.2f, energies outside: %d of %d' % (err.max(), (err/tol).max(), bad, len(E))
print('capped run (warnings %s): %s' % (w, score(P)))
fixed = {}
for k in (72, 108, 115, 162):                    # one grid each: a loose tolerance accepts the first repeat
    fixed[k], _ = run(rtol=1.0, atol=1.0, n_slabs=k, max_n_slabs=k)
    print('fixed %3d slabs: %s' % (k, score(fixed[k])))
print('returned == fixed 115: %s' % np.array_equal(P, fixed[115]))
d = np.abs(fixed[115] - fixed[108]); print('108 vs 115 differ by %.2e (max d/tol %.2f)' % (d.max(), (d/tol).max()))
d = np.abs(fixed[108] - fixed[72]);  print(' 72 vs 108 differ by %.2e (max d/tol %.2f)' % (d.max(), (d/tol).max()))
