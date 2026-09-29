"""2nu Sun, 10 MeV, 0.9 R_sun: the general ladder at 1e-4 (auto's route at 1e-3) under three floors."""
import time, warnings, numpy as np
import magnus.globaldefs as gd, magnus.oscprob as op
sth, Dm2 = np.sqrt(0.308), 7.5e-5
energy, L = 10.0*gd.UNIT_MEV, 0.9*gd.SUN_RADIUS*gd.UNIT_KM
Pref = np.asarray(op.osc_prob_2nu_sun(energy, L, 0.0, sth, Dm2, strategy='hybrid', rtol=1e-10, atol=1e-10, validate_input=False))
print('hybrid 1e-10 P_ee=%.10f' % Pref[0][0])
for label, floor in (('no floor', 1), ('traceless floor', 3500), ('trace floor', 18334)):
    with warnings.catch_warnings(record=True) as w, op._engine_probe(disabled=('ip_exp', 'hybrid')):
        warnings.simplefilter('always')
        t = time.perf_counter()
        P = np.asarray(op.osc_prob_2nu_sun(energy, L, 0.0, sth, Dm2, strategy='magnus', rtol=1e-4, atol=1e-4,
                                           min_n_slabs=floor, validate_input=False))
        dt = time.perf_counter() - t
    print('%-16s t=%.3fs err=%.1e %s' % (label, dt, abs(P[0][0] - Pref[0][0]), sorted({x.category.__name__ for x in w})))
