import time, warnings, numpy as np
import magnus.globaldefs as gd, magnus.oscprob as op, magnus.magnus as mg
sth, Dm2 = np.sqrt(0.308), 7.5e-5
energy, L = 10.0*gd.UNIT_MEV, 0.9*gd.SUN_RADIUS*gd.UNIT_KM
for strat in ('auto', 'hybrid'):
    info = {}
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter('always')
        t = time.perf_counter()
        P = np.asarray(op.osc_prob_2nu_sun(energy, L, 0.0, sth, Dm2, strategy=strat, validate_input=False, strategy_info=info))
        dt = time.perf_counter() - t
    tr = [x for x in info.get('trace', []) if 'estimated_phase' in x]
    print(strat, 'P_ee=%.8f t=%.2fs engine=%s' % (P[0][0], dt, info.get('engine')), tr, [(x.category.__name__, str(x.message)[:160]) for x in w])
