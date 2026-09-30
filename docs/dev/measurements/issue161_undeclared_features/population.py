"""#161 design population: how often each candidate trigger fires, on bad and on good calls."""
import warnings, time, json, sys, numpy as np
warnings.simplefilter('ignore')
import magnus.oscprob as o, magnus.globaldefs as gd, magnus.hamiltonians as hm, magnus.avgprob as ap, magnus.matter as matter
from scipy.integrate import solve_ivp
KM = gd.UNIT_KM; eps = np.finfo(float).eps
Vf = lambda rho: np.sqrt(2)*gd.GF*rho*gd.UNIT_G_PER_CM3/((gd.MASS_PROTON + gd.MASS_NEUTRON)/2)*0.5
e00 = np.diag([1., 0, 0])
p = gd.OSC_PARAMS_PREDEFINED['OSC_PARAMS_DEFAULT']; osc = {k: p[k] for k in ('s12','s23','s13','dCP','D21','D31')}
hv3 = np.asarray(hm.hamiltonian_3nu_vacuum_energy_independent(**osc))
# Record every accepted per-point ladder level
LOG = []
_orig = o.osc_prob
def _logged(*a, **k):
    ci = k.get('convergence_info')
    if ci is None:
        ci = {}; k['convergence_info'] = ci
    r = _orig(*a, **k)
    P = r[0] if isinstance(r, tuple) else r
    LOG.append(dict(ci, maxP=float(np.max(np.abs(P)))))
    return r
o.osc_prob = _logged

def dop(Hl, L, splits):
    psi = np.eye(3, dtype=complex); pts = [0.0] + sorted(s for s in splits if 0 < s < L) + [L]
    for a, b in zip(pts[:-1], pts[1:]):
        s = solve_ivp(lambda l, y: (-1j*Hl(l) @ y.reshape(3, 3)).ravel(), (a, b), psi.ravel(), method='DOP853', rtol=1e-12, atol=1e-14)
        psi = s.y[:, -1].reshape(3, 3)
    R = np.abs(psi)**2
    return R.T   # P[i, f]

rows = []
def run(tag, kind, fn, ref):
    LOG.clear()
    t = time.perf_counter(); P = np.asarray(fn()); dt = time.perf_counter() - t
    acc = [c for c in LOG if c.get('tolerance_achieved') is not None]
    c = acc[-1] if acc else {}
    gap = c.get('last_gap'); mp = c.get('maxP', 1.0)
    rows.append(dict(tag=tag, kind=kind, err=float(np.abs(P - ref).max()) if ref is not None else None, ms=1e3*dt,
                     n=c.get('n_slabs'), nprev=c.get('n_slabs_previous'), gap=gap,
                     eps_gap=bool(gap is not None and gap <= 1e3*eps*mp), n_ladder_calls=len(LOG)))

# P1: castle wall, osc_prob_energy_baseline, 40 baselines
E = 3e9
Hc = lambda l: hv3/E + Vf(3.3 if (l < 3000*KM or l > 9000*KM) else 11.5)*e00
for L in np.linspace(1000, 12000, 40)*KM:
    ref = dop(Hc, L, [3000*KM, 9000*KM])
    run('castle L=%.0f' % (L/KM), 'bad?', lambda: o.osc_prob_energy_baseline(lambda e, l: Hc(l), E, L), ref)
# P2: spikes via osc_prob_matter_std_potential (default auto)
for E in (0.5e9, 1e9, 3e9):
    for pos in (500, 1000, 2500):
        for w in (2, 5, 20, 50):
            rho = lambda x, pos=pos, w=w: 3 + 50*np.exp(-((np.asarray(x, float) - pos*KM)/(w*KM))**2)
            vcc = matter.vcc_func_from_rho_func(rho, density_matter_is_in_g_per_cm3=True)
            Hs = lambda l, E=E, vcc=vcc: hv3/E + float(vcc(l))*e00
            ref = dop(Hs, 4000*KM, [(pos - 4*w)*KM, (pos + 4*w)*KM])
            run('spike E=%.1f pos=%d w=%d' % (E/1e9, pos, w), 'bad?',
                lambda: o.osc_prob_matter_std_potential(3, rho, E, 4000*KM, osc, density_matter_is_in_g_per_cm3=True), ref)
# P3: negatives
for E in (0.5e9, 3e9, 10e9):
    for cz in (-1.0, -0.8, -0.5, -0.2):
        run('earth E=%.1f cz=%.1f' % (E/1e9, cz), 'good', lambda: o.osc_prob_3nu_earth(E, costhz=cz, L=__import__('magnus.earth',fromlist=['x']).distance_traveled_inside_earth(cz)*KM), None)
    run('exp E=%.1f' % (E/1e9), 'good', lambda: o.osc_prob_3nu_matter_exp_density(E, 3000*KM, 0.0, 10.0, 1000*KM), None)
    run('linear E=%.1f' % (E/1e9), 'good', lambda: o.osc_prob_matter_std_potential(3, lambda x: 2 + 8*np.asarray(x)/(5000*KM), E, 5000*KM, osc, density_matter_is_in_g_per_cm3=True), None)
    run('sine E=%.1f' % (E/1e9), 'good', lambda: o.osc_prob_matter_std_potential(3, lambda x: 5 + 3*np.sin(np.asarray(x)/(800*KM)), E, 6000*KM, osc, density_matter_is_in_g_per_cm3=True), None)
    run('const-callable E=%.1f' % (E/1e9), 'good', lambda: o.osc_prob_energy_baseline(lambda e, l: hv3/e + Vf(4.0)*e00, E, 4000*KM), None)
for E in (1e6, 5e6, 10e6):
    run('sun E=%.0fMeV' % (E/1e6), 'good', lambda: o.osc_prob_3nu_sun(E, gd.SUN_RADIUS*KM, 0.0), None)
json.dump(rows, open(sys.argv[1], 'w'), indent=0)
bad = [r for r in rows if r['kind'] == 'bad?']
wrong = [r for r in bad if r['err'] > 1e-3]
print('bad-population calls: %d, of which err>1e-3: %d' % (len(bad), len(wrong)))
print('  err>1e-3 with eps-gap acceptance: %d' % sum(r['eps_gap'] for r in wrong))
print('  err<=1e-3 with eps-gap acceptance: %d' % sum(r['eps_gap'] for r in bad if r['err'] <= 1e-3))
good = [r for r in rows if r['kind'] == 'good']
print('good calls: %d, eps-gap acceptances: %d, ladder-less: %d' % (len(good), sum(r['eps_gap'] for r in good), sum(r['n'] is None for r in good)))
for r in rows:
    if r['kind'] == 'good' and r['eps_gap']: print('  GOOD eps-gap:', r)
for r in wrong:
    if not r['eps_gap']: print('  MISSED:', r['tag'], 'err %.1e n %s<-%s gap %s' % (r['err'], r['n'], r['nprev'], r['gap']))
