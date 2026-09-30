r"""Exhaustive audit of the #160 checklist: one case per sub-case.

R: must raise a Magnus error ('Error in magnus' or argparse 'error:') whose text contains `name`.
W: must emit a warning whose class name + text contains `sub`.
OK: must not raise.  C: a boolean condition (docs, values).
Every group starts with a BASE case that must pass, so a broken harness call cannot pass.
Usage: python docs/dev/measurements/issue160_audit/audit160.py src [out.json]

Issue #160's checklist, one case per sub-case, as run for the bundled fix of issues #155, #160 and
others.  Cases marked "Revised" expect the behaviour decided there rather than the one first
asked for: each says why.  On main before that fix it reports 104 failures; on the fix, none.
"""
import sys, json, warnings, subprocess, os, inspect
sys.path.insert(0, sys.argv[1])
os.environ['MPLBACKEND'] = 'Agg'
import numpy as np
import magnus.oscprob as o, magnus.globaldefs as gd, magnus.earth as ea, magnus.matter as ma
import magnus.hamiltonians as h, magnus.avgprob as ap, magnus.adiabatic as ad, magnus.magnus as mm
import magnus.solarmodels as sm, magnus.oscprobstd as std, magnus.plotting as mp
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
ROOT = os.path.dirname(os.path.abspath(sys.argv[1]))
G = gd.UNIT_GEV; KM = gd.UNIT_KM; E = G; L = 1000*KM; RS = gd.SUN_RADIUS*KM
p = gd.load_nufit_params('NuFIT 6.1'); H0 = h.hamiltonian_3nu_vacuum_energy_independent(**p)
e00 = np.diag([1., 0, 0])
out = []

def _run(fn):
    msg, cls, ws = '', None, []
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter('always')
        try:
            fn()
        except SystemExit as e:
            cls, msg = 'SystemExit', str(e)
        except BaseException as e:
            cls, msg = type(e).__name__, str(e)
        ws = [(x.category.__name__ + ': ' + str(x.message)) for x in w]
    plt.close('all')
    return cls, msg, ws

def R(sec, cid, fn, name):
    cls, msg, ws = _run(fn)
    ok = cls is not None and ('Error in magnus' in msg or 'error:' in msg) and name in msg
    det = 'ACCEPTED' if cls is None else (cls + ': ' + msg[:160])
    out.append(dict(sec=sec, id=cid, want='refused naming ' + repr(name), ok=bool(ok), got=det))

def W(sec, cid, fn, sub):
    cls, msg, ws = _run(fn)
    ok = any(sub in x for x in ws)
    det = ('warned' if ok else 'no matching warning') + ('' if cls is None else '; raised ' + cls + ': ' + msg[:120])
    out.append(dict(sec=sec, id=cid, want='warns ' + repr(sub), ok=bool(ok), got=det))

def OK(sec, cid, fn, check=None):
    cls, msg, ws = _run(lambda: out.append(('__v', fn())))
    v = out.pop() if out and isinstance(out[-1], tuple) else None
    ok = cls is None and (check is None or check(v[1] if v else None))
    det = 'ok' if ok else ('raised ' + str(cls) + ': ' + msg[:150] if cls else 'returned but check failed')
    out.append(dict(sec=sec, id=cid, want='accepted' + ('' if check is None else ' and correct'), ok=bool(ok), got=det))

def C(sec, cid, cond, want):
    try:
        ok = bool(cond())
    except Exception as e:
        ok = False
    out.append(dict(sec=sec, id=cid, want=want, ok=ok, got='yes' if ok else 'no'))

def cli(*args):
    r = subprocess.run([sys.executable, '-m', 'magnus', 'prob', *args], capture_output=True, text=True,
                       env=dict(os.environ, PYTHONPATH=sys.argv[1]))
    if r.returncode != 0:
        raise ValueError('error: ' + r.stderr.strip().splitlines()[-1] if r.stderr.strip() else 'error: exit %d' % r.returncode)
    return r.stdout

vac = lambda **k: o.osc_prob_3nu_vacuum(k.pop('energy', E), k.pop('L', L), **k)
con = lambda **k: o.osc_prob_3nu_matter_constant_density(k.pop('energy', E), k.pop('L', L), k.pop('rho', 3.0), density_matter_is_in_g_per_cm3=k.pop('g', True), **k)
exd = lambda **k: o.osc_prob_3nu_matter_exp_density(k.pop('energy', E), k.pop('L', L), k.pop('L0', 0.0), k.pop('rho_central', 3.0), k.pop('l_scale', 300*KM), density_matter_is_in_g_per_cm3=True, **k)
sun = lambda **k: o.osc_prob_3nu_sun(k.pop('energy', 10e6), k.pop('L', RS), k.pop('L0', 0.0), **k)
eth = lambda **k: o.osc_prob_3nu_earth(k.pop('energy', E), **({'costhz': -0.5} if 'costhz' not in k and 'loc_ini' not in k else {}), **({'L': L} if 'L' not in k and 'loc_ini' not in k else {}), **k)
oeb = lambda Hf=None, **k: o.osc_prob_energy_baseline(H0/E if Hf is None else Hf, k.pop('energy', E), k.pop('L', L), **k)
ex2 = lambda **k: o.osc_prob_2nu_matter_exp_density(3*gd.UNIT_MEV, 6e5*KM, 0.0, rho_central=100*gd.UNIT_G_PER_CM3, l_scale=7e4*KM, sth=0.55, Dm2=7.5e-5, nu_i=0, nu_f=0, **k)
sun2 = lambda **k: o.osc_prob_2nu_sun(10e6, RS, 0.0, sth=0.55, Dm2=7.5e-5, nu_i=0, nu_f=0, **k)
liv3 = lambda **k: o.osc_prob_3nu_vacuum_liv(E, L, **{**dict(sxi12=0.3, sxi23=0.4, sxi13=0.2, dxiCP=0.1, b1=1e-23, b2=2e-23, b3=3e-23, Lambda=1e9, n_liv=0), **k})
nsi3 = lambda **k: o.osc_prob_3nu_matter_nsi_constant_density(E, L, 3.0, density_matter_is_in_g_per_cm3=True, **k)
nsi4 = lambda **k: o.osc_prob_4nu_matter_nsi_constant_density(E, L, 3.0, density_matter_is_in_g_per_cm3=True, **k)
nsi5 = lambda **k: o.osc_prob_5nu_matter_nsi_constant_density(E, L, 3.0, density_matter_is_in_g_per_cm3=True, **k)
liv4 = lambda **k: o.osc_prob_4nu_vacuum_liv(E, L, **{**dict(b1=1e-23, b2=2e-23, b3=3e-23, b4=4e-23, Lambda=1e9, n_liv=0), **k})
liv5 = lambda **k: o.osc_prob_5nu_vacuum_liv(E, L, **{**dict(b1=1e-23, b2=2e-23, b3=3e-23, b4=4e-23, b5=5e-23, Lambda=1e9, n_liv=0), **k})
Harr = lambda l: np.broadcast_to(H0/E, np.shape(l) + (3, 3))

# ---------------- §1
s = '1'
for n, f in [('vac', vac), ('con', con), ('exd', exd), ('sun', sun), ('eth', eth), ('oeb', oeb), ('ex2', ex2), ('liv3', liv3), ('nsi3', nsi3), ('nsi4', nsi4), ('nsi5', nsi5), ('liv4', liv4), ('liv5', liv5)]:
    OK(s, 'BASE ' + n, f)
R(s, 'type named: s12="a"', lambda: vac(s12='a'), 's12')
for arg, val, f in [('energy', np.int64(10**9), vac), ('energy', np.float32(1e9), vac), ('L', np.float32(L), vac), ('L', np.array(L), vac),
                    ('s12', np.float32(0.55), vac), ('dCP', np.int64(3), vac), ('D21', np.float32(7.5e-5), vac), ('rho', np.float32(3.0), con), ('L0', np.int64(0), exd)]:
    OK(s, 'numpy scalar %s %s' % (arg, type(val).__name__), lambda arg=arg, val=val, f=f: f(**{arg: val}))
for arg, f in [('energy', vac), ('L', vac), ('L0', exd), ('s12', vac), ('dCP', vac), ('D21', vac), ('rho', con), ('electron_fraction', con),
               ('eps_ee', nsi3), ('b1', liv3), ('Lambda', liv3), ('nu_i', vac), ('nu_f', vac)]:
    R(s, 'True refused: ' + arg, lambda arg=arg, f=f: f(**({arg: True, 'nu_f': 0} if arg == 'nu_i' else {arg: True, 'nu_i': 0} if arg == 'nu_f' else {arg: True})), arg)
for arg, val, f in [('L', np.nan, vac), ('L', np.inf, vac), ('dCP', np.nan, vac), ('dCP', np.inf, vac), ('D21', np.nan, vac), ('D31', np.inf, vac),
                    ('eps_em', np.nan, nsi3), ('eps_em', np.inf, nsi3), ('b1', np.nan, liv3), ('b2', np.inf, liv3), ('Lambda', np.nan, liv3), ('Lambda', np.inf, liv3),
                    ('electron_fraction', np.nan, con), ('ratio_number_neutrons_to_protons', np.nan, con), ('ratio_number_neutrons_to_protons', np.inf, con)]:
    R(s, 'finite: %s=%s' % (arg, val), lambda arg=arg, val=val, f=f: f(**{arg: val}), arg)
R(s, 'energy nan own message', lambda: vac(energy=np.nan), 'finite')
for n, f in [('vac', vac), ('con', con), ('exd', exd)]:
    R(s, 'energy array inf entry: ' + n, lambda f=f: f(energy=np.array([1e9, np.inf]), L=np.array([L, L])), 'energy')
R(s, 'energy array inf entry: sun', lambda: sun(energy=np.array([1e7, np.inf]), L=np.array([RS, RS])), 'energy')
R(s, 'energy inf: sun average', lambda: sun(energy=np.inf, average=True), 'energy')
for n, f in [('vac', vac), ('con', con), ('exd', exd)]:
    R(s, 'L array nan entry: ' + n, lambda f=f: f(energy=np.array([1e9, 2e9]), L=np.array([L, np.nan])), 'L')
R(s, 'L array nan entry: sun', lambda: sun(energy=np.array([1e7, 2e7]), L=np.array([RS, np.nan])), 'L')
R(s, 'L array nan entry: earth', lambda: eth(energy=np.array([1e9, 2e9]), L=np.array([L, np.nan])), 'L')
for n, f in [('vac', vac), ('con', con)]:
    R(s, 'L array inf entry: ' + n, lambda f=f: f(energy=np.array([1e9, 2e9]), L=np.array([L, np.inf])), 'L')
R(s, 'L array negative entry named', lambda: vac(energy=np.array([1e9, 2e9]), L=np.array([L, -L])), 'L')
R(s, 'L < L0 named with values', lambda: exd(L=0.5*L, L0=L), 'L0')
R(s, 'L array: first offending index named', lambda: exd(energy=np.array([1e9, 2e9]), L=np.array([2*L, 0.5*L]), L0=L), 'entry 1')
R(s, 'Sun average=True with L0 > L', lambda: sun(L=0.1*RS, L0=0.5*RS, average=True), 'L')
R(s, 'Sun scan with some L < L0', lambda: sun(energy=np.array([1e7, 1e7]), L=np.array([1e5*KM, 3e5*KM]), L0=2e5*KM), 'L')
R(s, 'Sun L0 < 0', lambda: sun(L0=-1e5*KM), 'L0')
R(s, 'L0 array: error names the wrapper', lambda: sun(L0=np.array([0.0, 1.0])), 'osc_prob_3nu_sun')
C(s, 'L0 docstring says scalar; loop over production points', lambda: 'loop over' in (o.osc_prob_3nu_sun.__doc__ or ''), 'docstring note')
C(s, 'solar_models.rst says L0 is scalar', lambda: 'scalar' in open(os.path.join(ROOT, 'docs/source/solar_models.rst')).read() and 'L0' in open(os.path.join(ROOT, 'docs/source/solar_models.rst')).read(), 'docs note')
R(s, 'empty arrays refused by name', lambda: vac(energy=np.array([]), L=np.array([])), 'empty')
def _same_len1():
    a = _run(lambda: vac(energy=np.array([1e9]), L=np.array([L, 2*L])))[0]
    b = _run(lambda: oeb(energy=np.array([1e9]), L=np.array([L, 2*L])))[0]
    return (a is None) == (b is None)
C(s, 'length-1 energy vs several L: wrappers and oeb alike', _same_len1, 'same behaviour')
R(s, 'oeb L=True', lambda: oeb(L=True), 'L')
R(s, 'oeb 2-D L', lambda: oeb(L=np.ones((2, 2))*L), 'L')
R(s, 'oeb masked energy', lambda: oeb(energy=np.ma.masked_array([1e9, 2e9], [0, 1])), 'mask')
R(s, 'oeb (H, [], [])', lambda: oeb(energy=[], L=[]), 'empty')
R(s, 'oeb 2-D energy with scalar L', lambda: oeb(energy=np.ones((2, 2))*1e9), 'energy')
R(s, 'rho negative named', lambda: con(rho=-3.0), 'rho')
R(s, 'rho nan named', lambda: con(rho=np.nan), 'rho')
R(s, 'rho_central negative named', lambda: exd(rho_central=-3.0), 'rho_central')
R(s, 'l_scale 0 named', lambda: exd(l_scale=0.0), 'l_scale')
R(s, 'l_scale nan named', lambda: exd(l_scale=np.nan), 'l_scale')
R(s, 'L0 nan named', lambda: exd(L0=np.nan), 'L0')
R(s, 'costhz nan named', lambda: eth(costhz=np.nan, L=L), 'costhz')
R(s, 'density_matter_ocean -1 named', lambda: eth(density_matter_ocean=-1.0), 'density_matter_ocean')
R(s, 'CLI --rho -3 names --rho', lambda: cli('--environment', 'matter', '--energy', '1', '--baseline', '1000', '--rho', '-3'), '--rho')
R(s, 'electron_fraction 1.5', lambda: con(electron_fraction=1.5), 'electron_fraction')
R(s, "electron_fraction 'a'", lambda: con(electron_fraction='a'), 'electron_fraction')
for k in ('core', 'mantle', 'crust', 'ocean'):
    R(s, 'electron_fraction_%s 1.5' % k, lambda k=k: eth(**{'electron_fraction_' + k: 1.5}), 'electron_fraction_' + k)
R(s, 'electron_fraction_mantle 0 refused', lambda: eth(electron_fraction_mantle=0.0), 'electron_fraction_mantle')
R(s, 'ratio n/p negative', lambda: con(ratio_number_neutrons_to_protons=-1.0), 'ratio_number_neutrons_to_protons')
def _callable_ratio_earth():
    cls, msg, ws = _run(lambda: eth(ratio_number_neutrons_to_protons=lambda l: np.nan))
    return cls is not None and 'ratio_number_neutrons_to_protons' in msg
C(s, 'callable ratio on Earth: refused (or used)', _callable_ratio_earth, 'refused naming the argument')
R(s, 'num_flavors 6 without H', lambda: o.osc_prob_vacuum(6, E, L, p), 'h_vac_energy_indep')
R(s, 'num_flavors 3.0', lambda: o.osc_prob_vacuum(3.0, E, L, p), 'num_flavors')
R(s, 'num_flavors "3"', lambda: o.osc_prob_vacuum('3', E, L, p), 'num_flavors')
def _hvac():
    H4 = np.diag([0, 7e-5, 2.5e-3, 1.0]).astype(complex)
    cls, msg, ws = _run(lambda: o.osc_prob_vacuum(4, E, L, dict(p, s14=.1, d14=.3, s24=.15, d24=.5, s34=.2, D41=1.), h_vac_energy_indep=H4))
    return cls is not None and 'h_vac_energy_indep' in msg
C(s, 'h_vac_energy_indep at 4 flavors: used or refused', _hvac, 'refused naming it (or used)')
R(s, 'no effect: t_breakpoints on vacuum', lambda: vac(t_breakpoints=[0.5*L]), 't_breakpoints')
# Revised: cumulative=True is not ignored on a constant density; it runs the cumulative engine.
def _cumulative_served():
    info = {}
    con(cumulative=True, L0=0.0, strategy_info=info)
    return info.get('engine') == 'cumulative'
C(s, 'cumulative=True on constant density is served by the cumulative engine', _cumulative_served, 'engine cumulative')
R(s, 'no effect: electron_fraction on tabulated Sun', lambda: sun(density_profile='B16-GS98', electron_fraction=0.5), 'electron_fraction')
R(s, 'no effect: ratio n/p on tabulated Sun', lambda: sun(density_profile='B16-GS98', ratio_number_neutrons_to_protons=1.0), 'ratio_number_neutrons_to_protons')
R(s, "no effect: integration_method='foo' on vacuum", lambda: vac(integration_method='foo'), 'integration_method')
# Revised (decision): rtol/atol stay accepted in vacuum -- the CLI and shared calls forward one
# set to every wrapper -- and the vacuum docstrings say they have no effect there.
OK(s, 'rtol on vacuum accepted', lambda: vac(rtol=1e-6))
C(s, 'vacuum docstring: refinement keywords have no effect', lambda: 'have no effect' in ' '.join((o.osc_prob_3nu_vacuum.__doc__ or '').split()), 'documented')
R(s, 'no effect: average_n_samples without average', lambda: vac(average_n_samples=11), 'average_n_samples')
R(s, 'L0= on vacuum: standard message', lambda: vac(L0=0.0), 'L0')
R(s, 'L0= on Earth: standard message', lambda: eth(L0=0.0), 'L0')
W(s, 't_breakpoints in raw km warns', lambda: exd(t_breakpoints=[500.0]), 'BaselineUnitWarning')
W(s, 'energy 10 eV warns (Sun)', lambda: sun(energy=10.0), 'EnergyUnitWarning')

# ---------------- §2
s = '2'
R(s, 'dcp typo, did-you-mean dCP', lambda: o.osc_prob_vacuum(3, E, L, dict(p, dcp=0.5)), 'dCP')
R(s, 's14 at 3 flavors', lambda: o.osc_prob_vacuum(3, E, L, dict(p, s14=0.1)), 's14')
OK(s, 'name/description allowed', lambda: o.osc_prob_vacuum(3, E, L, dict(p, name='x', description='y')))
R(s, 'non-Mapping osc_params', lambda: o.osc_prob_vacuum(3, E, L, [1, 2]), 'osc_params')
R(s, 'NSI dict unknown key', lambda: o.osc_prob_matter_nsi(3, 3.0, E, L, p, dict(eps_ee=0.1, eps_xx=0.1), density_matter_is_in_g_per_cm3=True), 'eps_xx')
R(s, 'LIV dict unknown key', lambda: o.osc_prob_liv(3, E, L, p, dict(sxi12=0.3, sxi23=0.4, sxi13=0.2, dxiCP=0.1, b1=1e-23, b2=2e-23, b3=3e-23, Lambda=1e9, n_liv=0, bx=1)), 'bx')
for f, n, args in [(nsi4, '4', ['eps_em', 'eps_ss', 'eps_ts']), (nsi5, '5', ['eps_em', 'eps_s1s1', 'eps_s1s2'])]:
    for a in args:
        R(s, "'a' in %s at %s flavors" % (a, n), lambda f=f, a=a: f(**{a: 'a'}), a)
for f, n, args in [(liv4, '4', ['b4', 'sxi14', 'dxi24', 'Lambda', 'n_liv']), (liv5, '5', ['b5', 'sxi15', 'dxi35', 'Lambda', 'n_liv'])]:
    for a in args:
        R(s, "'a' in %s at %s flavors" % (a, n), lambda f=f, a=a: f(**{a: 'a'}), a)
for f, n, a in [(nsi3, '3', 'eps_ee'), (nsi4, '4', 'eps_ss'), (nsi5, '5', 'eps_s1s1'), (nsi5, '5', 'eps_s2s2')]:
    R(s, 'complex diagonal %s at %s' % (a, n), lambda f=f, a=a: f(**{a: 0.1 + 0.1j}), a)
OK(s, 'complex diagonal with zero imaginary part accepted', lambda: nsi3(eps_ee=0.1 + 0j))
for a in ('b1', 'sxi12', 'dxiCP'):
    R(s, 'complex LIV %s' % a, lambda a=a: liv3(**{a: 0.1 + 0.1j if a != 'b1' else 1e-23 + 1e-23j}), a)
for f, n, a in [(nsi4, '4', 'eps_es'), (nsi5, '5', 'eps_ts2'), (liv4, '4', 'b4'), (liv4, '4', 'dxi14'), (liv5, '5', 'Lambda'), (liv5, '5', 'b5')]:
    R(s, 'nan %s at %s' % (a, n), lambda f=f, a=a: f(**{a: np.nan}), a)
    R(s, 'inf %s at %s' % (a, n), lambda f=f, a=a: f(**{a: np.inf}), a)
R(s, 'sxi14 nan named at 4', lambda: liv4(sxi14=np.nan), 'sxi14')
R(s, 'Lambda 0', lambda: liv3(Lambda=0.0), 'Lambda')
for v in (-1, 1.5, np.nan, np.inf, 1j, True):
    R(s, 'n_liv=%r' % (v,), lambda v=v: liv3(n_liv=v), 'n_liv')

# ---------------- §3
s = '3'
for n in ('2', '3', '5'):
    f = getattr(o, 'osc_prob_%snu_earth' % n)
    kw = dict(sth=0.55, Dm2=7.5e-5) if n == '2' else {}
    OK(s, 'BASE earth %s' % n, lambda f=f, kw=kw: f(E, costhz=-0.5, L=L, **kw))
    R(s, 'costhz -1.5 with L at %s flavors' % n, lambda f=f, kw=kw: f(E, costhz=-1.5, L=L, **kw), 'costhz')
R(s, 'costhz 2.0 with L=1e4 km', lambda: eth(costhz=2.0, L=1e4*KM), 'costhz')
R(s, 'L beyond the chord: names L, costhz, chord', lambda: eth(costhz=-0.8, L=12000*KM), 'chord')
C(s, 'shorter L than chord: documented as a partial path', lambda: 'partial' in (o.osc_prob_3nu_earth.__doc__ or ''), 'docstring')
fer, hom = ea.loc_coords_dms['fermilab'], ea.loc_coords_dms['homestake']
OK(s, 'BASE loc pair', lambda: eth(loc_ini='fermilab', loc_fin='homestake'))
R(s, 'loc latitude 146', lambda: eth(loc_ini=((146, 0, 0), (0, 0, 0)), loc_fin=((10, 0, 0), (10, 0, 0))), 'lat')
R(s, 'loc decimal pair', lambda: eth(loc_ini=(41.8, -88.3), loc_fin=(44.3, -103.8)), 'loc_ini')
OK(s, 'BASE chord helper', lambda: ea.chord_length_inside_earth(fer['lat'], fer['lon'], hom['lat'], hom['lon']))
for lab, args, nm in [('lat 100', ((100, 0, 0), fer['lon'], hom['lat'], hom['lon']), 'lat'), ('lon 400', (fer['lat'], (400, 0, 0), hom['lat'], hom['lon']), 'lon'),
                      ('minutes 70', ((41, 70, 0), fer['lon'], hom['lat'], hom['lon']), 'minutes'), ('seconds -1', ((41, 0, -1), fer['lon'], hom['lat'], hom['lon']), 'lat1_dms'),
                      ('nan', ((np.nan, 0, 0), fer['lon'], hom['lat'], hom['lon']), 'finite'), ('(d, m) pair', ((41, 0), fer['lon'], hom['lat'], hom['lon']), 'lat')]:
    R(s, 'chord_length_inside_earth ' + lab, lambda args=args: ea.chord_length_inside_earth(*args), nm)
    R(s, 'costhz_between_points_on_surface ' + lab, lambda args=args: ea.costhz_between_points_on_surface(*args), nm)
OK(s, 'BASE dms_to_decimal', lambda: ea.dms_to_decimal(41, 50, 0))
for lab, args, nm in [('degrees 500', (500, 0, 0), 'degrees'), ('minutes 70', (10, 70, 0), 'minutes'), ('minutes -1', (10, -1, 0), 'minutes'), ('seconds 70', (10, 0, 70), 'seconds'),
                      ('nan', (np.nan, 0, 0), 'degrees'), ('inf', (np.inf, 0, 0), 'degrees'), ('True', (True, 0, 0), 'degrees')]:
    R(s, 'dms_to_decimal ' + lab, lambda args=args: ea.dms_to_decimal(*args), nm)
OK(s, 'BASE distance_traveled_inside_earth', lambda: ea.distance_traveled_inside_earth(-0.5))
for lab, v in [('-1.5', -1.5), ('1.0000001', 1.0000001), ('nan', np.nan), ("'a'", 'a')]:
    R(s, 'distance_traveled_inside_earth ' + lab, lambda v=v: ea.distance_traveled_inside_earth(v), 'costhz')
cls_, msg_, _ = _run(lambda: ea.distance_traveled_inside_earth(np.array([-0.5, -0.2])))
C(s, 'distance_traveled_inside_earth array: vectorized or refused by name', lambda: cls_ is None or 'costhz' in msg_, 'no ambiguous-truth error')
R(s, 'prem_layer_edges_along_chord -1.5', lambda: ea.prem_layer_edges_along_chord(-1.5), 'costhz')
R(s, 'prem_layer_edges_along_chord detector_depth True', lambda: ea.prem_layer_edges_along_chord(-0.5, detector_depth=True), 'detector_depth')
R(s, "prem_layer_edges_along_chord source_depth 'a'", lambda: ea.prem_layer_edges_along_chord(-0.5, source_depth='a'), 'source_depth')
OK(s, 'BASE earth_radial_distance_from_depth', lambda: ea.earth_radial_distance_from_depth(-0.5, 1000.0))
for lab, args, nm in [('costhz -1.5', (-1.5, 1000.0), 'costhz'), ('costhz nan', (np.nan, 1000.0), 'costhz'), ('l nan', (-0.5, np.nan), 'l'),
                      ('l negative', (-0.5, -10.0), 'l'), ('l True', (-0.5, True), 'l')]:
    R(s, 'earth_radial_distance_from_depth ' + lab, lambda args=args: ea.earth_radial_distance_from_depth(*args), nm)
R(s, 'coordinates_of_named_location None', lambda: ea.coordinates_of_named_location('x', None), 'loc')
R(s, 'coordinates_of_named_location 1', lambda: ea.coordinates_of_named_location('x', 1), 'loc')
R(s, 'canonical_name None', lambda: sm.canonical_name(None), 'name')
R(s, 'table_edge 1', lambda: sm.table_edge(1), 'name')
for fn in ('load_solar_model', 'solar_model_info', 'electron_density_profile', 'neutron_to_proton_ratio_profile'):
    R(s, fn + '(None) TypeError', lambda fn=fn: getattr(sm, fn)(None), 'name')
OK(s, "load_nufit_params case-insensitive", lambda: gd.load_nufit_params('nufit 6.1'))

# ---------------- §4
s = '4'
flags = [('nubar', vac), ('average', vac), ('density_matter_is_in_g_per_cm3', lambda **k: o.osc_prob_3nu_matter_constant_density(E, L, 3.0, **k)),
         ('density_is_of_number_of_electrons', lambda **k: o.osc_prob_3nu_matter_constant_density(E, L, 3.0, **k)),
         ('return_evolution_operator', exd), ('strict_convergence', exd), ('validate_input', vac), ('save_log', vac), ('close_file_log_upon_exit', vac),
         ('H_func_is_function_only_of_energy', oeb)]
for a, f in flags:
    for v in ('False', 2, None, np.array([True, False])):
        R(s, '%s=%r' % (a, v if not isinstance(v, np.ndarray) else 'array'), lambda a=a, f=f, v=v: f(**{a: v}), a)
for v in ('False', 2, np.array([True, False])):
    R(s, 'compute_matrix_multiplication=%r' % (v if not isinstance(v, np.ndarray) else 'array'), lambda v=v: h.hamiltonian_3nu_vacuum_energy_independent(**p, compute_matrix_multiplication=v), 'compute_matrix_multiplication')
for v in ('no', 2, None):
    R(s, 'set_color_output(%r)' % (v,), lambda v=v: gd.set_color_output(v), 'enabled')
gd.set_color_output(False)
OK(s, 'flags accept np.bool_', lambda: vac(nubar=np.bool_(True)))
for v in (1.0, True, 3, -1):
    R(s, 'nu_i=%r' % (v,), lambda v=v: vac(nu_i=v, nu_f=0), 'nu_i')
OK(s, 'nu_i np.int64 accepted', lambda: vac(nu_i=np.int64(1), nu_f=0))
for v in (True, 2.5, np.nan, 'a'):
    R(s, 'verbose=%r' % (v,), lambda v=v: vac(verbose=v), 'verbose')
for v in (0, -1, 2.5, 'a'):
    R(s, 'new_recursion_limit=%r' % (v,), lambda v=v: o.osc_prob(lambda l: H0/E, 0, L, new_recursion_limit=v), 'new_recursion_limit')
for v in (2, 1, None):
    R(s, 'cumulative=%r' % (v,), lambda v=v: exd(cumulative=v), 'cumulative')
R(s, "integration_method='foo' on exp profile", lambda: exd(integration_method='foo'), 'integration_method')
R(s, "average_initial_state='bogus'", lambda: vac(average=True, average_initial_state='bogus'), 'average_initial_state')
R(s, 'filename_log=3', lambda: vac(save_log=True, filename_log=3), 'filename_log')
R(s, 'strategy_info=[]', lambda: exd(strategy_info=[]), 'strategy_info')
R(s, 'hybrid_propagator info=[]', lambda: ad.hybrid_propagator(lambda l: H0/E, 0, L, info=[]), 'info')
R(s, 'default_osc_params_set_name=None', lambda: vac(default_osc_params_set_name=None), 'default_osc_params_set_name')
R(s, 'save_log to missing directory', lambda: vac(save_log=True, filename_log='/nonexistent_dir_xyz/out.log'), 'filename_log')
R(s, 'save_log to a directory', lambda: vac(save_log=True, filename_log='/tmp'), 'filename_log')
cls_, msg_, _ = _run(lambda: exd(strategy='foo'))
C(s, 'strategy error names the wrapper, not an internal function', lambda: 'osc_prob_3nu_matter_exp_density' in msg_ and 'strategy' in msg_, 'wrapper named')

# ---------------- §5
s = '5'
OK(s, 'BASE ex2', ex2); OK(s, 'BASE sun2', sun2)
for lab, kw in [('max_n_slabs=-1', dict(max_n_slabs=-1)), ('min_n_slabs=-5', dict(min_n_slabs=-5)), ('max_n_slabs=0', dict(max_n_slabs=0)), ('n_slabs=-3', dict(n_slabs=-3))]:
    nm = list(kw)[0]
    R(s, 'engine-independent %s on 2nu exp' % lab, lambda kw=kw: ex2(**kw), nm)
    R(s, 'engine-independent %s on 2nu Sun' % lab, lambda kw=kw: sun2(**kw), nm)
for args, nm in [(('--rtol', '0'), '--rtol'), (('--rtol', '-1'), '--rtol'), (('--rtol', 'nan'), '--rtol'), (('--atol', '0'), '--atol'), (('--atol', 'nan'), '--atol'),
                 (('--n-jobs', '0'), '--n-jobs'), (('--n-jobs', '-5'), '--n-jobs'), (('--magnus-exp-order', '11'), '--magnus-exp-order')]:
    R(s, 'CLI ' + ' '.join(args), lambda args=args: cli('--environment', 'matter', '--density-profile', 'exp', '--energy', '1', '--baseline', '1000', '--rho-central', '3', '--l-scale', '300', *args), nm)
for lab, kw, nm in [('rtol=nan', dict(rtol=np.nan), 'rtol'), ('atol=inf', dict(atol=np.inf), 'atol'), ('rtol=True', dict(rtol=True), 'rtol'), ('rtol=0, atol=0', dict(rtol=0, atol=0), 'rtol')]:
    R(s, lab, lambda kw=kw: exd(**kw), nm)
C(s, 'one-sided None documented', lambda: 'only one of the two is ``None``' in (o.osc_prob.__doc__ or ''), 'documented in osc_prob')
W(s, 'impossible tolerance warns (1e-14)', lambda: exd(rtol=1e-14, atol=1e-14, max_n_slabs=64), 'Warning')
for lab, kw, nm in [('n_slabs=2.5', dict(n_slabs=2.5), 'n_slabs'), ('n_slabs=nan', dict(n_slabs=np.nan), 'n_slabs'), ('n_slabs="10"', dict(n_slabs='10'), 'n_slabs'),
                    ('max_num_loops=2.5', dict(max_num_loops=2.5), 'max_num_loops'), ('n_tpts_per_slab=inf', dict(n_tpts_per_slab=np.inf), 'n_tpts_per_slab'),
                    ('max_n_slabs=0 off the ladder', dict(max_n_slabs=0, rtol=None, atol=None, n_slabs=4), 'max_n_slabs')]:
    R(s, lab, lambda kw=kw: exd(**kw), nm)
R(s, 'max_n_slabs=-1 names max_n_slabs as positive integer', lambda: exd(max_n_slabs=-1), 'max_n_slabs must')
C(s, 'n_slabs > max_n_slabs: documented clip', lambda: 'max_n_slabs' in (o.osc_prob.__doc__ or '') and 'clip' in (o.osc_prob.__doc__ or '').lower(), 'documented')
def _ntpts_clip():
    cls, msg, ws = _run(lambda: exd(integration_method='simpson', n_tpts_per_slab=601, max_n_tpts_per_slab=501, rtol=None, atol=None, n_slabs=4))
    return (cls is not None and 'n_tpts_per_slab' in msg) or 'clip' in (o.osc_prob.__doc__ or '').lower() and 'max_n_tpts_per_slab' in (o.osc_prob.__doc__ or '') and cls is None
C(s, 'n_tpts_per_slab > max_n_tpts_per_slab: documented clip or refused', _ntpts_clip, 'decided rule')
for lab, kw, nm in [('growth_factor_n_slabs=1.0', dict(growth_factor_n_slabs=1.0), 'growth_factor_n_slabs'), ('growth_factor_n_tpts_per_slab=1.0', dict(growth_factor_n_tpts_per_slab=1.0), 'growth_factor_n_tpts_per_slab'),
                    ('growth_factor_n_slabs=nan', dict(growth_factor_n_slabs=np.nan), 'growth_factor_n_slabs'), ('magnus_exp_order=2.5', dict(magnus_exp_order=2.5), 'magnus_exp_order'),
                    ('magnus_exp_order=True', dict(magnus_exp_order=True), 'magnus_exp_order'), ("gl with odd order 3", dict(integration_method='gl', magnus_exp_order=3), 'magnus_exp_order'),
                    ]:
    R(s, lab, lambda kw=kw: exd(**kw), nm)
# Revised (decision): points-per-slab settings under 'gl' warn that they do nothing, and are
# documented; they are not refused.
for lab, kw in [("gl with n_tpts_per_slab", dict(integration_method='gl', n_tpts_per_slab=600)),
                ("gl with max_n_tpts_per_slab", dict(integration_method='gl', max_n_tpts_per_slab=100)),
                ("gl with growth_factor_n_tpts_per_slab", dict(integration_method='gl', growth_factor_n_tpts_per_slab=2.0))]:
    W(s, lab + ' warns', lambda kw=kw: exd(**kw), 'IgnoredQuadratureSettingWarning')
OK(s, 'default integration method with default order still works', lambda: exd())
for v in (0, -5, 2.5, True, None, 'a'):
    R(s, 'n_jobs=%r' % (v,), lambda v=v: exd(energy=np.array([1e9, 2e9]), L=np.array([L, L]), n_jobs=v), 'n_jobs')
C(s, 'n_jobs docstring: loky pool persists until idle timeout', lambda: 'idle' in (o.osc_prob_energy_baseline.__doc__ or '').lower(), 'documented')
R(s, 'average=True with cumulative=True', lambda: exd(average=True, cumulative=True), 'cumulative')
cls_, msg_, _ = _run(lambda: exd(return_evolution_operator=True, nu_i=0, nu_f=0))
C(s, 'return_evolution_operator with nu_i/nu_f: refused by name or works', lambda: cls_ is None or ('return_evolution_operator' in msg_), 'no inhomogeneous-shape error')

# ---------------- §6
s = '6'
edges = lambda e: o.osc_prob(lambda l: H0/E, 0, L, t_slab_edges=e, rtol=None, atol=None)
OK(s, 'BASE t_slab_edges', lambda: edges([[0, 0.5*L], [0.5*L, L]]))
for lab, e in [('gap', [[0, 0.4*L], [0.5*L, L]]), ('overlap', [[0, 0.6*L], [0.5*L, L]]), ('short of L', [[0, 0.5*L]]), ('flat [0, 2L]', [0, 2*L]), ('nan', [[0, np.nan], [np.nan, L]]), ('[]', [])]:
    R(s, 't_slab_edges ' + lab, lambda e=e: edges(e), 't_slab_edges')
A = lambda t: -1j*H0/E
R(s, 't_slab_edges zero-width slab', lambda: edges([[0, 0.5*L], [0.5*L, 0.5*L], [0.5*L, L]]), 't_slab_edges')
OK(s, 'core accepts a zero-width slab (documented departure)', lambda: mm.magnus_expansion_multislab(A, [[0, 0.5*L], [0.5*L, 0.5*L], [0.5*L, L]]))
A = lambda t: -1j*H0/E
OK(s, 'BASE core', lambda: mm.magnus_expansion_multislab(A, [[0, L]]))
R(s, 'core t_slab_edges [0, 1, 2]', lambda: mm.magnus_expansion_multislab(A, [0, 1, 2]), 't_slab_edges')
R(s, 'core t_slab_edges zeros((0, 2))', lambda: mm.magnus_expansion_multislab(A, np.zeros((0, 2))), 't_slab_edges')
for lab, tb in [('nan', [np.nan]), ('2-D', [[0.5*L]]), ('outside [L0, L]', [2*L]), ('string', 'a')]:
    R(s, 't_breakpoints ' + lab, lambda tb=tb: exd(t_breakpoints=tb), 't_breakpoints')

# ---------------- §7
s = '7'
for lab, Hf in [('None', None), ("'abc'", 'abc'), ('1-D array', np.ones(3)), ('lambda: H0', lambda: H0)]:
    R(s, 'H_func ' + lab, lambda Hf=Hf: o.osc_prob_energy_baseline(Hf, E, L), 'H_func')
R(s, 'first sample list: oeb', lambda: oeb(Hf=lambda E_, l: (H0/E_).tolist()), 'H_func')
R(s, 'first sample list: osc_prob_sun', lambda: o.osc_prob_sun(lambda E_, l, V: (H0/E_).tolist(), 1e7, RS), 'H_func')
R(s, 'first sample object array', lambda: oeb(Hf=lambda E_, l: np.array(H0/E_, dtype=object)), 'H_func')
R(s, 'first sample NaN', lambda: oeb(Hf=lambda E_, l: H0/E_ + np.nan*e00), 'H_func')
R(s, 'hybrid_propagator NaN beyond l=5', lambda: ad.hybrid_propagator(lambda l: np.array([[l - 5, .3], [.3, 5 - l]], dtype=complex) if l < 5 else np.full((2, 2), np.nan), 0, 10), 'H_func')
R(s, 'non-Hermitian: triangular', lambda: oeb(Hf=np.triu(H0 + 1e-3)/E), 'Hermitian')
R(s, 'non-Hermitian: 1j*diag', lambda: oeb(Hf=1j*np.diag([1., 2., 3.])*1e-13), 'Hermitian')
R(s, 'non-Hermitian with return_evolution_operator', lambda: oeb(Hf=np.triu(H0 + 1e-3)/E, return_evolution_operator=True), 'Hermitian')
R(s, 'non-Hermitian in averaged_probabilities_constant_hamiltonian', lambda: ap.averaged_probabilities_constant_hamiltonian(np.triu(H0 + 1e-3)/E, L), 'Hermitian')
R(s, 'non-Hermitian in hybrid_propagator', lambda: ad.hybrid_propagator(lambda l: np.triu(H0 + 1e-3)/E, 0, L), 'Hermitian')
def _c64():
    f64 = lambda E_, l, V: H0/E_ + np.asarray(V)[..., None, None]*e00
    f32 = lambda E_, l, V: (H0/E_ + np.asarray(V)[..., None, None]*e00).astype(np.complex64)
    a = o.osc_prob_sun(f64, 5e6, RS, nu_i=0, nu_f=0); b = o.osc_prob_sun(f32, 5e6, RS, nu_i=0, nu_f=0)
    return abs(a - b) < 1e-3
C(s, 'complex64 H on the Sun within 1e-3', _c64, 'within rtol')
mid = 500*KM
OK(s, 'rho_func returning 0-d array accepted', lambda: o.osc_prob_matter_std_potential(3, lambda l: np.array(3.0), E, L, p, density_matter_is_in_g_per_cm3=True))
OK(s, 'rho_func returning np.float32 accepted', lambda: o.osc_prob_matter_std_potential(3, lambda l: np.float32(3.0), E, L, p, density_matter_is_in_g_per_cm3=True))
def _where_fast():
    cls, msg, ws = _run(lambda: o.osc_prob_matter_std_potential(3, lambda l: np.where(np.asarray(l) < mid, 3.0, 8.0), E, L, p, t_breakpoints=[mid], density_matter_is_in_g_per_cm3=True))
    return cls is None and not any('ScalarHamiltonianWarning' in x for x in ws)
C(s, 'vectorized np.where rho_func: fast path, no ScalarHamiltonianWarning', _where_fast, 'no warning')
W(s, 'scalar-only rho_func: warning names rho_func', lambda: o.osc_prob_matter_std_potential(3, lambda l: 3.0 if l < mid else 8.0, E, L, p, t_breakpoints=[mid], density_matter_is_in_g_per_cm3=True), 'rho_func')
for lab, rf in [('None', lambda l: None), ("'a'", lambda l: 'a'), ('complex', lambda l: 3.0 + 1j)]:
    R(s, 'rho_func returns ' + lab, lambda rf=rf: o.osc_prob_matter_std_potential(3, rf, E, L, p, density_matter_is_in_g_per_cm3=True), 'rho_func')
R(s, 'rho_func NaN partway', lambda: o.osc_prob_matter_std_potential(3, lambda l: np.nan if np.all(np.asarray(l) > 0.3*L) else 3.0, E, L, p, density_matter_is_in_g_per_cm3=True), 'rho_func')
cls_, msg_, ws_ = _run(lambda: oeb(Hf=lambda E_: H0/E_ + 1e-13*e00, energy=np.array([1e9, 2e9]), L=L))
C(s, 'H(E) read as H(l): warned or refused, pointing to the flag', lambda: 'H_func_is_function_only_of_energy' in msg_ or any('H_func_is_function_only_of_energy' in x for x in ws_), 'flag named')

# ---------------- §8
s = '8'
OK(s, 'BASE averaged_probabilities_numerically', lambda: ap.averaged_probabilities_numerically(lambda e: np.eye(3)*0.5, 1.0, 0.1, 5))
for lab, args, nm in [('energy 0', (lambda e: 0.5, 0.0, 0.1, 5), 'energy'), ('energy -1', (lambda e: 0.5, -1.0, 0.1, 5), 'energy'), ('energy nan', (lambda e: 0.5, np.nan, 0.1, 5), 'energy'),
                      ('energy inf', (lambda e: 0.5, np.inf, 0.1, 5), 'energy'), ('energy True', (lambda e: 0.5, True, 0.1, 5), 'energy'), ('n_samples 2.5', (lambda e: 0.5, 1.0, 0.1, 2.5), 'n_samples'),
                      ('n_samples nan', (lambda e: 0.5, 1.0, 0.1, np.nan), 'n_samples'), ('non-callable', (3, 1.0, 0.1, 5), 'prob_of_energy')]:
    R(s, 'averaged_probabilities_numerically ' + lab, lambda args=args: ap.averaged_probabilities_numerically(*args), nm)
OK(s, 'BASE averaged_probabilities_constant_hamiltonian', lambda: ap.averaged_probabilities_constant_hamiltonian(H0/E, L))
for lab, args, nm in [('NaN H', (H0/E + np.nan, L), 'hamiltonian'), ('0x0', (np.zeros((0, 0)), L), 'hamiltonian'), ('1-D', (np.ones(3), L), 'hamiltonian'), ("'a'", ('a', L), 'hamiltonian'),
                      ('non-Hermitian', (np.triu(H0 + 1e-3)/E, L), 'hamiltonian'), ('baseline -L', (H0/E, -L), 'baseline'), ('baseline nan', (H0/E, np.nan), 'baseline'),
                      ('baseline inf', (H0/E, np.inf), 'baseline'), ('baseline complex', (H0/E, 1j), 'baseline'), ('baseline True', (H0/E, True), 'baseline')]:
    R(s, 'averaged_probabilities_constant_hamiltonian ' + lab, lambda args=args: ap.averaged_probabilities_constant_hamiltonian(*args), nm)
OK(s, 'BASE phase_averaged', lambda: ap.phase_averaged_probabilities_constant_hamiltonian(H0/E, -H0/E, L, 0.1))
for lab, args, nm in [('non-Hermitian dH', (H0/E, np.triu(H0 + 1e-3)/E, L, 0.1), 'dH_dlnE'), ('NaN dH', (H0/E, H0/E + np.nan, L, 0.1), 'dH_dlnE'), ('baseline -L', (H0/E, -H0/E, -L, 0.1), 'baseline'),
                      ('spread nan', (H0/E, -H0/E, L, np.nan), 'spread'), ('spread inf', (H0/E, -H0/E, L, np.inf), 'spread'), ('spread True', (H0/E, -H0/E, L, True), 'spread'), ('spread -0.1', (H0/E, -H0/E, L, -0.1), 'spread')]:
    R(s, 'phase_averaged ' + lab, lambda args=args: ap.phase_averaged_probabilities_constant_hamiltonian(*args), nm)
C(s, 'phase_averaged spread 0: accepted (documented departure)', lambda: _run(lambda: ap.phase_averaged_probabilities_constant_hamiltonian(H0/E, -H0/E, L, 0.0))[0] is None, 'accepted')
R(s, 'averaged_probabilities_from_eigenbasis non-unitary', lambda: ap.averaged_probabilities_from_eigenbasis(np.ones((3, 3))), 'eigenvectors')
R(s, 'averaged_probabilities_from_eigenbasis NaN', lambda: ap.averaged_probabilities_from_eigenbasis(np.eye(3)*np.nan), 'eigenvectors')
for fn in ('coherence_blocks', 'coherence_report'):
    f = getattr(ap, fn)
    OK(s, 'BASE ' + fn, lambda f=f: f(np.array([1., 2., 3.]), 1.0))
    for lab, args, nm in [('NaN eigenvalue', (np.array([1, 2, np.nan]), 1.0), 'eigenvalues'), ('complex', (np.array([1, 2, 3j]), 1.0), 'eigenvalues'), ('2-D', (np.eye(3), 1.0), 'eigenvalues'),
                          ('phase_scale 0', (np.array([1., 2, 3]), 0.0), 'phase_scale'), ('phase_scale nan', (np.array([1., 2, 3]), np.nan), 'phase_scale'), ('phase_scale True', (np.array([1., 2, 3]), True), 'phase_scale')]:
        R(s, fn + ' ' + lab, lambda f=f, args=args: f(*args), nm)
    R(s, fn + ' decoherence_threshold 0', lambda f=f: f(np.array([1., 2, 3]), 1.0, decoherence_threshold=0.0), 'decoherence_threshold')
Hl = lambda l: np.array([[l - 5, .3], [.3, 5 - l]], dtype=complex)
OK(s, 'BASE averaged_probabilities_adiabatic', lambda: ap.averaged_probabilities_adiabatic(Hl, 0, 10))
for lab, kw, nm in [('n_points 1', dict(n_points=1), 'n_points'), ('l0 >= l1', dict(l0=10, l1=0), 'l0'), ('magnus_exp_order 0', dict(magnus_exp_order=0), 'magnus_exp_order'),
                    ('n_probe 0', dict(n_probe=0), 'n_probe'), ('l1 nan', dict(l1=np.nan), 'l1')]:
    R(s, 'averaged_probabilities_adiabatic ' + lab, lambda kw=kw: ap.averaged_probabilities_adiabatic(Hl, kw.pop('l0', 0), kw.pop('l1', 10), **kw), nm)
OK(s, 'BASE adiabatic_phase_differences', lambda: ap.adiabatic_phase_differences(Hl, 0, 10))
for lab, args, kw, nm in [('n_points 0', (0, 10), dict(n_points=0), 'n_points'), ('n_points 2.5', (0, 10), dict(n_points=2.5), 'n_points'), ('n_points True', (0, 10), dict(n_points=True), 'n_points'),
                          ("n_points 'a'", (0, 10), dict(n_points='a'), 'n_points'), ('l0 > l1', (10, 0), {}, 'l0'), ('l0 nan', (np.nan, 10), {}, 'l0'), ('l1 inf', (0, np.inf), {}, 'l1')]:
    R(s, 'adiabatic_phase_differences ' + lab, lambda args=args, kw=kw: ap.adiabatic_phase_differences(Hl, *args, **kw), nm)
R(s, 'adiabatic_phase_differences non-callable', lambda: ap.adiabatic_phase_differences(3, 0, 10), 'H_func')
OK(s, 'BASE level_crossing_matrix', lambda: ap.level_crossing_matrix(Hl, 0, 10))
for lab, kw, nm in [('threshold 0', dict(threshold=0.0), 'threshold'), ('threshold nan', dict(threshold=np.nan), 'threshold'), ('fd_step_frac 0', dict(fd_step_frac=0.0), 'fd_step_frac'),
                    ('fd_step_frac True', dict(fd_step_frac=True), 'fd_step_frac'), ('n_probe 0', dict(n_probe=0), 'n_probe'), ('magnus_exp_order 11', dict(magnus_exp_order=11), 'magnus_exp_order'),
                    ('magnus_exp_order 2.5', dict(magnus_exp_order=2.5), 'magnus_exp_order'), ("integration_method 'foo'", dict(integration_method='foo'), 'integration_method')]:
    R(s, 'level_crossing_matrix ' + lab, lambda kw=kw: ap.level_crossing_matrix(Hl, 0, 10, **kw), nm)
R(s, 'level_crossing_matrix l1 < l0', lambda: ap.level_crossing_matrix(Hl, 10, 0), 'l0')

# ---------------- §9
s = '9'
OK(s, 'BASE hybrid_propagator', lambda: ad.hybrid_propagator(Hl, 0, 10))
OK(s, 'hybrid l0 == l1 identity', lambda: ad.hybrid_propagator(Hl, 3, 3), lambda v: np.allclose(v[0], np.eye(2)))
for fn in ('hybrid_propagator', 'find_resonance_candidates', 'find_nonadiabatic_windows', 'adiabatic_propagator', 'oscillation_sampling'):
    f = getattr(ad, fn)
    R(s, fn + ' l0 > l1', lambda f=f: f(Hl, 10, 0), 'l0')
    R(s, fn + ' l1 nan', lambda f=f: f(Hl, 0, np.nan), 'l1')
    R(s, fn + ' l1 complex', lambda f=f: f(Hl, 0, 10j), 'l1')
    R(s, fn + ' non-callable', lambda f=f: f(3, 0, 10), 'H_func')
    if fn == 'oscillation_sampling':
        # Revised: a diagnostic returns a quiet empty report rather than raising (its documented
        # design, pinned by test_oscillation_sampling_refuses_quietly_rather_than_breaking_the_call).
        for lab, bad in [('NaN H', lambda l: np.full((2, 2), np.nan)), ('non-square H', lambda l: np.ones((2, 3))),
                         ('non-Hermitian H', lambda l: np.array([[0, 1], [0, 0]], dtype=complex))]:
            OK(s, fn + ' ' + lab + ': quiet empty report', lambda f=f, bad=bad: f(bad, 0, 10), lambda v: v == {})
        continue
    R(s, fn + ' NaN H', lambda f=f: f(lambda l: np.full((2, 2), np.nan), 0, 10), 'H_func')
    R(s, fn + ' non-square H', lambda f=f: f(lambda l: np.ones((2, 3)), 0, 10), 'H_func')
    R(s, fn + ' non-Hermitian H', lambda f=f: f(lambda l: np.array([[0, 1], [0, 0]], dtype=complex), 0, 10), 'Hermitian')
R(s, 'find_hidden_features l0 > l1', lambda: ad.find_hidden_features(lambda l: 1.0, 10, 0), 'l0')
for lab, kw, nm in [('fd_step_frac 0', dict(fd_step_frac=0.0), 'fd_step_frac'), ('fd_step_frac -1', dict(fd_step_frac=-1.0), 'fd_step_frac'), ('fd_step_frac nan', dict(fd_step_frac=np.nan), 'fd_step_frac'),
                    ('max_n_probe 1', dict(max_n_probe=1), 'max_n_probe'), ('n_probe0 2', dict(n_probe0=2), 'n_probe0'), ('n_probe0 0', dict(n_probe0=0), 'n_probe0'), ('n_points0 0', dict(n_points0=0), 'n_points0'),
                    ('rtol 0', dict(rtol=0.0), 'rtol'), ('rtol nan', dict(rtol=np.nan), 'rtol'), ('atol True', dict(atol=True), 'atol'), ('rtol None', dict(rtol=None), 'rtol'),
                    ('threshold0 0', dict(threshold0=0.0), 'threshold0'), ('min_threshold -1', dict(min_threshold=-1.0), 'min_threshold'), ('max_n_points 0', dict(max_n_points=0), 'max_n_points'),
                    ('max_n_points 2.5', dict(max_n_points=2.5), 'max_n_points'), ('max_iters 0', dict(max_iters=0), 'max_iters'), ('magnus_exp_order 2.5', dict(magnus_exp_order=2.5), 'magnus_exp_order'),
                    ('magnus_exp_order 0', dict(magnus_exp_order=0), 'magnus_exp_order'), ('magnus_exp_order 11', dict(magnus_exp_order=11), 'magnus_exp_order')]:
    R(s, 'hybrid_propagator ' + lab, lambda kw=kw: ad.hybrid_propagator(Hl, 0, 10, **kw), nm)
for lab, kw, nm in [('n_probe 0', dict(n_probe=0), 'n_probe'), ('n_probe None', dict(n_probe=None), 'n_probe'), ('n_probe 2.5', dict(n_probe=2.5), 'n_probe'), ('fd_step_frac nan', dict(fd_step_frac=np.nan), 'fd_step_frac'), ('info []', dict(info=[]), 'info')]:
    R(s, 'find_resonance_candidates ' + lab, lambda kw=kw: ad.find_resonance_candidates(Hl, 0, 10, **kw), nm)
R(s, 'find_nonadiabatic_windows threshold nan', lambda: ad.find_nonadiabatic_windows(Hl, 0, 10, threshold=np.nan), 'threshold')
for v in (0, -1, 2.5):
    R(s, 'adiabatic_propagator n_points=%r' % v, lambda v=v: ad.adiabatic_propagator(Hl, 0, 10, n_points=v), 'n_points')
R(s, 'adiabatic_propagator n_points=True', lambda: ad.adiabatic_propagator(Hl, 0, 10, n_points=True), 'n_points')
for v in (0, -1, 2.5, True):
    R(s, 'oscillation_sampling n_probe=%r' % (v,), lambda v=v: ad.oscillation_sampling(Hl, 0, 10, n_probe=v), 'n_probe')
R(s, 'oscillation_sampling baselines NaN', lambda: ad.oscillation_sampling(Hl, 0, 10, baselines=np.array([np.nan])), 'baselines')
R(s, 'oscillation_sampling baselines negative', lambda: ad.oscillation_sampling(Hl, 0, 10, baselines=np.array([-1.0])), 'baselines')
OK(s, 'BASE find_hidden_features', lambda: ad.find_hidden_features(lambda l: 1.0 + 0*np.asarray(l), 0, 10, n_ref=64))
R(s, 'find_hidden_features non-callable', lambda: ad.find_hidden_features(3, 0, 10), 'profile')
# Revised: documented quiet refusal ('hidden=False there means "not measured"').
OK(s, 'find_hidden_features NaN profile: quiet report', lambda: ad.find_hidden_features(lambda l: np.nan*np.asarray(l), 0, 10, n_ref=64), lambda v: v.get('hidden') is False)
for lab, kw, nm in [('n_ref 0', dict(n_ref=0), 'n_ref'), ('n_ref -1', dict(n_ref=-1), 'n_ref'), ('n_ref True', dict(n_ref=True), 'n_ref'), ('n_ref 2.5', dict(n_ref=2.5), 'n_ref'), ('n_sub 0', dict(n_sub=0), 'n_sub')]:
    R(s, 'find_hidden_features ' + lab, lambda kw=kw: ad.find_hidden_features(lambda l: 1.0 + 0*np.asarray(l), 0, 10, **kw), nm)

# ---------------- §10
s = '10'
for v in (2.5, True, '4'):
    R(s, 'multislab order=%r' % (v,), lambda v=v: mm.magnus_expansion_multislab(A, [[0, L]], order=v), 'order')
for v in (0, -1, True, 2.5, '3'):
    R(s, 'gl_nodes(%r)' % (v,), lambda v=v: mm.gl_nodes(v), 'order')
R(s, 'symmetric_over=(0,)', lambda: mm.magnus_expansion_multislab(A, [[0, L]], symmetric_over=(0,)), 'symmetric_over')
R(s, "A_eval_mode='zzz'", lambda: mm.magnus_expansion_multislab(A, [[0, L]], A_eval_mode='zzz'), 'A_eval_mode')
for lab, kw in [('trapezoid n_tpts 0', dict(integration_method='trapezoid', n_tpts_per_slab=0)), ('trapezoid n_tpts 1', dict(integration_method='trapezoid', n_tpts_per_slab=1)),
                ('simpson n_tpts 1', dict(integration_method='simpson', n_tpts_per_slab=1)), ('n_tpts 5.5', dict(integration_method='trapezoid', n_tpts_per_slab=5.5))]:
    R(s, 'multislab ' + lab, lambda kw=kw: mm.magnus_expansion_multislab(A, [[0, L]], **kw), 'n_tpts_per_slab')
for lab, a in [('np.eye(2)', np.eye(2)), ('scalar return', lambda t: 1.0), ('2x3 return', lambda t: np.ones((2, 3)))]:
    R(s, 'multislab A ' + lab, lambda a=a: mm.magnus_expansion_multislab(a, [[0, L]]), 'A')
R(s, 'ordered_product empty stack', lambda: mm.ordered_product(np.zeros((0, 2, 2))), 'U')
R(s, 'commutator mismatched shapes', lambda: mm.commutator(np.eye(2), np.eye(3)), 'X')

# ---------------- §11
s = '11'
for lab, v in [('0', 0.0), ('negative', -1e9), ('nan', np.nan), ('inf', np.inf), ('complex', 1e9 + 1j), ('True', True)]:
    R(s, 'builder energy ' + lab, lambda v=v: h.hamiltonian_3nu_vacuum(v, **p), 'energy')
for a in ('D21', 'dCP'):
    R(s, 'builder %s nan' % a, lambda a=a: h.hamiltonian_3nu_vacuum_energy_independent(**{**p, a: np.nan}), a)
R(s, 'builder eps_ee complex', lambda: h.hamiltonian_3nu_nsi(1e-13, 0.1j, 0, 0, 0, 0, 0), 'eps_ee')
R(s, 'builder Lambda nan', lambda: h.hamiltonian_3nu_liv_energy_independent(0.3, 0.4, 0.2, 0.1, 1e-23, 2e-23, 3e-23, np.nan, 0), 'Lambda')
for v in (2, 'no', np.array([True])):
    R(s, 'builder nubar=%r' % (v if not isinstance(v, np.ndarray) else 'array'), lambda v=v: h.hamiltonian_3nu_vacuum_energy_independent(**p, nubar=v), 'nubar')
for lab, kw in [('rad 2.0', dict(s12=2.0, angles='rad')), ('deg 120', dict(s12=120.0, s23=40., s13=8.6, dCP=200., angles='deg')), ('rad -2.0', dict(s12=-2.0, angles='rad'))]:
    R(s, 'angle beyond 90 degrees: ' + lab, lambda kw=kw: o.osc_prob_3nu_vacuum(E, L, **{**dict(s12=0.6, s23=0.7, s13=0.15, dCP=1.0, D21=7.5e-5, D31=2.5e-3), **kw}), 's12')
OK(s, 'angle -1.0 rad accepted', lambda: o.osc_prob_3nu_vacuum(E, L, s12=-1.0, s23=0.7, s13=0.15, dCP=1.0, D21=7.5e-5, D31=2.5e-3, angles='rad'))
OK(s, 'angle pi/2 rad accepted', lambda: o.osc_prob_3nu_vacuum(E, L, s12=np.pi/2, s23=0.7, s13=0.15, dCP=1.0, D21=7.5e-5, D31=2.5e-3, angles='rad'))
W(s, 'angle warning: 4 flavors, sines in active slots, 5-degree sterile angles', lambda: o.osc_prob_4nu_vacuum(E, L, s12=0.55, s23=0.75, s13=0.15, dCP=1.0, s14=5., s24=5., s34=5., d14=0., d24=0., D21=7.5e-5, D31=2.5e-3, D41=1., angles='deg'), 'MixingAngleConventionWarning')
def _once():
    cls, msg, ws = _run(lambda: h.hamiltonian_3nu_liv(1e9, 0.5, 0.5, 0.5, 0.0, 1e-23, 2e-23, 3e-23, 1e9, 0, angles='deg'))
    return sum('MixingAngleConventionWarning' in x for x in ws) == 1
C(s, 'angle warning once per call (LIV builder)', _once, 'exactly one')
W(s, 'angle warning text: false positive and how to silence', lambda: o.osc_prob_2nu_vacuum(E, L, sth=0.5, Dm2=7.5e-5, angles='deg'), 'filterwarnings')
OK(s, 'BASE exp_density_profile', lambda: ma.exp_density_profile(3.0, 1.0))
for lab, args, nm in [('central -3', (-3.0, 1.0), 'density_matter_central'), ('l_scale 0', (3.0, 0.0), 'l_scale'), ('l_scale -1', (3.0, -1.0), 'l_scale'), ('l_scale nan', (3.0, np.nan), 'l_scale'),
                      ('l_scale inf', (3.0, np.inf), 'l_scale'), ('l_scale True', (3.0, True), 'l_scale'), ("central 'a'", ('a', 1.0), 'density_matter_central')]:
    R(s, 'exp_density_profile ' + lab, lambda args=args: ma.exp_density_profile(*args), nm)
for lab, kw, nm in [('electron_fraction 1.5', dict(electron_fraction=1.5), 'electron_fraction'), ('ratio -1', dict(ratio_number_neutrons_to_protons=-1.0), 'ratio_number_neutrons_to_protons'),
                    ("g_per_cm3 'no'", dict(density_matter_is_in_g_per_cm3='no'), 'density_matter_is_in_g_per_cm3'), ('number_of_electrons 2', dict(density_is_of_number_of_electrons=2), 'density_is_of_number_of_electrons')]:
    R(s, 'vcc_func_from_rho_func ' + lab, lambda kw=kw: ma.vcc_func_from_rho_func(lambda l: 3.0, **kw), nm)
for k in ('core', 'mantle', 'crust', 'ocean'):
    for v in (-0.1, 1.5, np.nan, np.inf, True):
        R(s, 'electron_fraction_func_prem %s=%r' % (k, v), lambda k=k, v=v: ea.electron_fraction_func_prem(1000.0, **{'electron_fraction_' + k: v}), 'electron_fraction_' + k)
OK(s, 'BASE matter_potential_projector', lambda: ma.matter_potential_projector(3))
for v in (1, 0, 2.5, True, 'a'):
    R(s, 'matter_potential_projector num_flavors=%r' % (v,), lambda v=v: ma.matter_potential_projector(v), 'num_flavors')
for v in (-1.0, np.nan, np.inf, True):
    R(s, 'matter_potential_projector ratio=%r' % (v,), lambda v=v: ma.matter_potential_projector(3, v), 'ratio_number_neutrons_to_protons')
for v in (0.0, -0.1, 1.5, np.nan, np.inf, True):
    R(s, 'neutron_to_proton_ratio_from_electron_fraction %r' % (v,), lambda v=v: ea.neutron_to_proton_ratio_from_electron_fraction(v), 'electron_fraction')
dmf = lambda l: 3.0
OK(s, 'BASE num_density_e_func', lambda: ma.num_density_e_func(0.0, dmf))
for v in (1.5, -0.1, np.nan, np.inf, 1j, True):
    R(s, 'num_density_e_func electron_fraction=%r' % (v,), lambda v=v: ma.num_density_e_func(0.0, dmf, electron_fraction=v), 'electron_fraction')
R(s, 'num_density_e_func ratio -1', lambda: ma.num_density_e_func(0.0, dmf, ratio_number_neutrons_to_protons=-1.0), 'ratio_number_neutrons_to_protons')
R(s, 'num_density_e_func non-callable density', lambda: ma.num_density_e_func(0.0, 3), 'density_matter_func')
for fn, modname in [('density_matter_func_const', ma), ('density_matter_func_exp', ma), ('density_matter_func_prem', ea)]:
    f = getattr(modname, fn, None)
    C(s, 'hot-path note in %s docstring' % fn, lambda f=f: f is not None and 'hot path' in (f.__doc__ or ''), 'docstring note')
C(s, 'hot-path note in hamiltonian_3nu_matter docstring', lambda: 'hot path' in (h.hamiltonian_3nu_matter.__doc__ or ''), 'docstring note')
C(s, 'hot-path note in a *_td builder docstring', lambda: 'hot path' in (h.hamiltonian_3nu_matter_td.__doc__ or ''), 'docstring note')
OK(s, 'BASE std vacuum', lambda: std.osc_prob_2nu_vacuum_std(0.5, 7.5e-5, E, L))
for arg, i in [('energy', 2), ('L', 3), ('Dm2', 1)]:
    for lab, v in [('negative', -1.0), ('nan', np.nan), ('inf', np.inf), ('complex', 1j), ('True', True)]:
        if arg == 'Dm2' and lab == 'negative':
            continue
        def _std(i=i, v=v, fn='vacuum'):
            a = [0.5, 7.5e-5, E, L]; a[i] = v
            return std.osc_prob_2nu_vacuum_std(*a)
        R(s, 'std vacuum %s %s' % (arg, lab), _std, arg)
OK(s, 'std matter Dm2=0', lambda: std.osc_prob_2nu_matter_std(0.5, 0.0, 1e-13, E, L))

# ---------------- §12
s = '12'
X = np.linspace(0.5, 5, 50); Y = np.full(50, 0.1); CZ = np.linspace(-1, -0.1, 4); LG = np.linspace(9, 10, 4)
OK(s, 'BASE plot_probability_vs_energy', lambda: mp.plot_probability_vs_energy(X, [Y], nu_i=1, nu_f=0))
R(s, 'ax= names the routine called', lambda: mp.plot_probability_vs_energy(X, [Y], nu_i=1, nu_f=0, ax=None), 'plot_probability_vs_energy')
R(s, 'energy_unit on vs_baseline names the routine', lambda: mp.plot_probability_vs_baseline(X, [Y], energy_unit='MeV'), 'plot_probability_vs_baseline')
for k in ('nrows', 'ncols', 'figsize'):
    R(s, 'subplots_kw %s refused' % k, lambda k=k: mp.plot_probability_vs_energy(X, [Y], nu_i=1, nu_f=0, subplots_kw={k: 2 if k != 'figsize' else (3, 3)}), k)
for v in ('foo', 'parsec'):
    R(s, 'energy_unit %r' % v, lambda v=v: mp.plot_probability_vs_energy(X, [Y], nu_i=1, nu_f=0, energy_unit=v), 'energy_unit')
R(s, 'x_unit int', lambda: mp.plot_probability_with_profile(np.linspace(1, 1300, 50), [np.full(50, 0.4)], [[dict(y=np.linspace(0, 1, 50))]], x_unit='km', xscale='linear'), 'x_unit')
W(s, 'energies in eV with GeV unit warn', lambda: mp.plot_probability_vs_energy(X*1e9, [Y], nu_i=1, nu_f=0), 'energy')
R(s, 'nu_f outside num_flavors', lambda: mp.plot_probability_vs_energy(X, [Y], nu_i=1, nu_f=gd.NUS2, num_flavors=3), 'nu_f')
R(s, 'curve wrong length', lambda: mp.plot_curves(X, [np.ones(5)]), 'curves[0]')
R(s, 'curve 2-D', lambda: mp.plot_curves(X, [np.ones((50, 2))]), 'curves[0]')
R(s, 'x <= 0 on log axis', lambda: mp.plot_probability_vs_energy(np.array([0., 1., 2.]), [np.array([.1, .2, .3])], nu_i=1, nu_f=0), 'xscale')
R(s, 'all y <= 0 on log y axis', lambda: mp.plot_curves(np.linspace(1, 5, 50), [np.zeros(50)], yscale='log'), 'yscale')
R(s, 'all y <= 0 on log y axis (stacked)', lambda: mp.plot_curves_stacked(np.linspace(1, 5, 50), [[np.zeros(50)]], yscale='log'), 'yscale')
R(s, 'oscillogram 1x1 grid', lambda: mp.plot_oscillogram(np.array([-1.]), np.array([1.]), np.array([[0.5]]), nu_i=1, nu_f=1), 'costhz')
R(s, 'oscillogram transposed', lambda: mp.plot_oscillogram(CZ, np.linspace(9, 10, 5), np.full((4, 5), 0.5), nu_i=1, nu_f=1), 'probability')
R(s, 'probabilities above 1 in vs_energy', lambda: mp.plot_probability_vs_energy(X, [np.full(50, 1.5)], nu_i=1, nu_f=0), '[0, 1]')
R(s, 'NaN in probabilities', lambda: mp.plot_probability_vs_energy(X, [np.full(50, np.nan)], nu_i=1, nu_f=0), 'finite')
R(s, 'NaN abscissa', lambda: mp.plot_curves(np.r_[X[:-1], np.nan], [Y]), 'x')
R(s, 'string data', lambda: mp.plot_curves(X, [['a']*50]), 'curves[0]')
C(s, 'oscillogram accepts a difference map (documented departure)', lambda: _run(lambda: mp.plot_oscillogram(CZ, LG, np.full((4, 4), -0.2), nu_i=1, nu_f=0))[0] is None, 'accepted')
R(s, 'compute mode: dcp NaN', lambda: mp.plot_biprobability(dcp=[np.nan], energy=2.0, nu_i=1, nu_f=0, num_flavors=3, configurations=[dict(environment='vacuum', L=1300.0)]), 'dcp')
R(s, 'compute mode: dcp empty', lambda: mp.plot_biprobability(dcp=[], energy=2.0, nu_i=1, nu_f=0, num_flavors=3, configurations=[dict(environment='vacuum', L=1300.0)]), 'dcp')
R(s, 'compute mode: log10_energy [30, 31]', lambda: mp.plot_oscillogram(CZ, np.array([30., 31.]), nu_i=1, nu_f=0, num_flavors=3), 'log10_energy')
R(s, 'compute mode: bad physics named by the plotting routine', lambda: mp.plot_oscillogram(CZ, LG, nu_i=1, nu_f=0, num_flavors=3, osc_params=dict(p, dcp=1.0)), 'plot_oscillogram')
for lab, kw, nm in [('residual_height 1.5', dict(residual=Y, residual_height=1.5), 'residual_height'), ('residual_height 0', dict(residual=Y, residual_height=0.0), 'residual_height'),
                    ('figsize negative', dict(figsize=(-1, 3)), 'figsize'), ('xmajor 0', dict(xmajor=0.0), 'xmajor'), ('xmajor nan', dict(xmajor=np.nan), 'xmajor'),
                    ('xlim 3-tuple', dict(xlim=(1, 2, 3)), 'xlim'), ('xlim nan', dict(xlim=(np.nan, 2)), 'xlim'), ("xscale 'foo'", dict(xscale='foo'), 'xscale'),
                    ("legend 'yes'", dict(legend='yes'), 'legend'), ("grid 2", dict(grid=2), 'grid'), ("tight_layout 'no'", dict(tight_layout='no'), 'tight_layout')]:
    R(s, 'layout ' + lab, lambda kw=kw: mp.plot_curves(X, [Y], **kw), nm)
R(s, 'layout levels 0', lambda: mp.plot_oscillogram(CZ, LG, np.full((4, 4), 0.5), nu_i=1, nu_f=0, levels=0), 'levels')
R(s, 'layout legend_panel 1.5', lambda: mp.plot_curves_stacked(X, [[Y], [Y]], legend_panel=1.5), 'legend_panel')
R(s, 'return_probability 2', lambda: mp.plot_oscillogram(CZ, LG, np.full((4, 4), 0.5), nu_i=1, nu_f=0, return_probability=2), 'return_probability')
R(s, 'fontsize negative', lambda: mp.plot_oscillogram(CZ, LG, np.full((4, 4), 0.5), nu_i=1, nu_f=0, cbar_fontsize=-1), 'cbar_fontsize')

# ---------------- §13
s = '13'
base = ('--environment', 'vacuum', '--energy', '1', '--baseline', '1000')
OK(s, 'BASE CLI vacuum', lambda: cli(*base))
for args, nm in [(('--energy', 'inf'), '--energy'), (('--energy', 'nan'), '--energy'), (('--baseline', 'nan'), '--baseline'), (('--baseline', 'inf'), '--baseline'), (('--dcp', 'nan'), '--dcp'),
                 (('--dcp', 'inf'), '--dcp'), (('--dm21', 'nan'), '--dm21'), (('--dm21', 'inf'), '--dm21'), (('--precision', '-1'), '--precision'), (('--precision', '100'), '--precision')]:
    R(s, 'CLI ' + ' '.join(args), lambda args=args: cli(*[a for a in base if a not in ('--energy', '1', '--baseline', '1000')] + (['--energy', '1'] if args[0] != '--energy' else []) + (['--baseline', '1000'] if args[0] != '--baseline' else []) + list(args)), nm)
mat = ('--environment', 'matter', '--energy', '1', '--baseline', '1000', '--rho', '3')
OK(s, 'BASE CLI matter', lambda: cli(*mat))
for args, nm in [(('--electron-fraction', '1.5'), '--electron-fraction'), (('--electron-fraction', 'nan'), '--electron-fraction'), (('--ratio-n-to-p', 'nan'), '--ratio-n-to-p')]:
    R(s, 'CLI ' + ' '.join(args), lambda args=args: cli(*mat, *args), nm)
OK(s, 'BASE CLI nsi', lambda: cli(*mat, '--scenario', 'nsi', '--eps-ee', '0.1'))
R(s, 'CLI --eps-ee nan', lambda: cli(*mat, '--scenario', 'nsi', '--eps-ee', 'nan'), '--eps-ee')
R(s, 'CLI --eps-ee inf', lambda: cli(*mat, '--scenario', 'nsi', '--eps-ee', 'inf'), '--eps-ee')
liv = ('--environment', 'vacuum', '--scenario', 'liv', '--energy', '1', '--baseline', '1000', '--b1', '1e-23', '--liv-lambda', '1e9')
OK(s, 'BASE CLI liv', lambda: cli(*liv))
R(s, 'CLI --n-liv -1', lambda: cli(*liv, '--n-liv', '-1'), '--n-liv')
R(s, 'CLI --dxi13 with --dxicp', lambda: cli(*liv, '--dxi13', '1', '--dxicp', '2'), '--dxi13')
OK(s, 'BASE CLI earth', lambda: cli('--environment', 'earth', '--energy', '1', '--costhz', '-0.5', '--baseline', '100'))
R(s, 'CLI --costhz -1.5 --baseline 100', lambda: cli('--environment', 'earth', '--energy', '1', '--costhz', '-1.5', '--baseline', '100'), '--costhz')
R(s, 'CLI --costhz 2', lambda: cli('--environment', 'earth', '--energy', '1', '--costhz', '2', '--baseline', '100'), '--costhz')
for args, nm in [(('--rho', '3'), '--rho')]:
    R(s, 'CLI inapplicable --rho with vacuum', lambda: cli(*base, '--rho', '3'), '--rho')
R(s, 'CLI inapplicable --eps-ee with scenario std', lambda: cli(*mat, '--eps-ee', '0.1'), '--eps-ee')
R(s, 'CLI inapplicable --s14 with 3 flavors', lambda: cli(*base, '--s14', '0.1'), '--s14')
R(s, 'CLI inapplicable --rho-central with constant density', lambda: cli(*mat, '--rho-central', '3'), '--rho-central')
R(s, 'CLI inapplicable --costhz with matter', lambda: cli(*mat, '--costhz', '-0.5'), '--costhz')
R(s, 'CLI --rho -3 names --rho', lambda: cli('--environment', 'matter', '--energy', '1', '--baseline', '1000', '--rho', '-3'), '--rho')

# ---------------- §14
s = '14'
OK(s, 'BASE cross_check_strategies', lambda: o.cross_check_strategies(o.osc_prob_3nu_matter_exp_density, E, L, 0.0, 3.0, 300*KM, density_matter_is_in_g_per_cm3=True, nu_i=0, nu_f=0))
for lab, v in [("unknown 'zzz'", ['zzz']), ("bare string", 'hybrid'), ('int', 3), ('[]', [])]:
    R(s, 'engines ' + lab, lambda v=v: o.cross_check_strategies(o.osc_prob_3nu_matter_exp_density, E, L, 0.0, 3.0, 300*KM, density_matter_is_in_g_per_cm3=True, nu_i=0, nu_f=0, engines=v), 'engine')

res = [r for r in out if isinstance(r, dict)]
if len(sys.argv) > 2:
    json.dump(res, open(sys.argv[2], 'w'), indent=1)
bad_base = [r for r in res if r['id'].startswith('BASE') and not r['ok']]
print('cases:', len(res), ' pass:', sum(r['ok'] for r in res), ' fail:', sum(not r['ok'] for r in res), ' broken BASE:', len(bad_base))
for r in res:
    if not r['ok']:
        print('FAIL §%s | %s | want %s | got %s' % (r['sec'], r['id'], r['want'], r['got']))
