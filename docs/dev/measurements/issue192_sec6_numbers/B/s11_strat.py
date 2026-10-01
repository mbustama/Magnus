import numpy as np, warnings
import magnus.oscprob as oscprob
import magnus.earth as earth
import magnus.globaldefs as gd
L = earth.distance_traveled_inside_earth(-0.9)*gd.UNIT_KM
info = {}
P = oscprob.osc_prob_3nu_earth(10.0*gd.UNIT_GEV, costhz=-0.9, L=L, nu_i=gd.NUMU, nu_f=gd.NUMU, strategy_info=info)
print(sorted(info), len(info))
print(info['engine'], info['family'], info['certified'], info['declined'], info['hidden_feature'])
s = info['sampling']
print(s)
print(s['oscillation_length']/gd.UNIT_KM, s['cycles_over_trajectory'], s['nyquist_points'])
kw = dict(costhz=-0.9, L=L, nu_i=gd.NUMU, nu_f=gd.NUMU, strategy_info=info)
E = np.linspace(5.0, 15.0, 3)*gd.UNIT_GEV
oscprob.osc_prob_3nu_earth(E, **kw); print(info['engine'])
oscprob.osc_prob_3nu_earth(10.0*gd.UNIT_GEV, **kw, average=True); print(info['engine'])
oscprob.osc_prob_3nu_matter_constant_density(10.0*gd.UNIT_GEV, L, 4.0, density_matter_is_in_g_per_cm3=True, nu_i=gd.NUMU, nu_f=gd.NUMU, strategy_info=info); print(info['engine'])
for ls in (0.05, 50.0):
    info = {}
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter('always')
        P = oscprob.osc_prob_2nu_matter_exp_density(5.0*gd.UNIT_MEV, 12742.0*gd.UNIT_KM, 0.0, 3.0, ls*gd.UNIT_KM, sth=0.5, Dm2=2.5e-3, density_matter_is_in_g_per_cm3=True, nu_i=0, nu_f=1, strategy_info=info)
    print(ls, 'declined', info['declined'], 'engine', info['engine'], 'certified', info['certified'], 'nwarn', len(w), [x.category.__name__ for x in w])
for f in ('osc_prob_2nu_vacuum_std','osc_prob_3nu_vacuum_std','osc_prob_2nu_matter_std','osc_prob'):
    import inspect
    fn = getattr(oscprob, f, None)
    print(f, fn is not None and 'strategy_info' in inspect.signature(fn).parameters)
