"""Log the agreement at the clamped step for the #94 case and both #122 cases, and score the returns."""
import json, warnings, numpy as np
import magnus.oscprob as op, magnus.globaldefs as gd
S = '/tmp/claude-1000/-home-mbustamante-Research-magnus/93bacad5-e69e-4b42-9e80-36335d33ac78/scratchpad'
LOG = S + '/clamp122.jsonl'
def mark(name):
    open(LOG, 'a').write(json.dumps(dict(case=name)) + '\n')
osc = gd.load_nufit_params('NuFIT 6.1')
# the #94 case
rho = lambda l: 3e3*np.exp(-np.asarray(l, dtype=float)/(100.0*gd.UNIT_KM))
mark('94: 3nu cap 115')
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    op.osc_prob_matter_std_potential(3, rho, np.linspace(0.1, 0.4, 40)*gd.UNIT_GEV, 250.0*gd.UNIT_KM, osc, L0=0.0,
        density_matter_is_in_g_per_cm3=True, rtol=1e-6, atol=1e-6, strategy='magnus', max_n_slabs=115)
# the #122 cases: Listing 1's 5nu curve, 26 energies
ster = dict(osc, s14=np.sqrt(0.1), s24=np.sqrt(0.1), s15=np.sqrt(0.06), s25=np.sqrt(0.06), D41=1.0, D51=1.7)
Eg = np.logspace(np.log10(2.0), np.log10(20.0), 26)
refs = {}
for line in open(S + '/../../../../../home/mbustamante/Research/magnus/resources/handover_auto_tight/hp/b_fix.jsonl'.replace(S + '/../../../../..', '')):
    r = json.loads(line)
    if r['name'] == 'L1-5nu' and 'ref13' in r:
        refs[round(r['E_gev'], 12)] = (r['ref13'], abs(r['ref12'] - r['ref13']))
ref = np.array([refs[round(e, 12)][0] for e in Eg]); spread = np.array([refs[round(e, 12)][1] for e in Eg])
out = {}
for floor in (13, 19):
    mark('122: 5nu floor %d' % floor)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        P = np.asarray(op.osc_prob_5nu_matter_exp_density(Eg*gd.UNIT_GEV, L=25.0*gd.UNIT_KM, L0=0.0, rho_central=3.0e3,
            l_scale=10.0*gd.UNIT_KM, density_matter_is_in_g_per_cm3=True, nu_i=gd.NUE, nu_f=gd.NUE, strategy='magnus',
            magnus_exp_order=4, rtol=5e-13, atol=5e-15, n_slabs=floor, **ster), float)
    out[floor] = dict(P=P.tolist())
json.dump(dict(out=out, ref=ref.tolist(), spread=spread.tolist()), open(S + '/clamp_cases_P.json', 'w'))
print('max DOP853 spread %.1e' % spread.max())
