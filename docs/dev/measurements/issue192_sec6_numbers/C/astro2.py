import warnings
with warnings.catch_warnings(record=True) as C:
    warnings.simplefilter('always')
    exec(open('astro_listing.py').read().split('for k, v in f_earth.items()')[0])
print('listing warnings:', sorted({(w.category.__name__) for w in C}))
for w in C: print('   ', w.category.__name__, str(w.message)[:150])
import magnus.avgprob as avgprob
E1 = 100.0*gd.UNIT_TEV
for dm2 in (1.0e-13, 1.0e-19):
    H = ham.hamiltonian_pseudo_dirac_vacuum(E1, U, m2, {1: dm2})
    lam = np.linalg.eigvalsh(H)
    b = avgprob.coherence_blocks(lam, L_src)
    print(b)
for dm2 in (1e-19, 1e-16, 1e-13):
    print('dm2 %.0e phase %.3g rad' % (dm2, dm2*L_src/(2*E1)))
# scan for block transitions
for dm2 in np.logspace(-20,-14,13):
    H = ham.hamiltonian_pseudo_dirac_vacuum(E1, U, m2, {1: dm2}); print('%.1e'%dm2, avgprob.coherence_blocks(np.linalg.eigvalsh(H), L_src))
# 1e8 km
L8 = 1e8*gd.UNIT_KM
print('E where 21 phase D21 L/2E = 2pi at 1e8 km: %.2f TeV' % (osc['D21']*L8/(4*np.pi)/gd.UNIT_TEV))
with warnings.catch_warnings(record=True) as C2:
    warnings.simplefilter('always')
    P8 = oscprob.osc_prob_3nu_vacuum(E, L8, average=True, **osc)
print('1e8 km warnings', sorted({w.category.__name__ for w in C2}))
for w in C2[:3]: print('   ', str(w.message)[:300])
# 3+1 ternary marks
for s2 in (0.1,0.2,0.3):
    P = oscprob.osc_prob_4nu_vacuum(E[:1], L_src, average=True, s14=np.sqrt(s2), s24=np.sqrt(s2), D41=1.0, **osc)
    print('3+1 s2=%.1f' % s2, at_earth(np.asarray(P)[0]).round(4))
# LIV scan over E for ternary: energies 32, 100, 316
for e in (31.6,100.,316.):
    r = 2*b3*(e*gd.UNIT_TEV)**2/osc['D31']; print('ratio new/vac at %g TeV: %.3f' % (e, r))
