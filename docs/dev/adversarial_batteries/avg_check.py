# -*- coding: utf-8 -*-
"""Does the 5 MeV silent miss survive PHASE AVERAGING?

For solar neutrinos the observable is the averaged survival probability: the vacuum phase over
1 AU is ~1e10 cycles, the 8B production region is extended, and detector resolution washes the
rest.  An error in the instantaneous probability at one baseline can be almost entirely a phase
error, which averages away and never reaches a user.

So the P1 criterion used all session -- instantaneous P within 1e-3 at a single baseline -- may
be the wrong observable for exactly the family that produced the headline.  This measures both:
the instantaneous error at the nominal baseline, and the error in the average over a window of
several oscillation lengths.
"""
import warnings
import numpy as np
import harness as H
import magnus.oscprob as op
import magnus.globaldefs as gd
import physical_profiles as pp

def main():
    fam = {k: next(f for f in pp.families() if f['label'] == 'BS05(AGS,OP) ' + k)
           for k in ('cubic', 'linear')}
    d, E = 2, 5.0e6
    params = H.params_for(d)
    p = gd.OSC_PARAMS_PREDEFINED['OSC_PARAMS_DEFAULT']
    L_osc = 4.0*np.pi*E/p['D21']
    print('# oscillation length at %.0f MeV: %.3e eV^-1' % (E/1e6, L_osc))
    for kind, f in fam.items():
        l0, l1 = f['l0'], f['l1']
        print('# trajectory = %.3e = %.0f oscillation lengths' % (l1 - l0, (l1 - l0)/L_osc))
        # A window of ~6 oscillation lengths ending at the nominal baseline.
        Ls = np.linspace(l1 - 6.0*L_osc, l1, 121)
        Hf = H.H_factory(d, params, H.vcc_of(f['ne']), E)
        Pref = np.array([H.P_of(U) for U in H.exact_U_many(Hf, l0, Ls, d)])
        Pget = []
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            for L in Ls:
                Pget.append(np.asarray(op.osc_prob_matter_std_potential(
                    d, f['ne'], E, float(L), params, L0=l0,
                    density_is_of_number_of_electrons=True)))
        Pget = np.array(Pget)
        inst = float(np.max(np.abs(Pget[-1] - Pref[-1])))
        inst_worst = float(np.max(np.abs(Pget - Pref)))
        avg_err = float(np.max(np.abs(Pget.mean(axis=0) - Pref.mean(axis=0))))
        print('\n%s:' % kind)
        print('  instantaneous error at the nominal baseline : %.3e' % inst)
        print('  worst instantaneous error over the window   : %.3e' % inst_worst)
        print('  error in the AVERAGE over the window        : %.3e  <-- the solar observable'
              % avg_err)
        print('  averaging reduces the error by              : %.0fx' % (inst/max(avg_err, 1e-30)))
        print('  P_ee averaged: package %.6f   oracle %.6f'
              % (Pget.mean(axis=0)[0, 0], Pref.mean(axis=0)[0, 0]))
    return 0

if __name__ == '__main__':
    raise SystemExit(main())
