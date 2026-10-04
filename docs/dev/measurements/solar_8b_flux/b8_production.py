r"""The 8B production profile of the B16-GS98 solar model, from its structure.

The B16 release published the radial distribution of neutrino production on a page that
is no longer online (Vinyoles et al. 2017, footnote 1).  This script rebuilds the 8B one
from the model's structure: the rate per volume of 7Be(p,gamma)8B, with 7Be in
equilibrium between 3He(4He,gamma)7Be and its destruction by electron and proton capture.
Non-resonant rates come from the Gamow-peak formula with the SFII S-factors (Adelberger
et al. 2011, Rev. Mod. Phys. 83, 195), Salpeter weak screening, and the 7Be
electron-capture rate of that review.  Run from this directory:

    python b8_production.py

It first checks the method on BS2005-AGS,OP, whose 8B distribution Bahcall published,
and then writes b16_gs98_b8_production.csv.
"""
import pathlib

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
BATTERIES = HERE.parent.parent/'adversarial_batteries'


def nacsv(T9, Z1, Z2, A1, A2, S_MeVb):
    """N_A <sigma v> in cm^3 mol^-1 s^-1, non-resonant, S in MeV b."""
    mu = A1*A2/(A1 + A2)
    tau = 4.2487*(Z1**2*Z2**2*mu/T9)**(1/3)
    return 7.83e9*(Z1*Z2/(mu*T9**2))**(1/3)*S_MeVb*np.exp(-tau)*(1 + 5/(12*tau))


def screen(Z1, Z2, rho, T, X, Y):
    """Salpeter weak screening, for a hydrogen and helium plasma."""
    zeta = np.sqrt(2.0*X + 1.5*Y)
    return np.exp(0.188*Z1*Z2*zeta*np.sqrt(rho)*(T/1e6)**-1.5)


def b8_profile(r, T, rho, X, Y, X3):
    """8B production per unit radius, normalized to 1 over r."""
    T9, T6 = T/1e9, T/1e6
    n_p, n_4, n_3 = rho*X/1.0078, rho*Y/4.0026, rho*X3/3.0160          # mol cm^-3
    s34 = nacsv(T9, 2, 2, 3, 4, 0.56e-3)*screen(2, 2, rho, T, X, Y)
    s17 = nacsv(T9, 4, 1, 7, 1, 20.8e-6)*screen(4, 1, rho, T, X, Y)
    lam_e = 5.60e-9*rho*(1 + X)/2/np.sqrt(T6)*(1 + 0.004*(T6 - 16))    # s^-1
    n_7 = n_3*n_4*s34/(lam_e + n_p*s17)
    dP = 4*np.pi*r**2*n_p*n_7*s17
    return dP/np.trapezoid(dP, r)


def quantiles(r, p):
    c = np.cumsum(p*np.gradient(r))
    c /= c[-1]
    return r[np.argmax(p)], *(r[np.searchsorted(c, q)] for q in (0.1, 0.5, 0.9))


def _rows(path, n):
    out = []
    for line in open(path):
        cols = line.split()
        if len(cols) == n:
            try:
                out.append([float(c) for c in cols])
            except ValueError:
                pass
    return np.array(out)


def check_bs05():
    """The method against Bahcall's published 8B distribution for BS2005-AGS,OP."""
    s = _rows(BATTERIES/'bs05_agsop.dat', 12)
    p = b8_profile(s[:, 1], s[:, 2], s[:, 3], s[:, 6], s[:, 7], s[:, 8])
    f = _rows(BATTERIES/'bs2005agsopflux.csv', 13)
    rb, pb = f[:, 0], f[:, 6]/np.trapezoid(f[:, 6], f[:, 0])
    pi = np.interp(rb, s[:, 1], p)
    print('BS05 published: peak %.4f, 10%% %.4f, median %.4f, 90%% %.4f' % quantiles(rb, pb))
    print('BS05 computed:  peak %.4f, 10%% %.4f, median %.4f, 90%% %.4f' % quantiles(rb, pi))
    print('largest difference: %.3f of the peak' % (np.max(np.abs(pi - pb))/pb.max()))


def main():
    check_bs05()
    d = np.loadtxt(HERE/'b16_gs98_struct_columns.csv', delimiter=',', comments='#')
    p = b8_profile(*d.T)
    print('B16-GS98:       peak %.4f, 10%% %.4f, median %.4f, 90%% %.4f' % quantiles(d[:, 0], p))
    np.savetxt(HERE/'b16_gs98_b8_production.csv', np.c_[d[:, 0], p], delimiter=',', fmt='%.6e',
               header='8B production per unit radius in B16-GS98, normalized to 1 over r in R_sun;\n'
                      'computed by b8_production.py from b16_gs98_struct_columns.csv\n'
                      'r_over_r_sun,dP_dr')


if __name__ == '__main__':
    main()
