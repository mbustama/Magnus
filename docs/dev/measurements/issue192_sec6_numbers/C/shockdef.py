import numpy as np
import magnus.globaldefs as gd
KM = gd.UNIT_KM
MEAN_NUCLEON = 0.5*(gd.MASS_PROTON + gd.MASS_NEUTRON)*1.783e-33/gd.CONV_EV_TO_G
R_CONTACT_KM, R_FORWARD_KM = 12348.0, 30323.0
R0_KM, R1_KM = 1.0e4, 8.0e4
L0, L1 = R0_KM*KM, R1_KM*KM
def smoothstep(u):
    u = np.clip(np.asarray(u, dtype=float), 0.0, 1.0); return u*u*(3.0 - 2.0*u)
def rarefaction(r_km, r_shock_km):
    u = np.clip(1.0 - np.asarray(r_km, dtype=float)/r_shock_km, 0.0, 1.0)
    return np.exp((0.28 - 0.69*np.log(r_shock_km))*np.arcsin(u)**1.1)
def sn_shock_ne(width_frac, contact_jump=2.5, y_e=0.5):
    w_km = float(width_frac)*(R1_KM - R0_KM)
    def ne(l):
        r = np.asarray(l, dtype=float)/KM
        rho = 1.0e14*r**(-2.4)
        shocked = smoothstep((R_FORWARD_KM + 0.5*w_km - r)/w_km)
        factor = 1.0 + shocked*(10.0*rarefaction(r, R_FORWARD_KM) - 1.0)
        inside = smoothstep((R_CONTACT_KM + 0.5*w_km - r)/w_km)
        factor = factor*(1.0 + inside*(contact_jump - 1.0))
        out = rho*factor*gd.UNIT_G_PER_CM3/MEAN_NUCLEON*y_e
        return out[()] if np.ndim(out) == 0 else out
    return ne
def undisturbed_ne(l):
    r = np.asarray(l, dtype=float)/KM
    out = 1.0e14*r**(-2.4)*gd.UNIT_G_PER_CM3/MEAN_NUCLEON*0.5
    return out[()] if np.ndim(out) == 0 else out
def edges_for(w_km):
    w = w_km*gd.UNIT_KM
    return [r + s*w/2 for r in (12348.0*gd.UNIT_KM, 30323.0*gd.UNIT_KM) for s in (-1.0, 1.0)]
