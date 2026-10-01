import numpy as np, warnings
warnings.simplefilter('ignore')
import magnus.globaldefs as gd, magnus.hamiltonians as hams, magnus.oscprob as op
o = gd.load_nufit_params('NuFIT 6.1')
U = hams.pmns_mixing_matrix(o['s12'], o['s23'], o['s13'], o['dCP']); m2=[0.0,o['D21'],o['D31']]
TEV=1e3*gd.UNIT_GEV
for delta in (1e-12, 1e-14, 1e-16, 1e-18):
  ew=eg=0
  for state in (0,1,2):
    for phase in (0.5,3.0,20.0):
        L = phase*2.0*TEV/delta
        P = np.asarray(op.osc_prob_pseudo_dirac_vacuum(TEV, L, {state: delta}, mixing_matrix=U, mass_squared=m2, average=True, average_spread=0.3))
        H0 = hams.hamiltonian_pseudo_dirac_vacuum_energy_independent(U, m2, {state: delta})
        Pg = np.asarray(op.osc_prob_energy_baseline(lambda e: H0/e, TEV, L, 0.0, None, None, True, average=True, average_spread=0.3))
        A = np.abs(U)**2; w = np.exp(-0.5*(0.3*phase)**2)
        ex = A @ A.T - np.outer(A[:, state], A[:, state])*(1.0 - w*np.cos(phase))/2.0
        ew=max(ew, np.abs(P[:3,:3]-ex).max()); eg=max(eg, np.abs(Pg[:3,:3]-ex).max())
  print('delta/m2max %.0e  wrapper err %.1e  generic err %.1e'%(delta/o['D31'], ew, eg))
