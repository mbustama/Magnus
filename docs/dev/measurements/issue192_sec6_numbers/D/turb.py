import numpy as np
import magnus.globaldefs as gd
import magnus.matter as matter
import magnus.hamiltonians as ham

OSC = gd.load_nufit_params('NuFIT 6.1')
E = 5.0*gd.UNIT_GEV    # Neutrino energy
RHO0 = 4.0             # Mean density, g/cm3
C = 0.03               # Mode amplitude

def one_mode(q):
    """Density along the ray: the mean, plus
    one mode of wavenumber q."""
    def rho(l):
        x = np.asarray(l, dtype=float)
        return RHO0*(1.0 + C*np.cos(q*x))
    return rho
ne0 = matter.num_density_e_func(
    0.0, lambda l: RHO0,
    density_matter_is_in_g_per_cm3=True)
V0 = matter.VCC_func(0.0, lambda l: ne0)

# H at the mean density: vacuum plus matter
H0 = np.array(ham.hamiltonian_3nu_vacuum(
    E, **OSC), dtype=complex)
H0 += ham.hamiltonian_3nu_matter(V0)

# Eigenvalues w, matter eigenvectors Vm, and
# the gap the mode has to match
w, Vm = np.linalg.eigh(H0)
d32 = abs(w[2] - w[1])
LOSC = 2*np.pi/d32     # 11 030 km
print("LOSC km", LOSC/gd.UNIT_KM)
from magnus.oscprob import (
  compute_evolution_operator_multiple_slabs
  as chain)

L = LOSC             # The region, one L_osc
rho = one_mode(d32)  # Mode tuned to the gap

def H(l):
    """H at position l, with the mode on."""
    h = np.array(ham.hamiltonian_3nu_vacuum(
        E, **OSC), dtype=complex)
    # V_CC follows the density
    return h + ham.hamiltonian_3nu_matter(
        V0*rho(l)/RHO0)

# One operator per slab, then the ordered
# product over the chain
edges = np.linspace(0.0, L, 501)
slabs = np.stack([edges[:-1], edges[1:]], 1)
Us = chain(H, slabs, 9, 6)
U = Us[0]
for u in Us[1:]:     # Earliest slab first
    U = u @ U

# Rotate to the matter basis, then read the
# transition between levels 2 and 3
P32 = abs((Vm.conj().T @ U @ Vm)[2, 1])**2
print("P32", repr(P32))
import magnus.oscprob as oscprob

kw = dict(
    nu_i=gd.NUMU, nu_f=gd.NUE,
    density_matter_is_in_g_per_cm3=True)

# Two flavors: the 1-3 sector
OSC2 = dict(sth=OSC['s13'], Dm2=OSC['D31'])
P2 = oscprob.osc_prob_matter_std_potential(
    2, one_mode(d32), E, L, OSC2, **kw)

# Three flavors
P3 = oscprob.osc_prob_matter_std_potential(
    3, one_mode(d32), E, L, OSC, **kw)

# The same call with a sterile state added
OSC4 = dict(OSC, s14=np.sqrt(0.10), D41=1.0,
            s24=np.sqrt(0.10), s34=0.0,
            d14=0.0, d24=0.0)
P4 = oscprob.osc_prob_matter_std_potential(
    4, one_mode(d32), E, L, OSC4, **kw)
print('P2 P3 P4', P2, P3, P4)
# gaps at 2 and 4 flavours
import magnus.hamiltonians as hh
h2 = None
try:
    H2 = np.array(hh.hamiltonian_2nu_vacuum(E, **OSC2), dtype=complex) + hh.hamiltonian_2nu_matter(V0)
    w2=np.linalg.eigvalsh(H2); print('2nu gap rel', abs(w2[1]-w2[0])/d32 - 1)
except Exception as e: print('2nu err', e)
try:
    H4 = np.array(hh.hamiltonian_4nu_vacuum(E, **{k:v for k,v in OSC4.items()}), dtype=complex) + hh.hamiltonian_4nu_matter(V0)
    w4=np.linalg.eigvalsh(H4); print('4nu gaps rel', [abs(w4[i+1]-w4[i])/d32-1 for i in range(3)])
except Exception as e: print('4nu err', e)
