"""Handover 2.7: does the gl_mirror prototype reproduce the shipped path?

Transcribed from the handover and exercised at even AND odd slab counts.  Only two things
differ from the text as written there: compound statements have been split across lines to
satisfy `ruff check docs/` in CI, which changes nothing semantically.

THE ODD-COUNT NUMBERS ARE NOT REPRODUCIBLE RUN TO RUN, and that is the finding rather than a
flaw in the test.  `Om = np.empty(...)` leaves the middle slab holding whatever was in that
memory, so the error there is whatever the previous occupant happened to be: 3.0e-01 on one
run, 1.8e-15 on the next, with no change to the code.  A single quiet run is therefore not
evidence the prototype is sound -- which is exactly how the bug survived, since every EVEN
count is genuinely exact to 1e-15.  Read the `worst slab` column, not the magnitude.
"""
import numpy as np
import common
import magnus.magnus as mg

L, A, H, vcc = common.chord_setup()
proj = np.zeros((3, 3))
proj[0][0] = 1.0
V0 = 4.0e-13
hv = H(0.0) - vcc(0.0)*np.diag([1.0, 0, 0])


def A_sym(l):
    l = np.asarray(l, dtype=float)
    v = V0*(1.0 + 0.5*np.cos(2.0*np.pi*l/L))          # v(l) == v(L-l)
    return -1j*(hv + v[..., None, None]*proj)


def gl_mirror(A, edges, order):
    """From HANDOVER 2.7; only compound statements split, for the linter."""
    n = edges.shape[0]
    m = n//2
    nodes = mg.gl_nodes(order)
    a, b = edges[:m, 0], edges[:m, 1]
    h = b - a
    ts = a[:, None] + np.outer(h, nodes)
    An = A(ts.ravel()).reshape(m, len(nodes), 3, 3)     # half the H evaluations
    hh = h[:, None, None]
    if order <= 2:                       # midpoint is its own mirror
        Om_f = hh*An[:, 0]
        Om_b = Om_f
    elif order <= 4:                     # nodes swap
        A1, A2 = An[:, 0], An[:, 1]
        S = 0.5*hh*(A1 + A2)
        K = (np.sqrt(3.)/12.)*(hh**2)*mg.commutator(A2, A1)
        Om_f, Om_b = S + K, S - K
    else:                                # order 6
        Om_f = mg._magnus_gl(An, h, order)
        Om_b = mg._magnus_gl(An[:, ::-1], h, order)
    Om = np.empty((n, 3, 3), dtype=complex)
    Om[:m] = Om_f
    Om[n-m:] = Om_b[::-1]
    return mg._expm_stack(Om)


print("=" * 74)
print("gl_mirror vs the shipped magnus_expansion_multislab, EVEN slab counts")
print("=" * 74)
print(f"  {'n_slabs':>8} | {'order':>5} | {'max|U_proto - U_shipped|':>26}")
print("  " + "-" * 48)
for n in (32, 64, 128):
    edges, w = common.symmetrised_edges(L, n)
    for order in (2, 4, 6):
        Us = mg.magnus_expansion_multislab(A_sym, edges, order=order,
                                           integration_method='gl',
                                           validate_input=False)
        Up = gl_mirror(A_sym, edges, order)
        print(f"  {n:>8} | {order:>5} | {np.max(np.abs(Up - Us)):>26.3e}")

print()
print("=" * 74)
print("ODD slab counts -- the middle slab has no mirror partner")
print("=" * 74)
print(f"  {'n_slabs':>8} | {'order':>5} | {'max|U_proto - U_shipped|':>26} | note")
print("  " + "-" * 68)
for n in (31, 63, 129):
    e = np.linspace(0.0, L, n + 1)
    wid = np.diff(e)
    wid = 0.5*(wid + wid[::-1])
    a = np.concatenate([[0.0], np.cumsum(wid)])
    edges = np.stack([a[:-1], a[1:]], axis=1)
    for order in (2, 4, 6):
        Us = mg.magnus_expansion_multislab(A_sym, edges, order=order,
                                           integration_method='gl',
                                           validate_input=False)
        Up = gl_mirror(A_sym, edges, order)
        diff = np.abs(Up - Us)
        mid = n//2
        note = f"worst slab = {int(np.argmax(diff.max(axis=(1,2))))} (middle={mid})"
        print(f"  {n:>8} | {order:>5} | {np.max(diff):>26.3e} | {note}")
