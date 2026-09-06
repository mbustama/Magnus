import numpy as np
import numba as nb

@nb.njit(cache=True, fastmath=False, nogil=True)
def _antiherm_scale_dev2(Om):
    # squared magnitudes; sqrt once at the end.  max is order-preserved.
    nB, d, _ = Om.shape
    s2 = 0.0
    d2 = 0.0
    for b in range(nB):
        for i in range(d):
            for j in range(d):
                z = Om[b, i, j]
                a = z.real*z.real + z.imag*z.imag
                if a > s2:
                    s2 = a
                w = z + np.conj(Om[b, j, i])
                a = w.real*w.real + w.imag*w.imag
                if a > d2:
                    d2 = a
    return np.sqrt(s2), np.sqrt(d2)

@nb.njit(cache=True, fastmath=False, nogil=True)
def _gl4_omega_from_At(At, h, out):
    # At: (nB, 2, d, d) contiguous; out = GL4 Omega.  Returns 1 if every
    # slab's two node samples were bit-identical (constant-A fast path taken:
    # out = h*A1), else 0 with the full expression.
    nB, m, d, _ = At.shape
    c1 = np.sqrt(3.0)/12.0
    # early-exit equality scan (cheap on smooth profiles: exits immediately)
    identical = True
    for b in range(nB):
        for i in range(d):
            for j in range(d):
                if At[b, 0, i, j] != At[b, 1, i, j]:
                    identical = False
                    break
            if not identical:
                break
        if not identical:
            break
    if identical:
        for b in range(nB):
            hb = h[b]
            for i in range(d):
                for j in range(d):
                    out[b, i, j] = hb*At[b, 0, i, j]
        return 1
    for b in range(nB):
        hb = h[b]
        f1 = 0.5*hb
        f2 = c1*hb*hb
        for i in range(d):
            for j in range(d):
                s = 0.0 + 0.0j
                for k in range(d):
                    s += At[b, 1, i, k]*At[b, 0, k, j] - At[b, 0, i, k]*At[b, 1, k, j]
                out[b, i, j] = f1*(At[b, 0, i, j] + At[b, 1, i, j]) + f2*s
    return 0
