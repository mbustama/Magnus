import numpy as np
import numba as nb

@nb.njit(cache=True, fastmath=False, nogil=True)
def _antiherm_scale_dev(Om):
    # scale = max |Om|_elementwise (== max|K| for K = i*Om),
    # dev   = max |Om + Om^H|_elementwise (== max|K - K^H|).  Same hypot arithmetic
    # numpy uses, no temporaries.
    nB, d, _ = Om.shape
    scale = 0.0
    dev = 0.0
    for b in range(nB):
        for i in range(d):
            for j in range(d):
                z = Om[b, i, j]
                a = np.abs(z)
                if a > scale:
                    scale = a
                w = z + np.conj(Om[b, j, i])
                a = np.abs(w)
                if a > dev:
                    dev = a
    return scale, dev

@nb.njit(cache=True, fastmath=False, nogil=True)
def _gl4_omega(A1, A2, h, out):
    # out = 0.5*h*(A1+A2) + (sqrt(3)/12)*h^2*(A2@A1 - A1@A2), one pass.
    nB, d, _ = A1.shape
    c1 = np.sqrt(3.0)/12.0
    for b in range(nB):
        hb = h[b]
        f1 = 0.5*hb
        f2 = c1*hb*hb
        for i in range(d):
            for j in range(d):
                s = 0.0 + 0.0j
                for m in range(d):
                    s += A2[b, i, m]*A1[b, m, j] - A1[b, i, m]*A2[b, m, j]
                out[b, i, j] = f1*(A1[b, i, j] + A2[b, i, j]) + f2*s
    return out
