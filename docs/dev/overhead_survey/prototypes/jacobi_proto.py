# Prototype: batched complex-Hermitian Jacobi eigensolver + expm for d=4,5.
# Same contract as expmkernels._chN_core: (K, out, lam) -> sev (always 0.0).
import numpy as np
import numba as nb

@nb.njit(cache=True, fastmath=False, nogil=True)
def _jacobi_expm_core(K, out, lam):
    nB, d, _ = K.shape
    A = np.empty((d, d), dtype=np.complex128)
    V = np.empty((d, d), dtype=np.complex128)
    for b in range(nB):
        # copy, and Frobenius scale for the stopping test
        fro2 = 0.0
        for i in range(d):
            for j in range(d):
                A[i, j] = K[b, i, j]
                fro2 += K[b, i, j].real**2 + K[b, i, j].imag**2
                V[i, j] = 0.0 + 0.0j
            V[i, i] = 1.0 + 0.0j
        if fro2 == 0.0:
            for i in range(d):
                lam[b, i] = 0.0
                for j in range(d):
                    out[b, i, j] = 0.0 + 0.0j
                out[b, i, i] = 1.0 + 0.0j
            continue
        thr2 = (1.0e-16)**2 * fro2
        for _sweep in range(30):
            off2 = 0.0
            for p in range(d - 1):
                for q in range(p + 1, d):
                    apq = A[p, q]
                    g2 = apq.real**2 + apq.imag**2
                    off2 += g2
                    if g2 <= thr2:
                        continue
                    g = np.sqrt(g2)
                    eph = apq/g                       # e^{i phi}
                    tau = (A[q, q].real - A[p, p].real)/(2.0*g)
                    if tau >= 0.0:
                        t = 1.0/(tau + np.sqrt(1.0 + tau*tau))
                    else:
                        t = -1.0/(-tau + np.sqrt(1.0 + tau*tau))
                    c = 1.0/np.sqrt(1.0 + t*t)
                    s = t*c
                    se = s*eph                        # s e^{i phi}
                    sec = np.conj(se)                 # s e^{-i phi}
                    # columns p,q of A  (A <- A J)
                    for k in range(d):
                        akp = A[k, p]; akq = A[k, q]
                        A[k, p] = c*akp - sec*akq
                        A[k, q] = se*akp + c*akq
                    # rows p,q of A     (A <- J^H A)
                    for k in range(d):
                        apk = A[p, k]; aqk = A[q, k]
                        A[p, k] = c*apk - se*aqk
                        A[q, k] = sec*apk + c*aqk
                    # exact zero by construction; enforce to kill drift
                    A[p, q] = 0.0 + 0.0j
                    A[q, p] = 0.0 + 0.0j
                    A[p, p] = complex(A[p, p].real, 0.0)
                    A[q, q] = complex(A[q, q].real, 0.0)
                    # accumulate V <- V J
                    for k in range(d):
                        vkp = V[k, p]; vkq = V[k, q]
                        V[k, p] = c*vkp - sec*vkq
                        V[k, q] = se*vkp + c*vkq
            if off2 <= thr2:
                break
        # eigenvalues (diagonal), insertion-sorted ascending with column swaps
        for i in range(d):
            lam[b, i] = A[i, i].real
        for i in range(1, d):
            key = lam[b, i]
            # save column i of V
            for k in range(d):
                A[k, 0] = V[k, i]      # reuse A col 0 as scratch
            j = i - 1
            while j >= 0 and lam[b, j] > key:
                lam[b, j + 1] = lam[b, j]
                for k in range(d):
                    V[k, j + 1] = V[k, j]
                j -= 1
            lam[b, j + 1] = key
            for k in range(d):
                V[k, j + 1] = A[k, 0]
        # out = V diag(exp(-i lam)) V^H
        for i in range(d):
            for j in range(d):
                acc = 0.0 + 0.0j
                for m in range(d):
                    f = np.exp(-1j*lam[b, m])
                    acc += V[i, m]*f*np.conj(V[j, m])
                out[b, i, j] = acc
    return 0.0
