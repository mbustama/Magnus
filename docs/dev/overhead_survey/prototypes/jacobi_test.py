import sys, time
import numpy as np
sys.path.insert(0, '.')
from jacobi_proto import _jacobi_expm_core

rng = np.random.default_rng(1)
for d in (4, 5):
    nB = 1104
    X = rng.normal(size=(nB, d, d)) + 1j*rng.normal(size=(nB, d, d))
    K = (X + np.conj(np.swapaxes(X, -1, -2)))/2
    K *= 2.0    # norm ~ a few, like a Magnus slab
    out = np.empty_like(K); lam = np.empty((nB, d))
    _jacobi_expm_core(K, out, lam)     # compile
    # reference: eigh
    lam_e, V = np.linalg.eigh(K)
    Vh = np.conj(np.swapaxes(V, -1, -2))
    U_e = (V*np.exp(-1j*lam_e)[..., None, :]) @ Vh
    print("d=%d  max|U-U_eigh| = %.3e   max|lam-lam_eigh| = %.3e" %
          (d, np.max(np.abs(out - U_e)), np.max(np.abs(lam - lam_e))))
    I = np.eye(d)
    print("      unitarity jacobi %.3e   eigh %.3e" % (
        np.max(np.abs(np.conj(np.swapaxes(out,-1,-2)) @ out - I)),
        np.max(np.abs(Vh.swapaxes(-1,-2)*0 + np.conj(np.swapaxes(U_e,-1,-2)) @ U_e - I))))
    # degenerate spectra
    lam_fix = np.array([1.0, 1.0, 1.0 + 1e-9, 2.0, 3.0][:d])
    Q = np.linalg.qr(rng.normal(size=(nB, d, d)) + 1j*rng.normal(size=(nB, d, d)))[0]
    Kd = (Q*lam_fix[None, None, :]) @ np.conj(np.swapaxes(Q, -1, -2))
    Kd = (Kd + np.conj(np.swapaxes(Kd, -1, -2)))/2
    outd = np.empty_like(Kd); lamd = np.empty((nB, d))
    _jacobi_expm_core(Kd, outd, lamd)
    lam_e2, V2 = np.linalg.eigh(Kd)
    U_e2 = (V2*np.exp(-1j*lam_e2)[..., None, :]) @ np.conj(np.swapaxes(V2, -1, -2))
    print("      degenerate: max|U-U_eigh| = %.3e  unitarity %.3e" % (
        np.max(np.abs(outd - U_e2)),
        np.max(np.abs(np.conj(np.swapaxes(outd,-1,-2)) @ outd - I))))
