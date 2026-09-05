# -*- coding: utf-8 -*-
# SPDX-License-Identifier: GPL-3.0-only
# Copyright (C) 2026 Mauricio Bustamante
r"""expmkernels.py

Compiled kernels for :math:`\exp(-iK)`, K Hermitian.

This module is the ``'numba'`` backend of ``magnus.magnus._expm_stack``.  For
2x2 and 3x3 matrices it computes the matrix exponential of a stack of small
Hermitian matrices without an eigenvector solver, by applying to :math:`K` the
polynomial that interpolates :math:`\exp(-i\lambda)` on the spectrum of
:math:`K`:

.. math::

   \exp(-iK) = a_0 I + a_1 K + a_2 K^2 .

Cayley-Hamilton guarantees such a polynomial exists (degree :math:`d-1` for a
:math:`d \times d` matrix); the eigenvalues are obtained in closed form.

Why this is worth a compiled kernel
-----------------------------------

``np.linalg.eigh`` costs about 1.25 us per 3x3 *regardless of stack size* --
1, 108 or 4096 matrices, the per-matrix cost is flat -- because it loops over
LAPACK internally rather than vectorizing over the stack.  On a 108-slab Magnus
pass that single call is roughly a quarter of the total.

The same algebra written in pure numpy does not help below stacks of about a
hundred: it is some twenty numpy calls, each paying dispatch overhead on
arithmetic that is otherwise trivial, so it *loses* at small stacks and wins
only mildly at large ones.  Only a compiled kernel removes the dispatch, which
is why numba is used here and why a numpy version of these formulas is not
offered as a third backend.

Dimensions 4 and 5: Jacobi, not a closed form
---------------------------------------------

There is no practical closed form for the eigenvalues of a 4x4 or 5x5 Hermitian
matrix, and for a long time that sentence ended "so 4nu and 5nu keep the
``eigh`` path".  The conclusion did not follow: what made those dimensions slow
was never the missing closed form but ``eigh``'s fixed per-matrix LAPACK
overhead, about 2.3 us on a 4x4 -- two thirds of a whole d=4 Magnus pass.  So
4x4 and 5x5 stacks go to :func:`_jacobi_expm_core` instead, a batched cyclic
Jacobi eigensolver that warm-starts each matrix from its predecessor's
eigenvectors and re-orthonormalizes that basis at every step; see its docstring
for the scheme and for which of its details are load-bearing.  Unlike the
closed forms it is iterative, so this backend replacement is *not* bit-identical
to what it replaces -- it is held to the same accuracy class as ``eigh``
instead (within 5.5x at every norm, clustering and degeneracy measured, against
the same yardstick that admits the 3x3 closed form at up to 10x).
:func:`supports_dim` is the single place that decides the routing.

Why the interpolation form is safe at a degeneracy
--------------------------------------------------

Coincident eigenvalues are the whole numerical risk in a Cayley-Hamilton
scheme, because the interpolation coefficients divide by eigenvalue
differences.  Two facts remove it here, and both are load-bearing enough to
state:

*A Hermitian matrix is never defective.*  Its minimal polynomial has simple
roots even when its characteriztic polynomial does not, so a polynomial
matching :math:`\exp(-i\lambda)` on the *distinct* eigenvalues already
reproduces the function exactly.  The confluent (Hermite) form, which matches
derivatives as well and is unavoidable for a general matrix, is not needed for
this one.  Nothing in this module differentiates anything.

*The ill-conditioned coefficient multiplies a correspondingly small matrix.*
Write the interpolant in Newton form on eigenvalues sorted ascending, with the
spectrum shifted so the median eigenvalue sits at zero:

.. math::

   \exp(-iZ) = f[z_0] I + f[z_0, z_1](Z - z_0 I)
               + f[z_0, z_1, z_2](Z - z_0 I)(Z - z_1 I) ,

with :math:`Z = K - \lambda_1 I` and :math:`z_0 \le z_1 = 0 \le z_2`.  The
first divided difference is evaluated as
:math:`f[a, b] = -i e^{-i(a+b)/2}\, \mathrm{sinc}((a-b)/2)`, which is
cancellation-free for every pair including :math:`a = b` (see ``_sinc``).
The second, :math:`(f[z_0,z_1] - f[z_1,z_2])/(z_0 - z_2)`, does lose digits as
the nodes coalesce -- its absolute error grows like
:math:`\epsilon/(z_2 - z_0)` -- but the matrix it multiplies has norm at most
:math:`(z_2 - z_0)^2`, so the product's error is bounded by
:math:`\epsilon\,(z_2 - z_0)` and *vanishes* with the gap.  Sorting is what
makes this true: it is what guarantees that a small :math:`z_2 - z_0` means all
three eigenvalues are close, rather than one unlucky pair out of three.

So there is no tolerance, no crossover, and no near-degenerate branch to place
correctly.  The only guard is for :math:`z_0 = z_2` exactly, where the term is
multiplied by the zero matrix and is simply dropped.

Measured against ``scipy.linalg.expm``, the error is 1e-16 at splittings of
1e-2, 1e-6, 1e-10, 1e-14 and exactly zero alike.  The closed-form
*eigenvalues* are much worse than that near a degeneracy -- they degrade to
~1e-9, because :math:`\arccos` has infinite derivative at the ends of its
range, which is where a repeated root sits -- and at
:math:`\lVert K \rVert \sim 1` it does not matter: the interpolation error is
*second* order in the displacement of a coalescing node, so a node that is 1e-9
off contributes 1e-18.  A test asserts both halves of that.

**That argument has a range of validity, and an earlier version of this
paragraph did not say so.**  The eigenvalue error scales with the norm, so the
second-order suppression is fighting a term that grows: where a clustered
spectrum meets a large norm the closed form reaches 2.7e-07 against ``eigh``'s
3.0e-11, a factor of 7440.  Neither a sweep over norms at generic separation
(ratios 0.4-2.5) nor a sweep over separations at norm 1 (0.4-1.0) visits that
corner, which is how the unqualified claim came to be written and believed.
:data:`SEV_TOL` is the gate that keeps it out of reach, and
``tests/test_expm_backend.py`` now crosses the two axes so the corner cannot go
unmeasured again.

A note on the determinant
-------------------------

The cross term of :math:`\det X` for Hermitian X is
:math:`2\,\mathrm{Re}(X_{01} X_{12} \overline{X_{02}})`.  Moving that conjugate
to either of the other two factors produces matrices that are still nearly
unitary and still plausible, and wrong by O(1).  An earlier prototype of this
kernel shipped that exact transposition and reported a 6x speed-up while
returning garbage; ``_ch3_core`` carries the term explicitly and
``tests/test_expm_backend.py`` pins it against ``np.linalg.det``.

.. versionadded:: 1.0.0

Routine listings
----------------

    * HAVE_NUMBA - Whether the compiled kernels are available
    * SEV_TOL - Conditioning above which eigh answers instead
    * supports_dim - Whether a given matrix dimension has a kernel
    * expm_herm_stack - exp(-iK) and the eigenvalues of K, for a stack
"""

__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"


import math

import numpy as np

try:
    import numba as nb
except ImportError:                                     # pragma: no cover
    nb = None

HAVE_NUMBA = nb is not None
r"""bool: Whether numba imported, and so whether the compiled kernels exist.

False leaves every backend decision to ``eigh``; nothing else in the package
changes.  numba is an optional dependency (``pip install magnuspy[fast]``).
"""


_TWO_PI_3 = 2.0*np.pi/3.0

# Below this the two-term Taylor series for sin(x)/x is already accurate to
# better than 1e-17 relative (the first neglected term is x^6/5040), and above
# it the direct quotient has nothing to cancel.  There is no accuracy notch at
# the crossover, which is why the threshold does not need tuning.
_SINC_SMALL = 1.e-4


def _sinc(x: float) -> float:
    r"""Returns :math:`\sin(x)/x`, including at :math:`x = 0`.

    Not :func:`numpy.sinc`, which carries a factor :math:`\pi`.  This is the
    unnormalized cardinal sine, and it is the reason the divided differences of
    :math:`\exp` in this module need no degenerate branch: the first divided
    difference of :math:`e^{-i\lambda}` over nodes a and b is
    :math:`-i e^{-i(a+b)/2}\,\mathrm{sinc}((a-b)/2)`, an expression in which
    :math:`a \to b` is an ordinary limit rather than a division by zero.

    Parameters
    ----------
    x : float
        Argument.

    Returns
    -------
    float
        :math:`\sin(x)/x`; 1.0 at x = 0.
    """
    if abs(x) < _SINC_SMALL:
        x2 = x*x
        return 1.0 - x2/6.0*(1.0 - x2/20.0)
    return math.sin(x)/x


SEV_TOL = 1.0e4
r"""float: Above this :math:`m = \mathrm{tr}(X^2)/6`, a 3x3 is handed back to ``eigh``.

The closed-form solve is used only in the range of :math:`\lVert K \rVert` where it is
*measured* to behave, and ``eigh`` answers beyond it.  What the gate admits, worst over two
spectrum families x 6 separations x 40 random bases per rung
(``docs/dev/calibrate_sev_tol.py``):

.. code-block:: text

    m         worst |closed - expm|    eigh, same cells
    4.4e1     1.0e-14                  1.4e-15
    4.0e2     8.1e-14                  2.0e-14
    1.1e3     2.1e-13                  3.9e-14
    4.4e3     8.5e-13                  7.8e-14
    1.0e4     2.0e-12                  6.1e-14

The worst measured is 2.0e-12 just under the gate and 2.3e-13 in the
:math:`m \le 1.1\times10^3` corner, so the guarantee is stated with headroom at **5e-12
absolute across everything the gate admits**, and **5e-13** in that corner.  The headroom is
deliberate: these are worst-over-random-bases quantities and more sampling keeps finding
slightly worse ones, which is exactly how the previous claim came to be false.  Both bounds
are far below any tolerance this package is asked for, and far below where the closed form
runs away past the gate: cells at :math:`m \ge 1.1\times10^5` reach 131x ``eigh``, and at
:math:`m \sim 4\times10^9`, 7440x.

An earlier version of this docstring claimed 2e-13 across the admitted range.  That was
measured on one corner of it and is not true of the rest; at :math:`m = 1.1\times10^3` itself
about 1% of random bases exceed it (11 of 1200, worst 2.3e-13).  No value of this constant
could have rescued that claim, because ``test_sev_tol_sits_inside_its_calibrated_window`` pins
:math:`m = 1.1\times10^3` as a cell that must stay on the kernel, so the gate cannot be
lowered past the point where the claim already fails.  The number was corrected instead.

**Read m, never "spectral scale".**  Two calibrations of this constant appeared to contradict
each other -- one finding the first unsafe cell at :math:`1.1\times10^5`, the other at
:math:`4.4\times10^3` -- purely because they used different spectrum families and both called
the result "scale :math:`10^2`".  :math:`[-s, -s(1-d), s]` spans :math:`2s` and gives
:math:`m \simeq 0.44\,s^2`; :math:`[0, d, S]` spans :math:`S` and gives
:math:`m \simeq 0.11\,S^2` -- a factor of four in :math:`m` at the same nominal scale.
Compared at equal :math:`m` the two families agree to within their sampling scatter, and the
disagreement dissolves.  :math:`m = \mathrm{tr}(X^2)/6` is a spectral invariant; the word
"scale" is not, and is what made this look like a contradiction.

**Why the scale and not the clustering, which is the actual mechanism.**  The damage needs a
clustered spectrum *and* a large norm together: :math:`\arccos` has infinite derivative at
:math:`u = \pm 1`, so clustering turns rounding in :math:`u` into an eigenvalue error
:math:`\sim\sqrt{\epsilon}\,\lVert K \rVert`, which only matters once the norm is large.  But
the clustering half cannot be gated on, because the danger is a **band** rather than a tail:
at *exact* degeneracy the pair comes out bit-identical and the answer is fine (measured 0.3-1.8x
of ``eigh`` at :math:`u = \pm 1` exactly), the damage sits at intermediate separations, and
:math:`1/(1-u^2)` is largest exactly where there is no problem.  A one-sided threshold on it
therefore cannot work -- verified by calibration, which found no separating value.

So this gate is deliberately **conservative rather than tight**: it also declines large-norm
spectra that are *not* clustered and would have been fine (measured 0.6x of ``eigh``).  That
costs speed on those, never accuracy, and it costs nothing where the speed comes from -- a
Magnus slab has :math:`\lVert\Omega\rVert \lesssim \pi` by construction, so slab chains are
never declined, and an ordinary 3nu constant-density or vacuum call measures
:math:`\lVert K \rVert \approx 4`.  What it does decline is the large accumulated phase of an
eV-scale sterile splitting, where accuracy is worth more than the microsecond.

.. versionadded:: 1.0.0
"""


def _ch2_core(K, out, lam):
    r"""Fills ``out`` with :math:`\exp(-iK)` and ``lam`` with K's eigenvalues, d = 2.

    Removing the trace leaves :math:`X = K - \mu I` traceless, and a traceless
    Hermitian 2x2 satisfies :math:`X^2 = r^2 I` with
    :math:`r = \lVert X \rVert`, so no eigenvalue solver is needed at all:

    .. math::

       \exp(-iK) = e^{-i\mu}\left[\cos r\, I
                   - i\, \mathrm{sinc}(r)\, X\right] .

    Uniformly stable, with no branch beyond the one inside :func:`_sinc`: the
    degenerate case is r = 0, and :math:`\mathrm{sinc}(0) = 1` covers it.

    Only the *lower* triangle of each input matrix is read, which is what
    ``np.linalg.eigh`` reads too (its ``UPLO`` defaults to ``'L'``).  That is
    deliberate rather than incidental: the caller admits input that is
    anti-Hermitian only to 1e-12 relative, so a kernel reading the other
    triangle would be exponentiating a different matrix from the backend it
    replaces, and the two would disagree by ~2e-12 on such input.  Reading the
    same triangle makes the backends interchangeable instead.

    Parameters
    ----------
    K : np.ndarray
        Hermitian matrices, shape (n, 2, 2), complex, C-contiguous.
    out : np.ndarray
        Output buffer for :math:`\exp(-iK)`, shape (n, 2, 2), complex.
    lam : np.ndarray
        Output buffer for the eigenvalues, ascending, shape (n, 2), float.

    Returns
    -------
    float
        Always 0.0.  A 2x2 has no characteriztic cubic and no ``arccos``, so it has nothing
        to be ill-conditioned about: measured against a 60-digit reference it tracks ``eigh``
        at every scale from 1 to 1e5 and at every eigenvalue separation down to exact
        degeneracy (ratios 0.4-1.0).  The return value exists only so both kernels share a
        calling convention with ``_ch3_core``, which does need it.
    """
    for i in range(K.shape[0]):
        k00 = K[i, 0, 0].real
        k11 = K[i, 1, 1].real
        k01 = K[i, 1, 0].conjugate()        # lower triangle, as eigh reads it

        mu = 0.5*(k00 + k11)
        xd = k00 - mu                       # X[0,0]; X[1,1] = -xd by tracelessness
        r = math.sqrt(xd*xd + k01.real*k01.real + k01.imag*k01.imag)

        lam[i, 0] = mu - r
        lam[i, 1] = mu + r

        c = math.cos(r)
        s = _sinc(r)
        ph = complex(math.cos(mu), -math.sin(mu))       # e^{-i mu}

        out[i, 0, 0] = ph*complex(c, -s*xd)
        out[i, 1, 1] = ph*complex(c, s*xd)
        w = ph*complex(0.0, -s)                         # -i sinc(r) e^{-i mu}
        out[i, 0, 1] = w*k01
        out[i, 1, 0] = w*k01.conjugate()

    return 0.0


def _ch3_core(K, out, lam):
    r"""Fills ``out`` with :math:`\exp(-iK)` and ``lam`` with K's eigenvalues, d = 3.

    Eigenvalues come from the trigonometric solution of the characteriztic
    cubic of the traceless part :math:`X = K - (\mathrm{tr}K/3) I`, whose
    characteriztic polynomial is :math:`y^3 - 3my - 2n` with
    :math:`m = \mathrm{tr}(X^2)/6` and :math:`n = \det(X)/2`.  Substituting
    :math:`y = 2\sqrt{m}\cos\theta` turns it into
    :math:`\cos 3\theta = n/m^{3/2}`, so

    .. math::

       y = 2\sqrt{m}\,\cos\!\left(\theta + \tfrac{2\pi k}{3}\right) ,
       \qquad \theta = \tfrac{1}{3}\arccos\frac{n}{m^{3/2}} ,
       \qquad k = 0, 1, 2 .

    Because :math:`\arccos` returns :math:`[0, \pi]`, the angle :math:`\theta`
    lies in :math:`[0, \pi/3]`, which confines each root to its own interval and
    so fixes the order in advance.  Written the way the code below writes them:
    :math:`\theta` gives :math:`\cos \in [\tfrac12, 1]`, the **largest**;
    :math:`\theta - 2\pi/3` gives :math:`\cos \in [-\tfrac12, \tfrac12]`, the
    **middle**; and :math:`\theta + 2\pi/3` gives :math:`\cos \in [-1,
    -\tfrac12]`, the **smallest**.  So they can be stored ascending without a
    sort.

    Take the three angles as written, not the index k: with the
    :math:`+2\pi k/3` form above it is :math:`k = 1` that is smallest and
    :math:`k = 2` that is the middle root, and an earlier draft of this
    paragraph had that backwards.  Anyone reconciling the code to the prose
    should reconcile it to the three displayed angles.

    The ordering is asserted against ``np.linalg.eigvalsh`` by the tests rather
    than enforced here by a sort, so that a slip in the formula shows up as a
    failure instead of being tidied away.

    The exponential is then the Newton interpolant of :math:`e^{-i\lambda}` on
    those three nodes, shifted to put the median at the origin; see the module
    docstring for why that shift plus the ascending order is what makes the
    scheme stable through a degeneracy.  Collecting the Newton form in powers
    of :math:`Z = K - \lambda_1 I` leaves

    .. math::

       \exp(-iZ) = (d_0 + d_1 g_1) I + (d_1 + d_2 g_1) Z + d_2 Z^2 ,

    so only :math:`Z^2` is formed, and only its upper triangle, since Z is
    Hermitian.  The result is not, so all nine entries are written.

    Only the *lower* triangle of each input matrix is read, which is the
    triangle ``np.linalg.eigh`` reads as well; see :func:`_ch2_core` for why
    that matters.

    Parameters
    ----------
    K : np.ndarray
        Hermitian matrices, shape (n, 3, 3), complex, C-contiguous.
    out : np.ndarray
        Output buffer for :math:`\exp(-iK)`, shape (n, 3, 3), complex.
    lam : np.ndarray
        Output buffer for the eigenvalues, ascending, shape (n, 3), float.

    Returns
    -------
    float
        The largest :math:`m = \mathrm{tr}(X^2)/6` over the stack, i.e. the spectral scale of
        the worst-conditioned matrix in it.  When this exceeds :data:`SEV_TOL` the caller
        recomputes the stack with ``eigh``; see that constant for the calibration, and for why
        the gate is on the scale even though the mechanism is the clustering.

        Costs one compare per matrix on a quantity already in hand, is reduced to one number
        here so the caller needs no ``numpy`` reduction, and is *returned* rather than written
        to a buffer so the common path allocates nothing.
    """
    sev_max = 0.0
    for i in range(K.shape[0]):
        k00 = K[i, 0, 0].real
        k11 = K[i, 1, 1].real
        k22 = K[i, 2, 2].real
        p01 = K[i, 1, 0].conjugate()        # lower triangle, as eigh reads it
        p02 = K[i, 2, 0].conjugate()
        p12 = K[i, 2, 1].conjugate()

        c0 = (k00 + k11 + k22)/3.0
        a = k00 - c0                        # the traceless part's diagonal
        b = k11 - c0
        c = k22 - c0
        n01 = p01.real*p01.real + p01.imag*p01.imag
        n02 = p02.real*p02.real + p02.imag*p02.imag
        n12 = p12.real*p12.real + p12.imag*p12.imag

        m = (a*a + b*b + c*c + 2.0*(n01 + n02 + n12))/6.0
        # det(X).  The cross term's conjugate belongs on p02 and nowhere else;
        # see the module docstring on what moving it costs.
        det = (a*b*c - a*n12 - b*n02 - c*n01
               + 2.0*(p01*p12*p02.conjugate()).real)

        if m > 0.0:
            sqm = math.sqrt(m)
            den = m*sqm
            # den underflows only for ||X|| below ~1e-100, where the
            # eigenvalues are negligible against exp() whatever this returns;
            # the guard is here so it cannot become a NaN.
            u = 0.5*det/den if den > 0.0 else 0.0
            if u > 1.0:
                u = 1.0
            elif u < -1.0:
                u = -1.0
            th = math.acos(u)/3.0
            y_hi = 2.0*sqm*math.cos(th)
            y_md = 2.0*sqm*math.cos(th - _TWO_PI_3)
            y_lo = 2.0*sqm*math.cos(th + _TWO_PI_3)
            # Report the spectral scale, which is what SEV_TOL gates on.  The mechanism is
            # the arccos derivative blowing up at u = +/-1, but that half cannot be gated on:
            # see SEV_TOL for why the clustering forms a band rather than a tail.
            if m > sev_max:
                sev_max = m
        else:
            y_hi = 0.0                      # X = 0: K is a multiple of I
            y_md = 0.0
            y_lo = 0.0

        lam[i, 0] = c0 + y_lo
        lam[i, 1] = c0 + y_md
        lam[i, 2] = c0 + y_hi

        g1 = y_md - y_lo                    # the two gaps, both >= 0
        g2 = y_hi - y_md

        # Divided differences of e^{-iz} on the shifted nodes (-g1, 0, g2).
        d0 = complex(math.cos(g1), math.sin(g1))                    # f[z0]
        h1 = 0.5*g1
        h2 = 0.5*g2
        d1 = complex(0.0, -_sinc(h1))*complex(math.cos(h1), math.sin(h1))
        dr = complex(0.0, -_sinc(h2))*complex(math.cos(h2), -math.sin(h2))
        tot = g1 + g2
        # tot == 0 means all three eigenvalues coincide, so Z is the zero
        # matrix and the term this coefficient multiplies is zero anyway.
        d2 = (d1 - dr)/(-tot) if tot > 0.0 else complex(0.0, 0.0)

        cI = d0 + d1*g1
        cZ = d1 + d2*g1
        cZ2 = d2

        za = a - y_md                       # Z = K - (c0 + y_md) I
        zb = b - y_md
        zc = c - y_md

        s00 = za*za + n01 + n02             # upper triangle of Z^2
        s11 = n01 + zb*zb + n12
        s22 = n02 + n12 + zc*zc
        s01 = p01*(za + zb) + p02*p12.conjugate()
        s02 = p02*(za + zc) + p01*p12
        s12 = p12*(zb + zc) + p01.conjugate()*p02

        mu = c0 + y_md
        ph = complex(math.cos(mu), -math.sin(mu))       # e^{-i mu}

        out[i, 0, 0] = ph*(cI + cZ*za + cZ2*s00)
        out[i, 1, 1] = ph*(cI + cZ*zb + cZ2*s11)
        out[i, 2, 2] = ph*(cI + cZ*zc + cZ2*s22)
        out[i, 0, 1] = ph*(cZ*p01 + cZ2*s01)
        out[i, 1, 0] = ph*(cZ*p01.conjugate() + cZ2*s01.conjugate())
        out[i, 0, 2] = ph*(cZ*p02 + cZ2*s02)
        out[i, 2, 0] = ph*(cZ*p02.conjugate() + cZ2*s02.conjugate())
        out[i, 1, 2] = ph*(cZ*p12 + cZ2*s12)
        out[i, 2, 1] = ph*(cZ*p12.conjugate() + cZ2*s12.conjugate())

    return sev_max


# The Jacobi sweep cap.  Cyclic Jacobi on a Hermitian matrix converges
# unconditionally, so this is a backstop rather than a tuning knob: hitting it
# makes _jacobi_expm_core report an infinite severity and the caller recomputes
# with eigh.  With the floor below calibrated above the rounding equilibrium,
# no input has been measured to reach it.
_JACOBI_MAX_SWEEPS = 30

# Rotation threshold for a single squared off-diagonal element, relative to
# the squared Frobenius norm of the matrix: an element below it is left alone,
# and a sweep that rotates nothing means every element is below it, which is
# the termination criterion.  Converged therefore means each off-diagonal
# element is at most 1e-15 of the matrix norm -- the same class as ``eigh``'s
# own backward error -- and termination is *guaranteed* well before the sweep
# cap, because every rotation removes at least twice this much from the
# squared off-norm while re-injecting only ~eps^2 = 4.9e-32 of it as rounding
# noise, a 40x margin.  Two prototype choices lost that guarantee, and both
# turned convergence into a lottery whose losers would each have silently
# sent their whole stack back to eigh through the severity flag: a threshold
# of (1e-16)^2 = 1e-32 sits inside the noise equilibrium itself (measured at
# up to 4.4e-32, so sub-threshold elements were unreachable), and checking
# the *sum* of squared elements against the same constant that gates each
# single one lets a sweep rotate nothing while the sum still fails -- a
# thrash with no exit.
_JACOBI_OFF2_REL = 1.0e-30


def _jacobi_expm_core(K, out, lam):
    r"""Fills ``out`` with :math:`\exp(-iK)` and ``lam`` with K's eigenvalues, d = 4 or 5.

    One kernel for both dimensions -- every loop bound comes from the shape --
    and nothing in it is specific to 4 or 5 beyond the guard in
    :func:`expm_herm_stack`.  Unlike :func:`_ch2_core` and :func:`_ch3_core`
    there is no closed form here: eigenvalues and eigenvectors come from a
    cyclic complex-Hermitian Jacobi iteration, and the exponential is
    reconstructed as :math:`U = \sum_j e^{-i\lambda_j} v_j v_j^\dagger`.  That
    is the same spectral reconstruction the ``eigh`` route uses, with the
    eigensolver swapped for one that has no LAPACK per-matrix overhead --
    which, at ~2.3 us per 4x4 call, is what made ``eigh`` two thirds of a d=4
    Magnus pass.

    Three details are load-bearing:

    * **Warm start.**  Each matrix's iteration starts from the previous
      matrix's converged eigenvectors (``A = V0^H K V0``).  Consecutive
      matrices are consecutive slabs of the same energy, so ``A`` arrives
      nearly diagonal and the sweep count drops by roughly half; at an energy
      boundary in the flattened stack the warm start is merely less effective
      for one matrix, and correctness never depends on it.
    * **Modified Gram-Schmidt on the saved basis, for every matrix.**  Without
      it, warm-start non-unitarity compounds linearly along the chain:
      measured 2.3e-11 against ``eigh`` after 13k chained matrices, against
      3.9e-14 with it.  The reconstruction reads the re-orthonormalized basis
      too, which is what keeps per-slab unitarity at ``eigh``'s level.
    * **Rotations act on pre-rotation values.**  Each 2x2 rotation reads both
      of the entries it updates before writing either, and the accumulating
      sums are left exactly as written -- reassociating them is what
      ``fastmath`` would license, and why it stays off.

    The eigenvalues are insertion-sorted ascending to honor the
    :func:`expm_herm_stack` contract; the basis needs no matching reorder,
    because the reconstruction above sums over all of its columns and is
    permutation-invariant.

    Only the *lower* triangle of each input matrix is read, which is what
    ``np.linalg.eigh`` reads too; see :func:`_ch2_core` for why the backends
    must agree on that.

    Parameters
    ----------
    K : np.ndarray
        Hermitian matrices, shape (n, d, d) with d = 4 or 5, complex,
        C-contiguous.
    out : np.ndarray
        Output buffer for :math:`\exp(-iK)`, shape (n, d, d), complex.
    lam : np.ndarray
        Output buffer for the eigenvalues, ascending, shape (n, d), float.

    Returns
    -------
    float
        0.0, or ``inf`` if any matrix hit :data:`_JACOBI_MAX_SWEEPS` without
        converging, in which case the caller recomputes the stack with
        ``eigh`` -- the same hook :data:`SEV_TOL` serves for the 3x3 kernel.
        There is no conditioning gate here because there is no conditioning
        cliff to gate: Jacobi is backward stable at every norm and clustering
        measured (worst 5.5x of ``eigh``, against the closed form's 7440x),
        so the only decline is the sweep cap -- a backstop no measured input
        reaches, and see :data:`_JACOBI_OFF2_REL` for the calibration that
        keeps it that way.
    """
    nB = K.shape[0]
    d = K.shape[1]
    A = np.empty((d, d), dtype=np.complex128)
    V = np.empty((d, d), dtype=np.complex128)
    V0 = np.empty((d, d), dtype=np.complex128)   # previous MGS'd eigenvectors
    T = np.empty((d, d), dtype=np.complex128)
    W = np.empty((d, d), dtype=np.complex128)
    f = np.empty(d, dtype=np.complex128)
    have_warm = False
    sev = 0.0
    for b in range(nB):
        # Hermitize from the lower triangle into T, accumulating the squared
        # Frobenius norm from the same reads.
        fro2 = 0.0
        for i in range(d):
            x = K[b, i, i].real
            T[i, i] = complex(x, 0.0)
            fro2 += x*x
            for j in range(i):
                pij = K[b, i, j]
                T[i, j] = pij
                T[j, i] = pij.conjugate()
                fro2 += 2.0*(pij.real*pij.real + pij.imag*pij.imag)
        if fro2 == 0.0:                          # K = 0: exp(-iK) = I exactly
            for i in range(d):
                lam[b, i] = 0.0
                for j in range(d):
                    out[b, i, j] = complex(0.0, 0.0)
                out[b, i, i] = complex(1.0, 0.0)
            continue
        if have_warm:
            # A = V0^H T V0, V = V0
            for i in range(d):
                for j in range(d):
                    acc = complex(0.0, 0.0)
                    for m in range(d):
                        acc += T[i, m]*V0[m, j]
                    W[i, j] = acc
            for i in range(d):
                for j in range(d):
                    acc = complex(0.0, 0.0)
                    for m in range(d):
                        acc += V0[m, i].conjugate()*W[m, j]
                    A[i, j] = acc
                    V[i, j] = V0[i, j]
            # re-hermitize the diagonal (rounding)
            for i in range(d):
                A[i, i] = complex(A[i, i].real, 0.0)
        else:
            for i in range(d):
                for j in range(d):
                    A[i, j] = T[i, j]
                    V[i, j] = complex(0.0, 0.0)
                V[i, i] = complex(1.0, 0.0)
        thr2 = _JACOBI_OFF2_REL*fro2
        converged = False
        for _sweep in range(_JACOBI_MAX_SWEEPS):
            rotated = False
            for p in range(d - 1):
                for q in range(p + 1, d):
                    apq = A[p, q]
                    g2 = apq.real*apq.real + apq.imag*apq.imag
                    if g2 <= thr2:
                        continue
                    rotated = True
                    g = math.sqrt(g2)
                    # X*(1.0/c), not X/c: numba divides a complex by a real
                    # componentwise where numpy multiplies by the reciprocal,
                    # and the uncompiled form of this source must agree with
                    # the compiled one bitwise.
                    eph = apq*(1.0/g)
                    tau = (A[q, q].real - A[p, p].real)/(2.0*g)
                    if tau >= 0.0:
                        t = 1.0/(tau + math.sqrt(1.0 + tau*tau))
                    else:
                        t = -1.0/(-tau + math.sqrt(1.0 + tau*tau))
                    c = 1.0/math.sqrt(1.0 + t*t)
                    s = t*c
                    se = s*eph
                    sec = se.conjugate()
                    for k in range(d):
                        akp = A[k, p]
                        akq = A[k, q]
                        A[k, p] = c*akp - sec*akq
                        A[k, q] = se*akp + c*akq
                    for k in range(d):
                        apk = A[p, k]
                        aqk = A[q, k]
                        A[p, k] = c*apk - se*aqk
                        A[q, k] = sec*apk + c*aqk
                    A[p, q] = complex(0.0, 0.0)
                    A[q, p] = complex(0.0, 0.0)
                    A[p, p] = complex(A[p, p].real, 0.0)
                    A[q, q] = complex(A[q, q].real, 0.0)
                    for k in range(d):
                        vkp = V[k, p]
                        vkq = V[k, q]
                        V[k, p] = c*vkp - sec*vkq
                        V[k, q] = se*vkp + c*vkq
            if not rotated:
                converged = True
                break
        if not converged:
            sev = np.inf
        # Save V as the next warm start, re-orthonormalized by modified
        # Gram-Schmidt so non-unitarity cannot compound along the chain.
        for i in range(d):
            for j in range(d):
                V0[i, j] = V[i, j]
        for j in range(d):
            for m in range(j):
                acc = complex(0.0, 0.0)
                for k in range(d):
                    acc += V0[k, m].conjugate()*V0[k, j]
                for k in range(d):
                    V0[k, j] -= acc*V0[k, m]
            nrm = 0.0
            for k in range(d):
                vkj = V0[k, j]
                nrm += vkj.real*vkj.real + vkj.imag*vkj.imag
            nrm = 1.0/math.sqrt(nrm)
            for k in range(d):
                V0[k, j] *= nrm
        have_warm = True
        # U = sum_m e^{-i lam_m} v_m v_m^H, from the MGS'd basis.  The phases
        # are computed straight-line rather than in a loop over m, and that is
        # load-bearing: LLVM vectorizes a trig loop into SVML's 1-ulp vector
        # routines where straight-line calls go to libm, and the uncompiled
        # form of this source must agree with the compiled one bitwise.  (The
        # arithmetic loops are safe either way: with fastmath off the compiler
        # may not reassociate their reductions, so it cannot vectorize them.)
        f[0] = complex(math.cos(A[0, 0].real), -math.sin(A[0, 0].real))
        f[1] = complex(math.cos(A[1, 1].real), -math.sin(A[1, 1].real))
        f[2] = complex(math.cos(A[2, 2].real), -math.sin(A[2, 2].real))
        f[3] = complex(math.cos(A[3, 3].real), -math.sin(A[3, 3].real))
        if d == 5:
            f[4] = complex(math.cos(A[4, 4].real), -math.sin(A[4, 4].real))
        for i in range(d):
            for j in range(d):
                acc = complex(0.0, 0.0)
                for m in range(d):
                    acc += V0[i, m]*f[m]*V0[j, m].conjugate()
                out[b, i, j] = acc
        for i in range(d):
            lam[b, i] = A[i, i].real
        for i in range(1, d):
            key = lam[b, i]
            j = i - 1
            while j >= 0 and lam[b, j] > key:
                lam[b, j + 1] = lam[b, j]
                j -= 1
            lam[b, j + 1] = key
    return sev


_sinc_py = _sinc
r"""callable: The uncompiled :func:`_sinc`, kept reachable for the tests.

The module global ``_sinc`` is rebound to its compiled form below, because that
is how the kernels come to see it; this name keeps the Python original for a
direct comparison of the two.
"""


if HAVE_NUMBA:
    # fastmath is deliberately off.  The stability argument in the module
    # docstring rests on cancellations happening exactly as written -- the
    # sinc form of the divided differences, the vanishing of the second
    # divided difference's error with the gap, and for the Jacobi kernel the
    # rotation updates and Gram-Schmidt subtractions that keep the basis
    # unitary -- and fastmath licenses the compiler to reassociate precisely
    # those expressions.
    #
    # parallel is off too: magnus is routinely called from inside joblib
    # workers, and a kernel that opens its own thread pool inside each of them
    # oversubscribes the machine rather than speeding anything up.
    #
    # Lazy compilation with cache=True, not an eager signature: eager would
    # move the ~0.7 s compile to import time and charge it to every caller,
    # including those who never exponentiate a 2x2 or 3x3.
    _jit = nb.njit(cache=True, fastmath=False, nogil=True)
    # _sinc has to be compiled as well: nopython mode cannot call back into a
    # pure-Python function.  Rebinding the module global is what makes the
    # compiled form visible inside the kernels, because numba resolves global
    # names at the moment it compiles them -- so this line must come first.
    _sinc = _jit(_sinc)
    _ch2_kernel = _jit(_ch2_core)
    _ch3_kernel = _jit(_ch3_core)
    _jacobi_kernel = _jit(_jacobi_expm_core)
else:                                                   # pragma: no cover
    _ch2_kernel = _ch2_core
    _ch3_kernel = _ch3_core
    _jacobi_kernel = _jacobi_expm_core


def supports_dim(d: int) -> bool:
    r"""Returns whether dimension ``d`` has a compiled kernel.

    True for 2 through 5.  Dimensions 2 and 3 have Cayley-Hamilton closed
    forms; 4 and 5 go to the batched Jacobi eigensolver of
    :func:`_jacobi_expm_core` instead.  An earlier version of this docstring
    reasoned that a 4x4 or 5x5 Hermitian eigenproblem has no practical closed
    form and concluded that 4nu and 5nu stay on ``eigh``.  The premise stands;
    the conclusion did not follow from it, because the missing closed form was
    never what made those dimensions slow -- ``eigh``'s fixed per-matrix LAPACK
    overhead (~2.3 us on a 4x4, two thirds of a d=4 Magnus pass) was, and an
    iterative solver with no such overhead removes it without any closed form.
    This is the one place that decision is made.

    Parameters
    ----------
    d : int
        Matrix dimension.

    Returns
    -------
    bool
        Whether a kernel exists for that dimension.
    """
    return 2 <= d <= 5


def expm_herm_stack(K: np.ndarray) -> tuple:
    r"""Returns :math:`\exp(-iK)` and the eigenvalues of K, for Hermitian K.

    The drop-in replacement for the ``eigh`` half of
    ``magnus.magnus._expm_stack``: it returns the eigenvalues alongside the
    exponential because the caller needs them anyway, for the slab-width
    convergence warning, and this way there is no second spectral computation.

    Parameters
    ----------
    K : np.ndarray
        Hermitian matrix or stack of them, shape (..., d, d), with d 2 through
        5 (see :func:`supports_dim`).

    Returns
    -------
    tuple
        :math:`\exp(-iK)`, shape (..., d, d); the eigenvalues in ascending order, shape
        (..., d); and a float conditioning severity for the whole stack.  A severity above
        :data:`SEV_TOL` means at least one matrix could not be answered at full accuracy --
        too ill-conditioned for the 3x3 closed form (see ``_ch3_core``), or at the Jacobi
        sweep cap for a 4x4/5x5 (see ``_jacobi_expm_core``) -- and the caller should
        recompute with ``eigh``.

    Raises
    ------
    ValueError
        If ``d`` is not 2 through 5.  Without this the ``else`` below handed unsupported
        input to the 3x3 kernel, which returned no exception, an error of 2.4 against
        ``scipy.linalg.expm``, a unitarity violation of 11.3, and uninitialized memory in the
        fourth eigenvalue -- and segfaulted at d=1 by indexing ``K[i,2,1]`` with numba's
        bounds checking off.  :func:`supports_dim` is documented as the single place that
        decides which dimensions are handled, and now actually is.
    """
    d = K.shape[-1]
    if not supports_dim(d):
        raise ValueError(
            "Error in magnus: magnus.expmkernels.expm_herm_stack: no compiled kernel for dimension "
            + str(d) + "; only 2 through 5 are supported (see supports_dim). Callers "
            "should use numpy.linalg.eigh for larger dimensions, as "
            "magnus.magnus._expm_stack does.")
    # dtype=complex, not just ascontiguousarray: a real-valued Hermitian K is a
    # perfectly good input, and without the cast np.empty_like below would make a
    # real output buffer that the kernel then cannot store a complex result into.
    flat = np.ascontiguousarray(K, dtype=complex).reshape(-1, d, d)
    out = np.empty_like(flat)
    lam = np.empty((flat.shape[0], d), dtype=float)
    if d == 2:
        sev = _ch2_kernel(flat, out, lam)
    elif d == 3:
        sev = _ch3_kernel(flat, out, lam)
    else:
        sev = _jacobi_kernel(flat, out, lam)
    return out.reshape(K.shape), lam.reshape(K.shape[:-1]), float(sev)


__all__ = [
    'HAVE_NUMBA',
    'SEV_TOL',
    'supports_dim',
    'expm_herm_stack',
]
