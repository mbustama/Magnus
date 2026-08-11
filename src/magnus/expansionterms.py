# -*- coding: utf-8 -*-
# SPDX-License-Identifier: GPL-3.0-only
# Copyright (C) 2026 Mauricio Bustamante
r"""expansionterms.py

Generate the terms of the Magnus expansion symbolically, to any order.

The numerical core in :mod:`magnus.magnus` evaluates the Magnus expansion with the
coefficients of each commutator group written out explicitly in Python, which is fast but
fixes the highest order at whatever was typed in.  This module derives those same terms
from the recursion itself, in exact rational arithmetic, for any order -- so the
hard-coded ones can be *checked* rather than trusted, and so an order beyond the
implemented ceiling can still be inspected on paper.

The recursion is the standard Bernoulli-number one [1]_ (in the :math:`B_1 = -1/2`
convention):

.. math::

   \Omega_1(t) &= \int_0^t A(s)\, ds \\
   \Omega_n(t) &= \sum_{j=1}^{n-1} \frac{B_j}{j!} \int_0^t S_n^{(j)}(s)\, ds ,

where the :math:`S_n^{(j)}` are themselves defined recursively,

.. math::

   S_n^{(1)} &= [\Omega_{n-1}, A] \\
   S_n^{(j)} &= \sum_{m=1}^{n-j} [\Omega_m, S_{n-m}^{(j-1)}] ,
   \qquad 2 \leq j \leq n-1 ,

so that every term of :math:`\Omega_n` is a nested commutator of lower-order
:math:`\Omega_m` with :math:`A`, carrying a rational coefficient.  Because
:math:`B_j = 0` for every odd :math:`j \geq 3`, whole groups drop out: only
:math:`j = 1, 2, 4, 6, \ldots` contribute.

This is the same form the numerical core implements, deliberately.  Writing
:math:`\Omega_n` instead as time-ordered multiple integrals of :math:`A` alone is the
other common presentation, but it does not correspond to anything the code evaluates, so
it would be of no use for checking the implementation.

The number of terms grows quickly -- 1 at order 1, 26 at order 6, 211 at order 8, 1918 at
order 10 -- which is why the implemented ceiling is a deliberate choice rather than an
oversight.  See :doc:`/expansion_terms` for the derivation, the expansion printed out,
and worked examples.

References
----------
.. [1] S. Blanes, F. Casas, J. A. Oteo & J. Ros, "The Magnus expansion
   and some of its applications", Phys. Rep. 470, 151 (2009).

Routine listings
----------------

    * bernoulli - Bernoulli number B_n as an exact Fraction
    * bernoulli_factor - The coefficient B_j / j! of a commutator group
    * omega_terms - Terms of Omega_n, as (coefficient, nested-commutator) pairs
    * magnus_terms - omega_terms for every order up to the one requested
    * format_term - One term as a readable string
    * print_magnus_terms - The expansion, printed order by order
    * count_terms - Number of terms in Omega_n, without building them
"""

__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"


from fractions import Fraction
from functools import lru_cache
from math import comb, factorial
from typing import Dict, List, Tuple, Union

# A "word" is a nested commutator built from two kinds of leaf:
#     'A'          -- the matrix function A(t) itself
#     ('Om', m)    -- the m-th Magnus term, Omega_m
# and one node type:
#     ('c', X, Y)  -- the commutator [X, Y]
# Words are plain tuples so they are hashable, and so like terms can be collected by
# using the word as a dictionary key.
Word = Union[str, Tuple]

# A term is a rational coefficient paired with a word.
Term = Tuple[Fraction, Word]


@lru_cache(maxsize=None)
def bernoulli(n: int) -> Fraction:
    r"""Returns the Bernoulli number :math:`B_n` as an exact :class:`fractions.Fraction`.

    Uses the :math:`B_1 = -1/2` convention, which is the one the Magnus recursion above is
    written in, and computes from the defining recursion

    .. math::

       B_m = -\frac{1}{m+1} \sum_{j=0}^{m-1} \binom{m+1}{j} B_j ,

    in exact rational arithmetic, so the result is not subject to rounding at any order.
    Every odd :math:`B_n` with :math:`n \geq 3` comes out exactly zero.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    n : int
        Index of the Bernoulli number; must be >= 0.

    Returns
    -------
    fractions.Fraction
        :math:`B_n`.

    Raises
    ------
    ValueError
        If ``n`` is negative.

    Examples
    --------
    .. jupyter-execute::

        import magnus.expansionterms as et

        [str(et.bernoulli(n)) for n in range(9)]
    """
    if n < 0:
        raise ValueError("magnus.expansionterms.bernoulli: n must be >= 0, not "
                         + str(n) + "Error in magnus: .")
    if n == 0:
        return Fraction(1)
    total = Fraction(0)
    for j in range(n):
        total += comb(n + 1, j) * bernoulli(j)
    return -total / (n + 1)


def bernoulli_factor(j: int) -> Fraction:
    r"""Returns :math:`B_j / j!`, the coefficient multiplying the :math:`j`-th commutator
    group in the Magnus recursion.

    These are the numbers the numerical core hard-codes: :math:`B_1/1! = -1/2`,
    :math:`B_2/2! = 1/12`, :math:`B_4/4! = -1/720`, :math:`B_6/6! = 1/30240`.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    j : int
        Index of the commutator group; must be >= 0.

    Returns
    -------
    fractions.Fraction
        :math:`B_j / j!`.

    Examples
    --------
    .. jupyter-execute::

        import magnus.expansionterms as et

        {j: str(et.bernoulli_factor(j)) for j in (1, 2, 4, 6)}
    """
    return bernoulli(j) / factorial(j)


def _collect(terms: List[Term]) -> List[Term]:
    """Sums the coefficients of repeated words and drops the ones that cancel.

    Two terms with the same nested-commutator structure are the same term, so they must be
    combined before the count means anything; without this, the recursion reports terms
    that are not there.
    """
    merged: Dict[Word, Fraction] = {}
    for coeff, word in terms:
        merged[word] = merged.get(word, Fraction(0)) + coeff
    return [(c, w) for w, c in merged.items() if c != 0]


@lru_cache(maxsize=None)
def _s_terms(n: int, j: int) -> Tuple[Term, ...]:
    r"""Terms of :math:`S_n^{(j)}`, memoized.

    Returns a tuple (rather than a list) so it can be cached.
    """
    if j == 1:
        return ((Fraction(1), ('c', ('Om', n - 1), 'A')),)
    out: List[Term] = []
    for m in range(1, n - j + 1):
        for coeff, word in _s_terms(n - m, j - 1):
            out.append((coeff, ('c', ('Om', m), word)))
    return tuple(_collect(out))


@lru_cache(maxsize=None)
def omega_terms(order: int) -> Tuple[Term, ...]:
    r"""Returns the terms of :math:`\Omega_n` for ``n = order``.

    Each term is a ``(coefficient, word)`` pair: an exact
    :class:`fractions.Fraction` and a nested commutator built from ``'A'`` and
    ``('Om', m)`` leaves via ``('c', X, Y)`` nodes.  The whole of :math:`\Omega_n` is the
    integral of the sum of these terms, as in the module docstring.

    Odd Bernoulli numbers above :math:`B_1` vanish, so the :math:`j = 3, 5, 7, \ldots`
    groups contribute nothing and are skipped rather than generated and discarded.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    order : int
        Order :math:`n` of the term; must be >= 1.

    Returns
    -------
    tuple of (fractions.Fraction, tuple)
        The terms of :math:`\Omega_n`.  For ``order=1`` this is the single term
        :math:`(1, A)`, the integrand of :math:`\Omega_1 = \int A`.

    Raises
    ------
    ValueError
        If ``order`` is less than 1.

    Examples
    --------
    .. jupyter-execute::

        import magnus.expansionterms as et

        for coeff, word in et.omega_terms(3):
            print(f"{str(coeff):>6s}  {et.format_term((coeff, word), with_coeff=False)}")
    """
    if order < 1:
        raise ValueError("magnus.expansionterms.omega_terms: order must be >= 1, not "
                         + str(order) + "Error in magnus: .")
    if order == 1:
        return ((Fraction(1), 'A'),)
    out: List[Term] = []
    for j in range(1, order):
        factor = bernoulli_factor(j)
        if factor == 0:
            continue
        for coeff, word in _s_terms(order, j):
            out.append((factor * coeff, word))
    return tuple(_collect(out))


def magnus_terms(max_order: int) -> Dict[int, Tuple[Term, ...]]:
    r"""Returns :func:`omega_terms` for every order from 1 up to ``max_order``.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    max_order : int
        Highest order to generate; must be >= 1.

    Returns
    -------
    dict
        Maps each order :math:`n` to the terms of :math:`\Omega_n`.

    Examples
    --------
    .. jupyter-execute::

        import magnus.expansionterms as et

        {n: len(terms) for n, terms in et.magnus_terms(8).items()}
    """
    if max_order < 1:
        raise ValueError("magnus.expansionterms.magnus_terms: max_order must be >= 1, not "
                         + str(max_order) + "Error in magnus: .")
    return {n: omega_terms(n) for n in range(1, max_order + 1)}


def count_terms(order: int) -> int:
    r"""Returns the number of terms in :math:`\Omega_n`.

    Builds the terms and counts them, so the count reflects the collection of like terms
    rather than the raw size of the recursion.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    order : int
        Order :math:`n`; must be >= 1.

    Returns
    -------
    int
        Number of distinct commutator terms in :math:`\Omega_n`.

    Examples
    --------
    .. jupyter-execute::

        import magnus.expansionterms as et

        {n: et.count_terms(n) for n in range(1, 11)}
    """
    return len(omega_terms(order))


def _format_word(word: Word) -> str:
    """Renders a word as nested square-bracket commutators."""
    if word == 'A':
        return 'A'
    if isinstance(word, tuple) and word[0] == 'Om':
        return 'Om_' + str(word[1])
    _, left, right = word
    return '[' + _format_word(left) + ', ' + _format_word(right) + ']'


def format_term(term: Term, with_coeff: bool = True) -> str:
    r"""Renders one term as a readable string.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    term : (fractions.Fraction, tuple)
        A term, as returned by :func:`omega_terms`.
    with_coeff : bool, optional
        If True (default), prefix the commutator with its coefficient.

    Returns
    -------
    str
        For example, ``'+1/12 [Om_1, [Om_1, A]]'``.

    Examples
    --------
    .. jupyter-execute::

        import magnus.expansionterms as et

        [et.format_term(t) for t in et.omega_terms(3)]
    """
    coeff, word = term
    body = _format_word(word)
    if not with_coeff:
        return body
    sign = '-' if coeff < 0 else '+'
    return sign + str(abs(coeff)) + ' ' + body


def print_magnus_terms(max_order: int, file=None) -> None:
    r"""Prints the Magnus expansion, order by order, up to ``max_order``.

    Each order is printed as the integrand of :math:`\Omega_n`, one term per line, with
    its exact rational coefficient.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    max_order : int
        Highest order to print; must be >= 1.
    file : file-like, optional
        Destination, forwarded to :func:`print`. Default: standard output.

    Returns
    -------
    None

    Examples
    --------
    .. jupyter-execute::

        import magnus.expansionterms as et

        et.print_magnus_terms(4)
    """
    for n, terms in magnus_terms(max_order).items():
        header = "Omega_" + str(n) + "  (" + str(len(terms)) + " term"
        header += "" if len(terms) == 1 else "s"
        header += ")"
        print(header, file=file)
        if n == 1:
            print("    int A", file=file)
        else:
            for term in terms:
                print("    int " + format_term(term), file=file)
        print(file=file)


__all__ = [
    # The two type aliases below appear in the signatures of everything else
    # here, so they have to be documented or those signatures render as dead
    # references.
    'Word',
    'Term',
    'bernoulli',
    'bernoulli_factor',
    'omega_terms',
    'magnus_terms',
    'count_terms',
    'format_term',
    'print_magnus_terms',
]
