r"""Argument checks shared by the public entry points (issue #160).

Every check here runs **once per call**, at the entry of a public function, and never inside a
quadrature, probe or slab loop.  Each raises :class:`TypeError` for a value of the wrong type
and :class:`ValueError` for a value of the right type but outside its range, and names both the
argument the caller passed and the public function the caller called.

The rules, stated once:

- **Reals** are Python or NumPy real scalars, 0-d arrays included; never ``bool``, never
  complex, always finite.
- **Integers** are ``int`` or ``np.integer``; never ``bool``, never a float such as 2.5.
- **Bools** are ``bool`` or ``np.bool_`` only.  Truth-testing is not a check: ``'False'`` is
  true.
- **Strings** come from a documented set.
- **Tolerances** ``rtol`` and ``atol`` are ``None`` or a finite real above zero.
- **Slab and point counts** are strictly positive integers.

.. versionadded:: 1.2.0

Routine listings
----------------

    * is_real_scalar - Whether a value is a single real number (not bool, not complex)
    * check_real - A finite real scalar, optionally signed or bounded
    * check_int - An integer, not bool, optionally bounded
    * check_bool - True or False only
    * check_choice - One of a set of values, compared by type as well
    * check_real_array - A real scalar or 1-D array, every entry finite
    * check_unit_fraction - A fraction in (0, 1], or an array of them
    * check_dict - None or a dict
    * check_slab_edges - A gap-free partition of the path into [start, end] pairs
    * check_hamiltonian_sample - A finite, square, Hermitian matrix (or stack)
    * check_refinement - The tolerance, slab and engine keywords, by one rule table
    * check_gl_order - Refuse an odd Magnus order on the Gauss-Legendre method
    * check_physics_params - A Hamiltonian builder's physics arguments, by name
    * validated - Decorator applying a rule table to calls from outside the package
    * r_real - Rule factory: a real number
    * r_int - Rule factory: an integer
    * r_bool - Rule: a bool
    * r_choice - Rule factory: one of a set
    * r_dict - Rule: None or a dict
    * r_callable - Rule: a callable
    * r_end_after - Rule factory: a real no smaller than another argument
    * r_hamiltonian_at - Rule factory: a callable Hamiltonian, sampled at another argument
    * r_hamiltonian - Rule: a Hamiltonian matrix
    * r_real_array - Rule factory: a real scalar or array
"""

import functools as _functools
import inspect as _inspect
import numbers
import sys as _sys

import numpy as np

__all__ = ['InputTypeError']


class InputTypeError(TypeError, ValueError):
    r"""An argument of the wrong type: a :class:`TypeError` that is also a :class:`ValueError`.

    Until 1.2.0 every input check raised :class:`ValueError`, wrong types included, and code
    written against that catches ``ValueError``.  The checks of issue #160 report a wrong type
    as a type error, and this class keeps both ``except`` clauses working.

    .. versionadded:: 1.2.0
    """

_BOOL_TYPES = (bool, np.bool_)
_INT_TYPES = (int, np.integer)


# The oscillation functions take arrays for energy and L only; said where it applies.
_SCAN_HINT = {'oscprob.': "  Only energy and L take arrays; to scan another parameter, make "
                          "one call per value."}


def _msg(where: str, text: str) -> str:
    # Imported here, on the failing path only: globaldefs imports magnus.magnus, which
    # imports this module.
    from magnus import globaldefs as gd
    return gd.ERROR_MSG_NO_COLOR + " " + where + ": " + text


def _show(x) -> str:
    r = repr(x)
    return r if len(r) <= 60 else r[:57] + '...'


def is_real_scalar(x) -> bool:
    r"""Whether ``x`` is a single real number: Python or NumPy, 0-d arrays included, not bool."""
    if isinstance(x, _BOOL_TYPES):
        return False
    if isinstance(x, (numbers.Real, np.floating, np.integer)):
        return True
    return isinstance(x, np.ndarray) and x.ndim == 0 and x.dtype.kind in 'iuf'


def check_real(name: str, x, where: str, *, positive: bool = False, nonnegative: bool = False,
               lo=None, hi=None, lo_open: bool = False, hi_open: bool = False,
               allow_none: bool = False, allow_inf: bool = False, what: str = None):
    r"""A single real number, finite, optionally signed or bounded; returned as ``float``.

    ``lo``/``hi`` are inclusive unless ``lo_open``/``hi_open``.  ``what`` replaces the default
    description of the expected value in the message.
    """
    if x is None:
        if allow_none:
            return None
        raise InputTypeError(_msg(where, name + " must be a real number, not None."))
    if type(x) is not float:
        if not is_real_scalar(x):
            if isinstance(x, _BOOL_TYPES):
                raise InputTypeError(_msg(where, name + " must be a real number, not a bool (" +
                                     _show(x) + ")."))
            if np.ndim(x) != 0:
                raise ValueError(_msg(where, name + " must be a single number, not an array "
                                      "or a list." + _SCAN_HINT.get(where[:8], "")))
            if isinstance(x, (complex, np.complexfloating)) or (
                    isinstance(x, np.ndarray) and x.dtype.kind == 'c'):
                raise InputTypeError(_msg(where, name + " must be real; got " + _show(x) + "."))
            raise InputTypeError(_msg(where, name + " must be a real number; got " +
                                 type(x).__name__ + " " + _show(x) + "."))
        x = float(x)
    if x != x or ((x == np.inf or x == -np.inf) and not allow_inf):
        raise ValueError(_msg(where, name + " must be finite; got " + repr(x) + "."))
    bad = ((positive and not x > 0.0) or (nonnegative and not x >= 0.0) or
           (lo is not None and (x <= lo if lo_open else x < lo)) or
           (hi is not None and (x >= hi if hi_open else x > hi)))
    if bad:
        if what is None:
            if positive:
                what = "positive"
            elif nonnegative:
                what = "non-negative"
            else:
                what = ("in " + ('(' if lo_open else '[') + str(lo if lo is not None else '-inf')
                        + ", " + str(hi if hi is not None else 'inf') + (')' if hi_open else ']'))
        raise ValueError(_msg(where, name + " must be " + what + "; got " + repr(x) + "."))
    return x


def check_int(name: str, x, where: str, *, lo=None, hi=None, allow_none: bool = False,
              what: str = None):
    r"""An integer (``int`` or ``np.integer``, not bool), bounded; returned as ``int``."""
    if x is None:
        if allow_none:
            return None
        raise InputTypeError(_msg(where, name + " must be an integer, not None."))
    if type(x) is not int:
        if isinstance(x, _BOOL_TYPES) or not isinstance(x, _INT_TYPES):
            if isinstance(x, np.ndarray) and x.ndim == 0 and x.dtype.kind in 'iu':
                x = int(x)
            elif not isinstance(x, (str, bytes)) and np.ndim(x) != 0:
                raise ValueError(_msg(where, name + " must be a single number, not an array "
                                      "or a list." + _SCAN_HINT.get(where[:8], "")))
            else:
                raise InputTypeError(_msg(where, name + " must be an integer; got " +
                                     type(x).__name__ + " " + _show(x) + "."))
        x = int(x)
    if (lo is not None and x < lo) or (hi is not None and x > hi):
        if what is None:
            if hi is None:
                what = ("a positive integer" if lo == 1 else "an integer >= " + str(lo))
            elif lo is None:
                what = "an integer <= " + str(hi)
            else:
                what = "an integer from " + str(lo) + " to " + str(hi)
        raise ValueError(_msg(where, name + " must be " + what + "; got " + str(x) + "."))
    return x


def check_bool(name: str, x, where: str, *, allow_none: bool = False):
    r"""``True`` or ``False`` (or ``np.bool_``); returned as ``bool``."""
    if x is True or x is False:
        return x
    if isinstance(x, np.bool_):
        return bool(x)
    if x is None and allow_none:
        return None
    raise InputTypeError(_msg(where, name + " must be True or False; got " + type(x).__name__ +
                         " " + _show(x) + ".  (Truth-testing is not used: 'False' and "
                         "'no' are non-empty strings, and would read as True.)"))


def check_choice(name: str, x, where: str, choices, *, allow_none: bool = False):
    r"""One of ``choices``, compared by equality and by type (so ``True`` is not ``1``)."""
    if x is None and allow_none:
        return None
    for c in choices:
        if type(x) is type(c) and x == c:
            return x
    raise ValueError(_msg(where, name + " must be one of " +
                          ", ".join(repr(c) for c in choices) + "; got " + _show(x) + "."))


def check_unit_fraction(name: str, x, where: str, *, allow_zero: bool = False):
    r"""A fraction in (0, 1] ([0, 1] with ``allow_zero``): a number or an array of them.

    For the helpers that run at every quadrature node (issue #160 §11): a plain float in range
    passes in one comparison; anything else, arrays included, gets the full check.  Returns
    ``x`` unchanged.
    """
    if type(x) is float and (0.0 <= x if allow_zero else 0.0 < x) and x <= 1.0:
        return x
    what = "in [0, 1]" if allow_zero else "in (0, 1]"
    if np.ndim(x) == 0:
        check_real(name, x, where, lo=0.0, lo_open=not allow_zero, hi=1.0, what=what)
        return x
    a = np.asarray(check_real_array(name, x, where, ndim=np.ndim(x)), dtype=float)
    bad = ~((a >= 0.0 if allow_zero else a > 0.0) & (a <= 1.0))
    if bad.any():
        i = int(np.argmax(bad.ravel()))
        raise ValueError(_msg(where, name + " must be " + what + "; entry " + str(i) + " is "
                              + repr(float(a.ravel()[i])) + "."))
    return x


def check_real_array(name: str, x, where: str, *, positive: bool = False,
                     nonnegative: bool = False, allow_scalar: bool = True, ndim: int = 1,
                     allow_empty: bool = False, allow_inf: bool = False):
    r"""A real scalar or a ``ndim``-dimensional real array, every entry checked.

    Returns the value unchanged when it is valid (arrays are not copied), so the caller's own
    handling of lists and arrays is untouched.
    """
    if np.ndim(x) == 0:
        if not allow_scalar:
            raise ValueError(_msg(where, name + " must be a " + str(ndim) + "-D array."))
        check_real(name, x, where, positive=positive, nonnegative=nonnegative,
                   allow_inf=allow_inf)
        return x
    if isinstance(x, np.ma.MaskedArray):
        raise InputTypeError(_msg(where, name + " must not be a masked array: the mask would be "
                             "ignored.  Pass the unmasked entries instead."))
    if not isinstance(x, (list, tuple, np.ndarray)):
        raise InputTypeError(_msg(where, name + " must be a number, a list or a NumPy array; got " +
                             type(x).__name__ + "."))
    a = np.asarray(x)
    if a.dtype.kind == 'b' or (a.dtype.kind == 'O' and any(isinstance(v, _BOOL_TYPES)
                                                          for v in a.ravel())):
        raise InputTypeError(_msg(where, name + " must hold real numbers, not bools."))
    if a.dtype.kind not in 'iuf':
        what = "complex numbers" if a.dtype.kind == 'c' else ("entries of type " + str(a.dtype))
        raise InputTypeError(_msg(where, name + " must hold real numbers; got " + what + "."))
    if a.ndim != ndim:
        raise ValueError(_msg(where, "if " + name + " is a list or NumPy array, it must be " +
                              str(ndim) + "-D; got " + str(a.ndim) + "-D."))
    if a.size == 0:
        if allow_empty:
            return x
        raise ValueError(_msg(where, name + " is empty."))
    af = a.astype(float, copy=False)
    fin = np.isfinite(af) if not allow_inf else ~np.isnan(af)
    if not fin.all():
        i = int(np.argmin(fin))
        raise ValueError(_msg(where, name + " must be finite; entry " + str(i) + " is " +
                              repr(float(af.ravel()[i])) + "."))
    if positive and not (af > 0.0).all():
        i = int(np.argmin(af > 0.0))
        raise ValueError(_msg(where, name + " must be positive; entry " + str(i) + " is " +
                              repr(float(af.ravel()[i])) + "."))
    if nonnegative and not (af >= 0.0).all():
        i = int(np.argmin(af >= 0.0))
        raise ValueError(_msg(where, name + " must be non-negative; entry " + str(i) + " is " +
                              repr(float(af.ravel()[i])) + "."))
    return x


def check_dict(name: str, x, where: str, *, allow_none: bool = True):
    r"""``None`` or a ``dict`` (filled in place by the callee)."""
    if (x is None and allow_none) or isinstance(x, dict):
        return x
    raise InputTypeError(_msg(where, name + " must be None or a dict; got " + type(x).__name__ +
                         "."))


def check_slab_edges(t_slab_edges, t_ini: float, t_fin: float, where: str,
                     name: str = 't_slab_edges', rel: float = 1e-12,
                     allow_zero_width: bool = False):
    r"""Slab edges as ``[[t0, t1], [t1, t2], ...]``: a gap-free partition of ``[t_ini, t_fin]``.

    Requires shape ``(n, 2)`` with ``n >= 1``, finite values, each pair increasing, each pair
    starting where the previous one ends, and the whole spanning ``[t_ini, t_fin]``, all to
    ``rel`` relative to the length of the interval.  Returns the edges as a float array.
    """
    try:
        e = np.asarray(t_slab_edges, dtype=float)
    except (TypeError, ValueError):
        raise InputTypeError(_msg(where, name + " must be a list of [start, end] pairs of real "
                             "numbers.")) from None
    if e.ndim != 2 or e.shape[1] != 2 or e.shape[0] < 1:
        raise ValueError(_msg(where, name + " must be a list of [start, end] pairs, "
                              "[[t0, t1], [t1, t2], ...]; got an array of shape " +
                              str(e.shape) + "."))
    if not np.isfinite(e).all():
        raise ValueError(_msg(where, name + " must be finite."))
    scale = max(abs(t_fin - t_ini), abs(t_ini), abs(t_fin), 1.0)
    tol = rel*scale
    widths = e[:, 1] - e[:, 0]
    bad = (widths < 0.0) if allow_zero_width else (widths <= 0.0)
    if bad.any():
        i = int(np.argmax(bad))
        raise ValueError(_msg(where, name + ": slab " + str(i) + " has end <= start (" +
                              repr(e[i].tolist()) + ")."))
    if e.shape[0] > 1:
        jumps = e[1:, 0] - e[:-1, 1]
        if (np.abs(jumps) > tol).any():
            i = int(np.argmax(np.abs(jumps) > tol))
            kind = "a gap" if jumps[i] > 0 else "an overlap"
            raise ValueError(_msg(where, name + ": there is " + kind + " between slab " +
                                  str(i) + " and slab " + str(i + 1) + " (" +
                                  repr(e[i].tolist()) + ", " + repr(e[i + 1].tolist()) +
                                  ").  The slabs must chain without gaps."))
    if abs(e[0, 0] - t_ini) > tol or abs(e[-1, 1] - t_fin) > tol:
        raise ValueError(_msg(where, name + " must span the whole path, from " + repr(t_ini) +
                              " to " + repr(t_fin) + "; it spans " + repr(float(e[0, 0])) +
                              " to " + repr(float(e[-1, 1])) + "."))
    return e


def check_hamiltonian_sample(name: str, H, where: str, *, at=None, hermitian: bool = True,
                             dim: int = None):
    r"""A sampled Hamiltonian: a finite square complex matrix (or a stack of them), Hermitian.

    Returns it as a ``complex128`` array.  ``at`` names the position sampled, for the message.
    """
    where_at = "" if at is None else " at " + str(at)
    try:
        A = np.asarray(H)
    except Exception:
        A = None
    if A is None or A.dtype.kind not in 'iufc':
        raise InputTypeError(_msg(where, name + " must be (or return) a square array of numbers; "
                             "got " + type(H).__name__ + where_at + "."))
    A = A.astype(np.complex128, copy=False)
    if A.ndim < 2 or A.shape[-1] != A.shape[-2] or A.shape[-1] == 0:
        raise ValueError(_msg(where, name + " must be (or return) a square matrix; got shape " +
                              str(A.shape) + where_at + "."))
    if dim is not None and A.shape[-1] != dim:
        raise ValueError(_msg(where, name + " must be " + str(dim) + "x" + str(dim) +
                              "; got " + str(A.shape[-2]) + "x" + str(A.shape[-1]) + where_at +
                              "."))
    if not np.isfinite(A).all():
        raise ValueError(_msg(where, name + " is not finite" + where_at + ": it contains " +
                              "NaN or inf."))
    if hermitian:
        norm = float(np.max(np.abs(A)))
        anti = float(np.max(np.abs(A - np.conj(np.swapaxes(A, -1, -2)))))
        if anti > 1e3*np.finfo(float).eps*max(norm, np.finfo(float).tiny):
            raise ValueError(_msg(where, name + " is not Hermitian" + where_at + ": its "
                                  "anti-Hermitian part is " + format(anti/max(norm, 1e-300),
                                                                      '.2g') +
                                  " of its largest entry.  A non-Hermitian Hamiltonian gives "
                                  "a non-unitary evolution and probabilities that do not sum "
                                  "to 1."))
    return A


# The rules for the refinement and engine keywords, shared by every scenario function and by
# osc_prob_energy_baseline.  Each maps a keyword to a check applied when the keyword is present
# (and not None where None is allowed).  Kept as a table so that the wrappers, the scenario
# functions and the CLI all refuse the same values with the same words.
def _tol(name, x, where):
    return check_real(name, x, where, positive=True, allow_none=True,
                      what="None or a finite number above 0")


def _count(lo=1, allow_none=True):
    return lambda name, x, where: check_int(name, x, where, lo=lo, allow_none=allow_none)


def _growth(name, x, where):
    return check_real(name, x, where, lo=1.0, lo_open=True,
                      what="a finite number above 1 (a factor of 1 never refines)")


def _order(name, x, where):
    return check_int(name, x, where, lo=1, hi=10, what="an integer from 1 to 10")


def _n_jobs(name, x, where):
    if x is None:
        return None
    x = check_int(name, x, where)
    if x == -1 or x >= 1:
        return x
    raise ValueError(_msg(where, name + " must be -1 (all cores) or a positive integer; got " +
                          str(x) + "."))


def _flag(name, x, where):
    return check_bool(name, x, where)


def _cumulative(name, x, where):
    if x is True or x is False or isinstance(x, np.bool_) or (isinstance(x, str) and
                                                             x == 'auto'):
        return x
    raise ValueError(_msg(where, name + " must be True, False, or 'auto'; got " + _show(x) +
                          "."))


REFINEMENT_RULES = {
    'rtol': _tol,
    'atol': _tol,
    'n_slabs': _count(1),
    'min_n_slabs': _count(1),
    'max_n_slabs': _count(1),
    'n_tpts_per_slab': _count(2),
    'min_n_tpts_per_slab': _count(2),
    'max_n_tpts_per_slab': _count(2),
    'max_num_loops': _count(1, allow_none=False),
    'growth_factor_n_slabs': _growth,
    'growth_factor_n_tpts_per_slab': _growth,
    'magnus_exp_order': _order,
    'n_jobs': _n_jobs,
    'integration_method': lambda name, x, where: check_choice(
        name, x, where, ('gl', 'trapezoid', 'simpson')),
    'strict_convergence': _flag,
    'cumulative': _cumulative,
    'new_recursion_limit': _count(1),
    'verbose': lambda name, x, where: check_int(name, x, where, lo=0, hi=2),
    'strategy_info': lambda name, x, where: check_dict(name, x, where),
    'convergence_info': lambda name, x, where: check_dict(name, x, where),
}


def check_refinement(where: str, values: dict) -> None:
    r"""Apply :data:`REFINEMENT_RULES` to the entries of ``values`` that it names.

    Also refuses the combinations no single rule sees: a floor above its ceiling.
    """
    for key, x in values.items():
        rule = REFINEMENT_RULES.get(key)
        if rule is not None:
            rule(key, x, where)
    # A floor above its ceiling is a contradiction.  n_slabs above max_n_slabs is not: it is
    # clipped to the cap, with ToleranceNotAchievedWarning, by design.
    if 'magnus_exp_order' in values:
        check_gl_order(values['magnus_exp_order'], values.get('integration_method', 'gl'), where)
    for lo_key, hi_key in (('min_n_slabs', 'max_n_slabs'),
                           ('min_n_tpts_per_slab', 'max_n_tpts_per_slab')):
        lo, hi = values.get(lo_key), values.get(hi_key)
        if lo is not None and hi is not None and lo > hi:
            raise ValueError(_msg(where, lo_key + " (" + str(lo) + ") must be <= " + hi_key +
                                  " (" + str(hi) + ")."))


def check_gl_order(order, integration_method, where: str) -> None:
    r"""Refuses an odd Magnus order on the Gauss-Legendre method (issue #160 §5).

    Each Gauss-Legendre scheme integrates to an even order: an odd request ran the scheme of the
    next even order and returned its result bit for bit, so ``magnus_exp_order=3`` was order 4
    under another name.  Refused, naming the order that was being computed.
    """
    if integration_method in (None, 'gl') and isinstance(order, (int, np.integer)) and \
            not isinstance(order, bool) and order % 2 == 1:
        raise ValueError(_msg(where, "magnus_exp_order=" + str(order) + " is odd, and the "
                              "Gauss-Legendre schemes of integration_method='gl' have even "
                              "orders only: this ran order " + str(order + 1) + ", bit for bit.  "
                              "Pass magnus_exp_order=" + str(order + 1) + ", or "
                              "integration_method='trapezoid' or 'simpson' for an odd order."))


# ---------------------------------------------------------------------------------------------
# Rule-table checks for the public helper functions (avgprob, adiabatic, magnus, earth, ...).
#
# The engines call several of these helpers themselves, once per energy point, with arguments
# they have already validated.  The checks below therefore run only when the caller is outside
# the package: one frame lookup decides, so an engine pays about a microsecond per call and a
# user gets every argument checked.
# ---------------------------------------------------------------------------------------------



def _called_from_inside(depth: int) -> bool:
    return _sys._getframe(depth).f_globals.get('__name__', '').startswith('magnus.')


def validated(spec: dict):
    r"""Decorate a public helper with argument rules, applied to calls from outside magnus.

    ``spec`` maps an argument name to ``rule(name, value, where, args)``, where ``args`` holds
    every argument of the call (defaults applied), so a rule can compare two of them.  A rule
    returns nothing; it raises to refuse.
    """
    def deco(func):
        sig = _inspect.signature(func)
        where = func.__module__.replace('magnus.', '', 1) + '.' + func.__name__
        # Keywords collected by a **kwargs parameter are checked by the same rules.
        var_kw = next((n for n, p in sig.parameters.items() if p.kind is p.VAR_KEYWORD), None)

        @_functools.wraps(func)
        def wrapper(*args, **kwargs):
            if not _called_from_inside(2):
                try:
                    bound = sig.bind(*args, **kwargs)
                except TypeError:
                    return func(*args, **kwargs)
                bound.apply_defaults()
                a = bound.arguments
                if var_kw is not None and a.get(var_kw):
                    a = dict(a, **a[var_kw])
                for name, rule in spec.items():
                    if name in a:
                        rule(name, a[name], where, a)
            return func(*args, **kwargs)
        return wrapper
    return deco


def r_real(**kw):
    return lambda name, x, where, a: check_real(name, x, where, **kw)


def r_int(**kw):
    return lambda name, x, where, a: check_int(name, x, where, **kw)


def r_bool(name, x, where, a):
    check_bool(name, x, where)


def r_choice(choices, allow_none=False):
    return lambda name, x, where, a: check_choice(name, x, where, choices,
                                                  allow_none=allow_none)


def r_dict(name, x, where, a):
    check_dict(name, x, where)


def r_callable(name, x, where, a):
    if not callable(x):
        raise InputTypeError(_msg(where, name + " must be callable; got " + type(x).__name__ +
                                  "."))


def r_end_after(start: str):
    r"""A finite real that is >= the argument named ``start``."""
    def rule(name, x, where, a):
        x = check_real(name, x, where)
        s = check_real(start, a[start], where)
        if x < s:
            raise ValueError(_msg(where, name + " (" + repr(x) + ") must be >= " + start +
                                  " (" + repr(s) + "); the interval runs from " + start +
                                  " to " + name + "."))
    return rule


def r_hamiltonian_at(position: str):
    r"""A callable Hamiltonian whose sample at the argument ``position`` is sound."""
    def rule(name, x, where, a):
        r_callable(name, x, where, a)
        at = a[position]
        check_hamiltonian_sample(name, x(at), where, at=position + ' = ' + format(at, '.6g'))
    return rule


def r_hamiltonian(name, x, where, a):
    check_hamiltonian_sample(name, x, where)


def r_real_array(**kw):
    return lambda name, x, where, a: (None if x is None else
                                      check_real_array(name, x, where, **kw))


# Diagonal NSI couplings are real by hermiticity; the off-diagonal ones may be complex.
NSI_DIAGONAL = frozenset(('eps_aa', 'eps_ee', 'eps_mm', 'eps_tt', 'eps_ss', 'eps_s1s1',
                          'eps_s2s2'))
_BUILDER_FLAGS = frozenset(('nubar', 'compute_matrix_multiplication'))
_INF = float('inf')


def check_physics_params(where: str, values: dict) -> None:
    r"""Check the physics arguments of a Hamiltonian builder, by name.

    Energies and ``Lambda`` are positive; ``n_liv`` is an integer >= 0; off-diagonal NSI
    couplings may be complex; every other sine, phase, splitting, coupling and LIV
    coefficient is a finite real; the two flags are bools.  A finite plain float passes in
    one comparison, so a builder called once per quadrature node inside a user Hamiltonian
    pays well under a microsecond (issue #160 §11).
    """
    for k, x in values.items():
        if k in _BUILDER_FLAGS:
            if x is not True and x is not False:
                check_bool(k, x, where)
            continue
        if k == 'energy' or k == 'Lambda':
            if not (type(x) is float and 0.0 < x < _INF):
                if k == 'energy' and np.ndim(x) != 0:
                    # The builders broadcast over an array of energies (issue #155 §2).
                    check_real_array(k, x, where, positive=True, ndim=np.ndim(x))
                else:
                    check_real(k, x, where, positive=True)
            continue
        if k == 'n_liv':
            check_int(k, x, where, lo=0, what="an integer >= 0 (the operator dimension "
                      "minus 3)")
            continue
        if k.startswith('eps_') and k not in NSI_DIAGONAL:
            if type(x) is float and -_INF < x < _INF:
                continue
            if isinstance(x, (complex, np.complexfloating)) and not isinstance(x, bool):
                if not (np.isfinite(x.real) and np.isfinite(x.imag)):
                    raise ValueError(_msg(where, k + " must be finite; got " + repr(x) + "."))
                continue
            check_real(k, x, where)
            continue
        if k[:1] in ('s', 'd', 'D', 'b') or k.startswith('eps_'):
            if type(x) is float and -_INF < x < _INF:
                continue
            if isinstance(x, (complex, np.complexfloating)) and x.imag == 0:
                continue
            if isinstance(x, (complex, np.complexfloating)):
                raise InputTypeError(_msg(where, k + " must be real" +
                    (" (a diagonal NSI coupling is real by hermiticity)"
                     if k in NSI_DIAGONAL else "") + "; got " + repr(x) + "."))
            check_real(k, x, where)
