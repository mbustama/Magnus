r"""Every Hamiltonian builder broadcasts over its one varying argument (issue #155 §2).

``energy``, ``VCC`` or ``l`` given as an array of shape ``s`` returns a stack of shape
``s + (d, d)`` whose entries equal the scalar calls exactly; given as a number, the builder
returns one ``(d, d)`` matrix, as before.  Before this, the vacuum, NSI and LIV builders either
raised a broadcasting error or, when the array length equalled ``d``, returned one wrong matrix
without a word.
"""

import inspect

import numpy as np
import pytest

import magnus.hamiltonians as hams
from magnus._validate import InputTypeError

NSI_DIAGONAL = {'eps_aa', 'eps_ee', 'eps_mm', 'eps_tt', 'eps_ss', 'eps_s1s1', 'eps_s2s2'}


def _value(name, k):
    """A valid, generic value for builder argument ``name``; ``k`` separates siblings."""
    if name in NSI_DIAGONAL:
        return 0.1 + 0.01*k
    if name.startswith('eps_'):
        return 0.05 + 0.02j*k
    if name in ('sth', 'sxi') or name[0] == 's':
        return 0.2 + 0.03*k
    if name[0] == 'd':
        return 0.7 + 0.1*k
    if name[0] == 'D':
        return 7.0e-5*(k + 1)**2
    if name[0] == 'b':
        return 1.0e-23*(k + 1)
    if name == 'Lambda':
        return 1.0e9
    if name == 'n_liv':
        return 1
    raise KeyError(name)


def _builders():
    for name in sorted(dir(hams)):
        m = __import__('re').fullmatch(r'hamiltonian_([2-5])nu_(vacuum|vacuum_td|'
                                       r'vacuum_energy_independent_td|nsi|nsi_td|liv|matter|'
                                       r'matter_td)', name)
        if m:
            yield name


BUILDERS = list(_builders())
VCC_FUNC = lambda l: 1.0e-13*np.exp(-np.asarray(l, dtype=float)/3.0e12)


def _call(name, varying):
    """Call builder ``name`` with its varying argument set to ``varying``."""
    f = getattr(hams, name)
    kw = {}
    for k, p in enumerate(inspect.signature(f).parameters.values()):
        if p.default is not inspect.Parameter.empty:
            continue
        if p.name in ('energy', 'VCC', 'l'):
            continue
        if p.name == 'VCC_func':
            kw[p.name] = VCC_FUNC
        else:
            kw[p.name] = _value(p.name, k)
    params = inspect.signature(f).parameters
    for arg in ('energy', 'VCC', 'l'):
        if arg in params:
            if arg == 'energy' and 'l' in params:      # *_vacuum_td: vary l, fix the energy
                kw['energy'] = 2.0e9
                continue
            kw[arg] = varying[arg]
    return f(**kw)


VARYING_ARRAY = dict(energy=np.array([1.0e9, 2.0e9, 3.5e9]),
                     VCC=np.array([0.0, 1.0e-13, 3.0e-13]),
                     l=np.array([0.0, 1.0e12, 5.0e12]))


def _entry(i):
    return {k: float(v[i]) for k, v in VARYING_ARRAY.items()}


def test_every_varying_builder_is_covered():
    """The flavor counts and kinds this file claims, so a builder added later is not missed."""
    assert len(BUILDERS) == 4*8


@pytest.mark.parametrize('name', BUILDERS)
def test_array_returns_a_stack_equal_to_the_scalar_calls(name):
    stack = np.asarray(_call(name, VARYING_ARRAY))
    singles = [np.asarray(_call(name, _entry(i))) for i in range(3)]
    d = singles[0].shape[-1]
    assert singles[0].shape == (d, d)
    assert stack.shape == (3, d, d)
    for i in range(3):
        assert np.array_equal(stack[i], singles[i])


@pytest.mark.parametrize('name', BUILDERS)
def test_length_equal_to_the_dimension_is_a_stack_not_one_mixed_matrix(name):
    """Two energies at two flavors used to give one 2x2 matrix mixing both columns."""
    d = np.asarray(_call(name, _entry(0))).shape[-1]
    varying = {k: np.resize(v, d) for k, v in VARYING_ARRAY.items()}
    assert np.asarray(_call(name, varying)).shape == (d, d, d)


@pytest.mark.parametrize('name', [n for n in BUILDERS if n.endswith(('_vacuum', '_liv'))])
def test_array_energy_is_checked_entry_by_entry(name):
    with pytest.raises(ValueError, match=r'energy must be positive; entry 1'):
        _call(name, dict(VARYING_ARRAY, energy=np.array([1.0e9, -1.0, 2.0e9])))
    with pytest.raises(ValueError, match=r'energy must be finite; entry 2'):
        _call(name, dict(VARYING_ARRAY, energy=np.array([1.0e9, 2.0e9, np.nan])))
    with pytest.raises(InputTypeError, match='energy'):
        _call(name, dict(VARYING_ARRAY, energy=np.array([1.0e9, 2.0e9j])))


def test_pseudo_dirac_vacuum_broadcasts_and_checks_the_energy():
    U = hams.pmns_mixing_matrix(0.55, 0.75, 0.15, 3.7)
    m2 = [0.0, 7.4e-5, 2.5e-3]
    pairs = {1: 1.0e-12}
    E = np.array([1.0e9, 2.0e9, 3.0e9])
    stack = hams.hamiltonian_pseudo_dirac_vacuum(E, U, m2, pairs)
    assert stack.shape == (3, 4, 4)
    for i in range(3):
        assert np.array_equal(stack[i], hams.hamiltonian_pseudo_dirac_vacuum(
            float(E[i]), U, m2, pairs))
    with pytest.raises(ValueError, match='energy must be positive'):
        hams.hamiltonian_pseudo_dirac_vacuum(0.0, U, m2, pairs)


def test_an_h_func_built_from_the_nsi_builder_takes_the_vectorized_path():
    """The case that motivated #155 §2: hv/E + hamiltonian_3nu_nsi(VCC(l), ...)."""
    import warnings

    import magnus.globaldefs as gd
    import magnus.oscprob as oscprob
    from magnus.magnus import ScalarHamiltonianWarning
    osc = gd.load_nufit_params('NuFIT 6.1')
    hv = hams.hamiltonian_3nu_vacuum_energy_independent(**osc)
    L = 1000.0*gd.UNIT_KM
    vcc = lambda l: 1.0e-13*(1.0 + 0.3*np.sin(np.asarray(l, dtype=float)/L))

    def H(E, l):
        return hv/E + hams.hamiltonian_3nu_matter(vcc(l)) + hams.hamiltonian_3nu_nsi(
            vcc(l), 0.1, 0.05, 0.0, 0.0, 0.0, 0.0)

    with warnings.catch_warnings():
        warnings.simplefilter('error', ScalarHamiltonianWarning)
        P = oscprob.osc_prob_energy_baseline(H, 2.0*gd.UNIT_GEV, L, n_slabs=16, rtol=None,
                                            atol=None)
    assert np.allclose(np.sum(P, axis=-1), 1.0)
