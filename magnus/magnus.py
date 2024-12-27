import numpy as np
import scipy as sp
# from scipy.integrate import quad
# from scipy.linalg import expm
# from scipy.special import factorial
from numba import njit, prange
# from joblib import Parallel, delayed


# Bernoulli numbers
B = {
    0:  1,
    1:  -0.5, # Use the negative-number convention
    2:  1./6.,
    3:  0,
    4:  -1./30.,
    5:  0,
    6:  1./42.,
    7:  0,
    8:  -1./30.,
    9:  0,
    10: 5./66.,
    11: 0,
    12: -691./2730.,
    13: 0,
    14: 7./6.,
    15: 0,
    16: -3617./510.,
    17: 0,
    18: 43867./798.,
    19: 0,
    20: -174611./330.
}

# Multiplicative factors
f1 = 1.0/12.0
f2 = -1.0/720.0


# Function to compute the Magnus expansion terms
def compute_magnus_terms(A, t0, t1, n_tpts=50, order=2, integration_method='trapezoid'):
    """
    Compute the Magnus expansion terms up to a given order over the range [t0, t1].
    """
    # Validate input
    valid_integration_methods = ['trapezoid', 'simpson']
    if not (integration_method in valid_integration_methods):
        raise ValueError("compute_magnus_terms: integration_method must be one of "+ \
            str(valid_integration_methods)+".")

    nA = A(t0).shape[0]
    zero_matrix = np.zeros((nA, nA), dtype=complex)

    # Precompute time points and weights
    if t0 > 0.0:
        times = np.logspace(np.log10(t0), np.log10(t1), n_tpts)
    else:
        times = np.linspace(t0, t1, n_tpts)
    delta = (t1 - t0) / (n_tpts - 1)

    def commutator(X, Y):
        return X @ Y - Y @ X

    def integral_cumulative_simpson(times, matrices):
        # We need to write our custom routine to compute cumulative matrix integrals because the
        # scipy routine `integrate.cumulative_simpson` does not handle complex numbers; it casts
        # them into real numbers.  This is not a problem of the `integrate.cumulative_trapezoid`, 
        # which we do use.
        result = []
        for i in range(len(times)):
            if (i == 0): # Integral from t0 to t0
                result.append(zero_matrix)
            elif (i == 1): # Integral from t0 to t1
                result.append(0.5*(times[1]-times[0])*sum(matrices[:2]))
            else:
                result.append(sp.integrate.simpson(matrices[0:i+1], x=times[0:i+1], axis=0))
        return np.array(result)

    magnus_terms = []

    # Precompute the A(t) terms
    At = np.array([A(t) for t in times])

    # Precompute the Omega_1(t) terms, integrating from t0 to t = t0, ..., t1
    if (integration_method == 'trapezoid'):
        o1t = sp.integrate.cumulative_trapezoid(At, x=times, axis=0, initial=0)
    elif (integration_method == 'simpson'):
        # o1t = sp.integrate.cumulative_simpson(At, x=times, axis=0, initial=0)
        o1t = integral_cumulative_simpson(times, At)
    magnus_terms.append(o1t[-1]) # Integral from t0 to t1

    # Precompute the Omega_2(t) terms, integrating from t0 to t = t0, ..., t1
    if order >= 2:
        o2t_integrand = -0.5 * np.stack([commutator(o1t[i], At[i]) for i in range(n_tpts)], axis=0)
        if (integration_method == 'trapezoid'):
            o2t = sp.integrate.cumulative_trapezoid(o2t_integrand, x=times, axis=0, initial=0)
        elif (integration_method == 'simpson'):
            # o2t = sp.integrate.cumulative_simpson(o2t_integrand, x=times, axis=0, initial=0)
            o2t = integral_cumulative_simpson(times, o2t_integrand)
        magnus_terms.append(o2t[-1])

    # Precompute the Omega_3(t) terms, integrating from t0 to t = t0, ..., t1
    if order >= 3:
        t1 = -0.5 * np.stack([commutator(o2t[i], At[i]) for i in range(n_tpts)], axis=0)
        t2 = f1 * np.stack(
            [commutator(o1t[i], commutator(o1t[i], At[i])) for i in range(n_tpts)], axis=0
        )
        o3t_integrand = t1 + t2
        if (integration_method == 'trapezoid'):
            o3t = sp.integrate.cumulative_trapezoid(o3t_integrand, x=times, axis=0, initial=0)
        elif (integration_method == 'simpson'):
            # o3t = sp.integrate.cumulative_simpson(o3t_integrand, x=times, axis=0, initial=0)
            o3t = integral_cumulative_simpson(times, o3t_integrand)
        # The term in the Magnus expansion is the last one, which goes from t0 to t1
        magnus_terms.append(o3t[-1])
        # ------------------------------------------------------------------------------------------
        # Deprecated (slower):
        # t1 = -0.5 * np.array([commutator(o2t[i], At[i]) for i in range(n_tpts)])
        # t2 = (1.0/12.0) * np.array([commutator(o1t[i], commutator(o1t[i], At[i])) for i in range(n_tpts)])
        # o3t_integrand = t1+t2
        # o3t = [matrix_integral(times[0:i+1], o3t_integrand[0:i+1]) for i in range(n_tpts)]
        # magnus_terms.append(o3t[-1])
        # ------------------------------------------------------------------------------------------
        # Deprecated (slower, not sure numba is handled properly):
        # t1, t2, o3t = [], [], []
        # for i in prange(n_tpts):
        #     t1.append(-0.5*commutator(o2t[i], At[i]))
        #     t2.append((1.0/12.0)*commutator(o1t[i], commutator(o1t[i], At[i])))
        # o3t_integrand = t1+t2
        # for i in prange(n_tpts):
        #     o3t.append(matrix_integral(times[0:i+1], o3t_integrand[0:i+1]))
        # magnus_terms.append(o3t[-1])
        # ------------------------------------------------------------------------------------------

    # Precompute the Omega_4(t) terms, integrating from t0 to t = t0, ..., t1
    if order >= 4:
        t1 = -0.5 * np.stack([commutator(o3t[i], At[i]) for i in range(n_tpts)], axis=0)
        t2 = f1 * np.stack(
            [commutator(o1t[i], commutator(o2t[i], At[i])) \
                + commutator(o2t[i], commutator(o1t[i], At[i])) for i in range(n_tpts)], axis=0
        )
        o4t_integrand = t1 + t2
        if (integration_method == 'trapezoid'):
            o4t = sp.integrate.cumulative_trapezoid(o4t_integrand, x=times, axis=0, initial=0)
        elif (integration_method == 'simpson'):
            # o4t = sp.integrate.cumulative_simpson(o4t_integrand, x=times, axis=0, initial=0)
            o4t = integral_cumulative_simpson(times, o4t_integrand) 
        magnus_terms.append(o4t[-1])

    # Precompute the Omega_5(t) terms, integrating from t0 to t = t0, ..., t1
    if order >= 5:
        t1 = -0.5 * np.stack([commutator(o4t[i], At[i]) for i in range(n_tpts)], axis=0)
        t2 = f1 * np.stack(
            [commutator(o1t[i], commutator(o3t[i], At[i])) \
            + commutator(o2t[i], commutator(o2t[i], At[i])) \
            + commutator(o3t[i], commutator(o1t[i], At[i])) for i in range(n_tpts)], axis=0
        )
        t3 = f2 * np.stack(
            [commutator(o1t[i], commutator(o1t[i], commutator(o1t[i], 
            commutator(o1t[i], At[i])))) for i in range(n_tpts)], axis=0
        )
        o5t_integrand = t1 + t2 + t3
        if (integration_method == 'trapezoid'):
            o5t = sp.integrate.cumulative_trapezoid(o5t_integrand, x=times, axis=0, initial=0)
        elif (integration_method == 'simpson'):
            # o5t = sp.integrate.cumulative_simpson(o5t_integrand, x=times, axis=0, initial=0)
            o5t = integral_cumulative_simpson(times, o5t_integrand)
        magnus_terms.append(o5t[-1])

    # Precompute the Omega_6(t) terms, integrating from t0 to t = t0, ..., t1
    if order >= 6:
        t1 = -0.5 * np.stack([commutator(o5t[i], At[i]) for i in range(n_tpts)], axis=0)
        t2 = f1 * np.stack(
            [commutator(o1t[i], commutator(o4t[i], At[i]))  \
            + commutator(o2t[i], commutator(o3t[i], At[i])) \
            + commutator(o3t[i], commutator(o2t[i], At[i])) \
            + commutator(o4t[i], commutator(o1t[i], At[i])) for i in range(n_tpts)], axis=0
        )
        t3 = f2 * np.stack(
            [commutator(o1t[i], commutator(o1t[i], commutator(o1t[i], commutator(o2t[i], At[i])))) \
            + commutator(o1t[i], commutator(o1t[i], commutator(o2t[i], commutator(o1t[i], At[i]))))\
            + commutator(o1t[i], commutator(o2t[i], commutator(o1t[i], commutator(o1t[i], At[i]))))\
            + commutator(o2t[i], commutator(o1t[i], commutator(o1t[i], commutator(o1t[i], At[i]))))\
            for i in range(n_tpts)], axis=0)
        o6t_integrand = t1 + t2 + t3
        if (integration_method == 'trapezoid'):
            o6t = sp.integrate.cumulative_trapezoid(o6t_integrand, x=times, axis=0, initial=0)
        elif (integration_method == 'simpson'):
            # o6t = sp.integrate.cumulative_simpson(o6t_integrand, x=times, axis=0, initial=0)
            o6t = matrix_integral_simpson(times, o6t_integrand)
        magnus_terms.append(o6t[-1])

    return magnus_terms

# Function to compute the matrix exponential using Magnus expansion
def magnus_expansion(A, t0, t1, n_tpts=50, order=2, integration_method='trapezoid'):
    """
    Compute the matrix exponential of A(t) from t0 to t1 using the Magnus expansion.
    """
    magnus_terms = compute_magnus_terms(A, t0, t1, n_tpts=n_tpts, order=order,
        integration_method=integration_method)
    # print(magnus_terms)
    Omega = sum(magnus_terms)  # Sum the Magnus terms
    # return sp.linalg.cosm(Omega) + 1j*sp.linalg.sinm(Omega)
    return sp.linalg.expm(Omega)


if __name__ == "__main__":
    def A(t):
        return np.array([[1.0*t, 2.0+3j*t],[2.0-3j*t, 2.0]])
    t0, t1 = 0.0, 1.0
    exp_Omega = magnus_expansion(A, t0, t1, n_tpts=100, order=1, integration_method='trapezoid')
    print(exp_Omega)
    exp_Omega = magnus_expansion(A, t0, t1, n_tpts=100, order=1, integration_method='simpson')
    print(exp_Omega)


