"""Library containing methods for calculations and
numerical analysis of molecular Hamiltonians
"""

import pandas as pd

import supp


def gen_zero_approximation(E, Emax, freqs, memo=None, n=None):
    """Recurrent method generating values of vibrational quantum numbers
    and energies for states with E0 < Emax

    Parameters
    ----------
    E: float
        Zero energy level. At first call should be equal to
        sum(frequencies)/2
    Emax: float
        Maximum value of energy for which calculations are performed
    freqs: list of floats
        List of normal vibrational frequencies of the molecule
    memo: list of ints
        List containing current quantum numbers. Length of it
        is equal to the length of freqs
    n: int
        Number of current frequency in recursion
    """
    if memo is None:
        memo = [0] * len(freqs)
    if n is None:
        n = 0
    if n >= len(freqs):
        yield (*memo, E)
    elif E <= Emax:
        yield from gen_zero_approximation(E, Emax, freqs, memo=memo, n=n+1)
        memo[n] += 1
        yield from gen_zero_approximation(E + freqs[n], Emax, freqs, memo=memo, n=n)
    else:
        memo[n] = 0


@supp.timing
def zero_approximation(Emax, freqs):
    """Generates Pandas DataFrame with states and corresponding
    vibrational quantum numbers and zero energies under condition E0 < Emax

    Parameters
    ----------
    Emax: float
        Maximum value of energy for which calculations are performed
    freqs: list of floats
        List of normal vibrational frequencies of the molecule

    Returns
    -------
    Pandas DataFrames with columns, labelled v1, v2, ..., vk, E
    where k is a number of normal vibrational modes in the molecule
    and E are zero order approximation to energies
    """
    labels = []
    for i in range(len(freqs)):
        labels.append('v' + str(i + 1))
    labels += ['E']
    zero_states = pd.DataFrame(data=gen_zero_approximation(sum(freqs) / 2, Emax, freqs),
                               columns=labels)

    return zero_states
