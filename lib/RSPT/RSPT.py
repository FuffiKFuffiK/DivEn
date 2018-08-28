"""Library containing methods for calculations and
numerical analysis of molecular Hamiltonians
"""

import pandas as pd
import numpy as np
import time
import sys

import supp
import harm_oscill


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

    zero_states.sort_values(by=['E'] + labels[:-1], inplace=True)
    zero_states.index = range(len(zero_states))

    return zero_states


@supp.timing
def fill_wmat(anh_coefs, zero_states):
    """Procedure to fill perturbation matrix as a
    system of anharmonic oscillators

    Parameters
    ----------
    anh_coefs: Pandas DataFrame
        Contains anharmonic coefficients and theit indices
        Note: number of indices equals length of the frequencies
        vector
    zero_states: Pandas DataFrame
        Contains zero order approximation (energies) and
        corresponding quantum numbers

    Returns
    -------
    Perturbation matrix of the system of anharmonic oscillators
    """

    N = len(zero_states)
    Nv = len(zero_states.columns) - 1
    Wmat = np.empty(shape=(N, N))
    for i in range(N):
        for j in range(i, N):
            Wmat[i, j] = 0
            for k, coef in enumerate(anh_coefs['k']):
                time1 = []
                time1.append(time.time())
                v1 = np.array([zero_states.iloc[i]['v' + str(l + 1)]
                               for l in range(Nv)], dtype=np.int32)
                time1.append(time.time())
                v2 = v1 if i == j else np.array([zero_states.iloc[j]['v' + str(l + 1)]
                                                 for l in range(Nv)], dtype=np.int32)
                time1.append(time.time())
                dv = abs(v2 - v1)
                time1.append(time.time())
                
                v = v1 if i == j else [max(v1[l], v2[l]) for l in range(Nv)]
                time1.append(time.time())
                n = np.array([anh_coefs.iloc[k]['ik' + str(l + 1)] for l in range(Nv)])
                time1.append(time.time())
                if not any(harm_oscill.ZeroEl(n[l], dv[l]) for l in range(Nv)):
                    Wmat[i, j] += np.prod([harm_oscill.MatElWeight(n[l], dv[l], v[l])
                                          for l in range(Nv)]) * coef
                time1.append(time.time())
            Wmat[j, i] = Wmat[i, j]
            time1.append(time.time())
        time2 = [time1[t] - time1[t - 1] for t in range(1, len(time1))]
        print(time2)
        sys.exit()
    return Wmat


def fill_wmat2(anh_coefs, zero_states):
    """Procedure to fill perturbation matrix as a
    system of anharmonic oscillators

    Parameters
    ----------
    anh_coefs: Pandas DataFrame
        Contains anharmonic coefficients and theit indices
        Note: number of indices equals length of the frequencies
        vector
    zero_states: Pandas DataFrame
        Contains zero order approximation (energies) and
        corresponding quantum numbers

    Returns
    -------
    Perturbation matrix of the system of anharmonic oscillators
    """

    v = []
    dv = []
    vmax = []
    N = len(zero_states)
    Nv = len(zero_states.columns) - 1
    Hmat = np.zeros(shape=(N, N))

    for i in range(Nv):
        v.append(np.array(zero_states['v' + str(i + 1)]))
        dv.append(abs(v[-1] - v[-1][:, np.newaxis]))
        vmax.append(np.maximum(v[-1], v[-1][:, np.newaxis]))

    for k, coef in enumerate(anh_coefs['k']):
        for i in range(Nv):
            delta = anh_coefs.iloc[k]['ik' + str(i + 1)]
            while delta >= 0:
                indices = dv[i].where(dv[i] == delta)
                Hmat[indices] += harm_oscill.MatElWeight(k, delta, v[i][indices]) * coef


    return 0