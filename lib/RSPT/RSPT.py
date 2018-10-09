"""Library containing methods for calculations and
numerical analysis of molecular Hamiltonians
"""

import pandas as pd
import numpy as np

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

    v = []
    dv = []
    vmat = []
    N = len(zero_states)
    Nv = len(zero_states.columns) - 1
    Wmat = np.zeros(shape=(N, N))

    nmax = anh_coefs[anh_coefs.columns[:-1]].max().max()
    vmax = zero_states[zero_states.columns[:-1]].max().max()

    Weights_dict = harm_oscill.MatElWeightDict(vmax, n=nmax)

    coef_ind = np.array(anh_coefs[anh_coefs.columns[:-1]])
    for i in range(Nv):
        v.append(np.array(zero_states['v' + str(i + 1)]))
        dv.append(abs(v[-1] - v[-1][:, np.newaxis]))
        vmat.append(np.maximum(v[-1], v[-1][:, np.newaxis]))

    for i in range(N):
        for j in range(i, N):
            for k, coef in enumerate(anh_coefs['k']):
                p = [Weights_dict.get((coef_ind[k][l], dv[l][i][j], vmat[l][i][j]))
                     for l in range(Nv)]
                if all(p):
                    Wmat[i, j] += np.prod(p) * coef

    i_lower = np.tril_indices(N, -1)
    Wmat[i_lower] = Wmat.T[i_lower]

    return Wmat


def add_zero_approx(Wmat, zero_states):
    """Procedure to add zero order approximation
    to the perturbation matrix

    Parameters
    ----------
    Wmat: Numpy 2d array
        Contains perturbation matrix
    zero_states: Pandas DataFrame
        Contains zero order approximation (energies) and
        corresponding quantum numbers

    Returns
    -------
    Hamiltonian matrix H = W + E0
    """

    Hmat = Wmat.copy()
    Hmat.flat[::Hmat.shape[0] + 1] += zero_states['E']
    return Hmat


@supp.timing
def diag_hmat(Hmat, zero_states, E0=None):
    """Procedure to add zero order approximation
    to the perturbation matrix

    Parameters
    ----------
    Hmat: Numpy 2d array
        Contains perturbation matrix
    zero_states: Pandas DataFrame
        Contains zero order approximation (energies) and
        corresponding quantum numbers

    Returns
    -------
    Pandas DataFrame, containing vibrational energies,
    retrieved by diagonalisation procedure
    """

    Eigh_values, Eigh_vectors = np.linalg.eigh(Hmat)

    #Creating and filling new DataFrame for vibrational states
    vib_states = zero_states
    #Assignment by largest coefficient in wavefunction
    new_index = np.square(np.absolute(Eigh_vectors)).argmax(axis=1)
    vib_states.index = new_index
    vib_states.sort_index(inplace=True)
    vib_states['E'] = Eigh_values
    vib_states.index = range(len(vib_states))

    if E0 is None:
        E0 = vib_states.iloc[0]['E']

    vib_states['E'] -= E0

    return vib_states, Eigh_vectors
