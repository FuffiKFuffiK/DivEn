"""Module containing input/output methods for DivEn project
"""

import pandas as pd
import numpy as np


def read_freqs(fname):
    """Method to read frequencies of normal vibrations of the molecule into Pandas
    DataFrame
    """

    return pd.read_csv(fname, delim_whitespace=True, header=None, names=['omega'])


def read_anh_coefs(fname):
    """Method to read anharmonic coefficents of the potential energy function
    of the molecule into Pandas DataFrame
    """

    f = open(fname)

    labels = []
    for i in range(len(f.readline().split()) - 1):
        labels.append('ik' + str(i + 1))
    labels.append('k')

    f.close()

    return pd.read_csv(fname, delim_whitespace=True, header=None, names=labels)


def write_states(fname, states):
    """Method to quickly read richmol format states file into a single pandas dataframe

    Parameters
    ----------
    fname : str
        Name of States file.
    states : Pandas DataFrame, containing information about states
    Possible options: zero_states, vib_states, etc.
    """

    fmt = '%4d ' * (len(states.columns) - 1)
    fmt = fmt + '%24.16f'
    np.savetxt(fname, states.values, fmt=fmt)
    return 1
