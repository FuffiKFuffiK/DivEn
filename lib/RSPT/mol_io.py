"""Module containing input/output methods for DivEn project
"""

import pandas as pd


def read_freqs(fname):
    """Method to read frequencies of normal vibrations of the molecule into Pandas
    DataFrame
    """

    return pd.read_csv(fname, delim_whitespace=True, header=None, names=['omega'])


def read_anh_coef(fname):
    """Method to read anharmonic coefficents of the potential energy function
    of the molecule into Pandas DataFrame
    """

    f = open(fname)

    labels = []
    for i in range(len(f.readline.split() - 1)):
        labels.append('k' + str(i))
    labels.append('k')

    return pd.read_csv(fname, delim_whitespace=True, header=None, names=labels)
