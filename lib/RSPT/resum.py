"""Module containing different resummation techniques for
convergent and divergent series
"""


def PH_approx(pt_series, orders):
    """Calculates coefficients of polynomials, representing Pade-Hermite
    approximants for the function, represented by power series with
    coefficients stored in pt_series list. The degree is equal to
    the length of orders list, orders of each polynomial are given
    by the values of this list.

    Parameters
    ----------
    pt_series: List of floats or complex
        Contains coefficients of the function's Taylor expansion.
    orders: List of int
        Contains list of desired orders of polynomials used for
    the representaions of Pade-Hermite approximant with the
    degree, equal to len(orders) - 1

    Implementation note: for correct solution of SLAU the following
    condition should be satisfied:
    len(pt_series) >= sum(orders) + len(orders) - 1
    """

    pass