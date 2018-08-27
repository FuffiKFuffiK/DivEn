"""Contains methods describing matrix elements of
harmonic oscillator (or system of harmonic oscillators)
"""


def ZeroEl(n, dv):
    """Quickly checks if matrix element equals zero

    Parameters
    ----------
    n: int
        Index of corresponding anharmonic coefficient
    dv: int
        Absolute difference of corresponding quantum numbers

    Returns:
    True if matrix element is zero, False otherwise

    Note: last condition is the limitation of the program
    It can be extended by implementing formulas for
    matrix element for harmonic oscillator for the case
    n (or dv) > 8 (anhatmonic funvtions of higher order)
    """

    return dv > n or (n - dv) % 2 or n > 8


def MatElWeight(n, dv, v):
    """Calculates matrix element of harmonic oscillator

    Parameters
    ----------
    n: int
        Index of corresponding anharmonic coefficient
    dv: int
        Absolute difference of corresponding quantum numbers
        dv = abs(v1 - v2)
    v: int
        Max value of quantum number pair
        v = max(v1, v2)

    Returns:
    Matrix element 'weight' of harmonic oscillator
    Note: for correct implementation, first call ZeroEl

    Note: last condition is the limitation of the program
    It can be extended by implementing formulas for
    matrix element for harmonic oscillator for the case
    n (or dv) > 8 (anharmonic funvtions of higher order)
    """

    if n == 0:
        return 1
    elif n == 1:
        return (v / 2) ** 0.5
    elif n == 2:
        return v + 0.5 if dv == 0 else 0.5 * (v * (v - 1)) ** 0.5
    elif n == 3:
        if dv == 1:
            return 0.75 * 2 ** 0.5 * v ** 1.5
        else:
            return 0.25 * (2 * v * (v - 1) * (v - 2)) ** 0.5
    elif n == 4:
        if dv == 0:
            return 1.5 * v ** 2 + 1.5 * v + 0.75
        elif dv == 2:
            return (v - 0.5) * (v * (v - 1)) ** 0.5
        else:
            return 0.25 * (v * (v - 1) * (v - 2) * (v - 3)) ** 0.5
    elif n == 5:
        if dv == 1:
            return (1.25 * v ** 2 + 0.625) * (2 * v) ** 0.5
        elif dv == 3:
            return (0.625 * v - 0.625) * (2 * v * (v - 1) * (v - 2)) ** 0.5
        else:
            return 0.125 * (2 * v * (v - 1) * (v - 2) * (v - 3) * (v - 4)) ** 0.5
    elif n == 6:
        if dv == 0:
            return 0.625 * (2 * v + 1) * (2 * v ** 2 + 2 * v + 3)
        elif dv == 2:
            return 1.875 * (v ** 2 - v + 1) * (v * (v - 1)) ** 0.5
        elif dv == 4:
            return (0.75 * v - 1.125) * (v * (v - 1) * (v - 2) * (v - 3)) ** 0.5
        else:
            return 0.125 * (v * (v - 1) * (v - 2) * (v - 3) * (v - 4) * (v - 5)) ** 0.5
    elif n == 7:
        if dv == 1:
            return 2.1875 * (v ** 2 + 2) * v * (2 * v) ** 0.5
        elif dv == 3:
            return 1.3125 * (v ** 2 - 2 * v + 2) * (2 * v * (v - 1) * (v - 2)) ** 0.5
        elif dv == 5:
            return 0.4375 * (v - 2) * (2 * v * (v - 1) * (v - 2) * (v - 3) * (v - 4)) ** 0.5
        else:
            return 0.0625 * (2 * v * (v - 1) * (v - 2) * (v - 3) *
                             (v - 4) * (v - 5) * (v - 6)) ** 0.5
    elif n == 8:
        if dv == 0:
            return 4.375 * (v ** 4 + 2 * v ** 3 + 5 * v ** 2 + 4 * v + 1.5)
        elif dv == 2:
            return 1.75 * (2 * v - 1) * (v ** 2 - v + 3) * (v * (v - 1)) ** 0.5
        elif dv == 4:
            return 0.875 * (2 * v ** 2 - 6 * v + 7) * (v * (v - 1) * (v - 2) * (v - 3)) ** 0.5
        elif dv == 6:
            return (0.5 * v - 1.25) * (v * (v - 1) * (v - 2) * (v - 3) * (v - 4) * (v - 5)) ** 0.5
        else:
            return 0.0625 * (v * (v - 1) * (v - 2) * (v - 3) *
                             (v - 4) * (v - 5) * (v - 6) * (v - 7)) ** 0.5
    else:
        return 0
