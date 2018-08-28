"""Contains methods describing matrix elements of
harmonic oscillator (or system of harmonic oscillators)
"""


def MatElWeight(n, dv, v):
    """Calculates matrix elements weights of harmonic oscillator

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


def MatElWeightGen(vmax):
    """Generates matrix element weights of harmonic oscillator

    Parameters
    ----------
    v: int
        Maximum value of quantum number v (for any frequency)

    Returns:
    Dictionary of non-zero matrix elements 'weight' of harmonic oscillator.
    Keys are tuples (n, dv, v)

    n >= 0, dv >= 0, n >= dv
    Note: last condition is the limitation of the program
    It can be extended by implementing formulas for
    matrix element for harmonic oscillator for the case
    n (or dv) > 8 (anharmonic functions of higher order)
    """

    #n = 0
    for _ in range(vmax + 1):
        yield 1
    #n = 1
    for v in range(vmax + 1):
        yield (v / 2) ** 0.5
    #n = 2
    for v in range(vmax + 1):
        yield v + 0.5
    for v in range(vmax + 1):
        yield 0.5 * (v * (v - 1)) ** 0.5
    #n = 3
    for v in range(vmax + 1):
        yield 0.75 * 2 ** 0.5 * v ** 1.5
    for v in range(vmax + 1):
        yield 0.25 * (2 * v * (v - 1) * (v - 2)) ** 0.5
    #n = 4
    for v in range(vmax + 1):
        yield 1.5 * v ** 2 + 1.5 * v + 0.75
    for v in range(vmax + 1):
        yield (v - 0.5) * (v * (v - 1)) ** 0.5
    for v in range(vmax + 1):
        yield 0.25 * (v * (v - 1) * (v - 2) * (v - 3)) ** 0.5
    #n = 5
    for v in range(vmax + 1):
        yield (1.25 * v ** 2 + 0.625) * (2 * v) ** 0.5
    for v in range(vmax + 1):
        yield (0.625 * v - 0.625) * (2 * v * (v - 1) * (v - 2)) ** 0.5
    for v in range(vmax + 1):
        yield 0.125 * (2 * v * (v - 1) * (v - 2) * (v - 3) * (v - 4)) ** 0.5
    #n = 6
    for v in range(vmax + 1):
        yield 0.625 * (2 * v + 1) * (2 * v ** 2 + 2 * v + 3)
    for v in range(vmax + 1):
        yield 1.875 * (v ** 2 - v + 1) * (v * (v - 1)) ** 0.5
    for v in range(vmax + 1):
        yield (0.75 * v - 1.125) * (v * (v - 1) * (v - 2) * (v - 3)) ** 0.5
    for v in range(vmax + 1):
        yield 0.125 * (v * (v - 1) * (v - 2) * (v - 3) * (v - 4) * (v - 5)) ** 0.5
    #n = 7
    for v in range(vmax + 1):
        yield 2.1875 * (v ** 2 + 2) * v * (2 * v) ** 0.5
    for v in range(vmax + 1):
        yield 1.3125 * (v ** 2 - 2 * v + 2) * (2 * v * (v - 1) * (v - 2)) ** 0.5
    for v in range(vmax + 1):
        yield 0.4375 * (v - 2) * (2 * v * (v - 1) * (v - 2) * (v - 3) * (v - 4)) ** 0.5
    for v in range(vmax + 1):
        yield 0.0625 * (2 * v * (v - 1) * (v - 2) * (v - 3) * (v - 4) * (v - 5) * (v - 6)) ** 0.5
    #n = 8
    for v in range(vmax + 1):
        yield 4.375 * (v ** 4 + 2 * v ** 3 + 5 * v ** 2 + 4 * v + 1.5)
    for v in range(vmax + 1):
        yield 1.75 * (2 * v - 1) * (v ** 2 - v + 3) * (v * (v - 1)) ** 0.5
    for v in range(vmax + 1):
        yield 0.875 * (2 * v ** 2 - 6 * v + 7) * (v * (v - 1) * (v - 2) * (v - 3)) ** 0.5
    for v in range(vmax + 1):
        yield (0.5 * v - 1.25) * (v * (v - 1) * (v - 2) * (v - 3) * (v - 4) * (v - 5)) ** 0.5
    for v in range(vmax + 1):
        yield 0.0625 * (v * (v - 1) * (v - 2) * (v - 3) * (v - 4)
                        * (v - 5) * (v - 6) * (v - 7)) ** 0.5


def MatElWeightDict(vmax, n=8):
    """Generates dictionary of non-zero matrix elements 'weight' of harmonic oscillator.
    Keys are tuples (n, dv, v)

    Parameters
    ----------
    n: int
        Maximal index of corresponding anharmonic coefficient
    v: int
        Maximum value of quantum number v (for any frequency)

    Note: n <= 8 is the limitation of the program
    It can be extended by implementing formulas for
    matrix element for harmonic oscillator for the case
    n (or dv) > 8 (anharmonic functions of higher order)
    """

    Harm_weight_dict = {}
    gen = MatElWeightGen(vmax)
    for i in range(n + 1):
        for dv in range(i % 2, i + 1, 2):
            for v in range(vmax + 1):
                Harm_weight_dict[(i, dv, v)] = next(gen)

    return Harm_weight_dict
