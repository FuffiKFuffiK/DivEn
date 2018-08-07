"""Contains supplementary function to handle errors and exceptions and
    methods to compute computational time
"""

import linecache
import sys
import os
import types
from functools import wraps
from time import time


def PrintException():
    """Prints filename, line number, line and type of the exception in a nice format
    """
    exc_type, exc_obj, tb = sys.exc_info()
    f = tb.tb_frame
    lineno = tb.tb_lineno
    filename = os.path.basename(f.f_code.co_filename)
    linecache.checkcache(filename)
    line = linecache.getline(filename, lineno, f.f_globals)
    print('EXCEPTION IN ({}, LINE {} "{}"): {}'.format(filename, lineno, line.strip(), exc_obj))

    return(filename, lineno, line.strip(), exc_obj, exc_type)


def PrintSysExitInfo(back=0):

    frame = sys._getframe(back+1)
    filename = os.path.basename(frame.f_code.co_filename)
    lineno = frame.f_lineno - 1
    line = linecache.getline(filename, lineno)
    while not line.strip() and lineno:
        lineno -= 1
        line = linecache.getline(filename, lineno)
    print('SYSTEM EXIT IN ({}, LINE {} "{}")'.format(filename, lineno, line.strip()))

    return(filename, lineno, line.strip())


def timing(f):
    """Decorator to calculate the computational time of the call of the
    method f

    Parameters
    ----------
    f: function
    """
    @wraps(f)
    def wrap(*args, **kw):
        ts = time()
        result = f(*args, **kw)
        te = time()
        print('function: \'{}\' took: {:2.4f} sec'.format(f.__name__, te-ts))
        return result
    return wrap


class Reader():
    """Class to read file objects (and not only) with the help of generators
    """
    def __init__(self, g):
        self.g = g

    @property
    def g(self):
        return self.__g

    @g.setter
    def g(self, g):
        if isinstance(g, types.GeneratorType):
            self.__g = g
        else:
            PrintSysExitInfo()
            raise TypeError('An input to Reader class should be a generator function')

    def read(self, n=0):
        try:
            return next(self.g)
        except StopIteration:
            return ''
