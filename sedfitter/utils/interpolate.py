from __future__ import print_function, division

import numpy as np


class check_bounds(object):
    """
    Decorator to add standard interpolation bounds checking.

    This decorator incurs an overhead of 30-50 percent on the runtime of
    the interpolation routine.
    """

    def __init__(self, f):
        self.f = f

    def __call__(self, x, y, xval, bounds_error=True, fill_value=np.nan):
        if np.isscalar(xval):  # xval is a scalar
            if xval < x[0] or xval > x[-1]:  # the value is out of bounds
                if bounds_error:
                    raise Exception("x value is out of interpolation bounds")
                else:
                    return np.nan
            else:  # the value is in the bounds
                return self.f(x, y, xval)
        else:  # xval is an array
            inside = (xval >= x[0]) & (xval <= x[-1])
            outside = ~inside
            if np.any(outside):  # some values are out of bounds
                if bounds_error:
                    raise Exception("x values are out of interpolation bounds")
                else:
                    if np.any(inside):
                        yval = np.zeros(xval.shape)
                        yval[inside] = self.f(x, y, xval[inside])
                        yval[outside] = fill_value
                        return yval
                    else:
                        return np.repeat(fill_value, xval.shape)
            else:  # all values are in the bounds
                return self.f(x, y, xval)


@check_bounds
def interp1d_fast(x, y, xval):
    """On-the-fly linear interpolator"""
    if len(x) != len(y):
        raise Exception("x and y should have the same length")
    ipos = np.searchsorted(x, xval)
    return (xval - x[ipos - 1]) \
        / (x[ipos] - x[ipos - 1]) \
        * (y[ipos] - y[ipos - 1]) \
        + y[ipos - 1]
