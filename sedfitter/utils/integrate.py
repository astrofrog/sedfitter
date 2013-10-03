from __future__ import print_function, division

import numpy as np

from .interpolate import interp1d_fast


def integrate_subset(x, y, xmin, xmax):
    """
    Perform trapezium integration of a set of points (x,y) between bounds xmin
    and xmax. The interpolation between the points is done in linear space, so
    this is designed for functions that are piecewise linear in linear space.
    """

    # Swap arrays if necessary
    if x[-1] < x[0]:
        x = x[::-1]
        y = y[::-1]

    # Swap limits if necessary
    if xmin > xmax:
        xmin, xmax = xmax, xmin
    elif xmin == xmax:
        return 0.

    # Find the subset of points to use and the value of the function at the
    # end-points of the integration

    if xmin == x[0]:
        i1 = 1
        ymin = y[0]
    else:
        i1 = np.searchsorted(x, xmin)
        ymin = interp1d_fast(x[i1 - 1:i1 + 1], y[i1 - 1:i1 + 1], xmin)

    if xmax == x[-1]:
        i2 = -2
        ymax = y[-1]
    else:
        i2 = np.searchsorted(x, xmax)
        ymax = interp1d_fast(x[i2 - 1:i2 + 1], y[i2 - 1:i2 + 1], xmax)

    # Construct sub-arrays of the relevant data
    x = np.hstack([xmin, x[i1:i2], xmax])
    y = np.hstack([ymin, y[i1:i2], ymax])

    # Call function to integrate the whole subset
    return integrate(x, y)


def integrate(x, y):
    """
    Perform trapezium integration of a set of points (x,y). The interpolation
    between the points is done in linear space, so this is designed for
    functions that are piecewise linear in linear space.
    """

    # Fix NaN values
    y[np.isnan(y)] = 0.

    # Find the integral of all the chunks
    integrals = 0.5 * (x[1:] - x[:-1]) * (y[1:] + y[:-1])

    # Sum them all up
    integral = np.sum(integrals)

    # Check if the integral is NaN or infinity
    if np.isnan(integral) or np.isinf(integral):
        raise Exception("Integral is NaN or Inf")

    return integral
