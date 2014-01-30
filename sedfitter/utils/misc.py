import numpy as np


def order_to_match(array, reference):
    """
    Given an array ``array``, return the index array needed to make it the same as ``reference``
    """
    return np.argsort(array)[np.argsort(np.argsort(reference))]
