from __future__ import print_function, division

try:
    import cPickle as pickle
except ImportError:
    import pickle

from atpy import Table, TableSet

from ..fit_info import FitInfo
from ..extinction import Extinction


def fitinfo2tables(file_in):

    # Open input file
    fin = open(file_in, 'rb')

    # Read in header of input file
    model_dir = pickle.load(fin)
    filters = pickle.load(fin)
    extinction = pickle.load(fin)

    # Intialize TableSet
    ts = TableSet()

    while True:

        # Read in next fit
        try:
            info = pickle.load(fin)
        except EOFError:
            break

        # Create Table
        t = Table(name=info.source.name)
        t.add_column('model_name', info.model_name)
        t.add_column('chi2', info.chi2)
        t.add_column('av', info.sc)
        t.add_column('sc', info.sc)

        # Append to TableSet
        ts.append(t)

    # Close input file
    fin.close()

    return ts
