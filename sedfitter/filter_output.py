from __future__ import print_function, division

try:
    import cPickle as pickle
except ImportError:
    import pickle

from .fit_info import FitInfo
from .extinction import Extinction


def filter_output(input=None, output_good='auto', output_bad='auto', chi=None,
                  cpd=None):

    if input[-8:] != '.fitinfo':
        raise Exception("Extension of input file should be .fitinfo")

    fin = open(input, 'rb')

    if output_good == 'auto':
        fout_good = open(input.replace('.fitinfo', '_good.fitinfo'), 'wb')
    else:
        fout_good = open(output_good, 'wb')

    if output_bad == 'auto':
        fout_bad = open(input.replace('.fitinfo', '_bad.fitinfo'), 'wb')
    else:
        fout_bad = open(output_bad, 'wb')

    # Read in header of output file
    model_dir = pickle.load(fin)
    filters = pickle.load(fin)
    extinction_law = pickle.load(fin)

    # Output header to good fits file
    pickle.dump(model_dir, fout_good, 2)
    pickle.dump(filters, fout_good, 2)
    pickle.dump(extinction, fout_good, 2)

    # Output header to bad fits file
    pickle.dump(model_dir, fout_bad, 2)
    pickle.dump(filters, fout_bad, 2)
    pickle.dump(extinction_law, fout_bad, 2)

    while True:

        # Read in next fit
        try:
            info = pickle.load(fin)
        except EOFError:
            break

        bestchi = info.chi2[0]
        bestcpd = info.chi2[0] / float(info.source.n_data)

        if (chi and bestchi < chi) or (cpd and bestcpd < cpd):
            pickle.dump(info, fout_good, 2)
        else:
            pickle.dump(info, fout_bad, 2)

    # Close input and output files
    fin.close()
    fout_good.close()
    fout_bad.close()
