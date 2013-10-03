from __future__ import print_function, division

import os
try:
    import cPickle as pickle
except ImportError:
    import pickle

import numpy as np
from astropy.table import Table

from .fit_info import FitInfo
from .extinction import Extinction
from .models import load_parameter_table


def write_parameters(input_file, output_file, select_format=("N", 1), additional={}):

    # Open input and output file
    fin = open(input_file, 'rb')
    fout = open(output_file, 'w')

    # Read in header of output file
    model_dir = pickle.load(fin)
    filters = pickle.load(fin)
    extinction = pickle.load(fin)

    # Read in table of parameters for model grid
    t = load_parameter_table(model_dir)

    t['MODEL_NAME'] = np.char.strip(t['MODEL_NAME'])
    t.sort('MODEL_NAME')

    # First header line
    fout.write("source_name".center(30) + ' ')
    fout.write("n_data".center(10) + ' ')
    fout.write("n_fits".center(10) + ' ')

    fout.write('\n')

    # Second header line

    fout.write('fit_id'.center(10) + ' ')
    fout.write('model_name'.center(30) + ' ')
    fout.write('chi2'.center(10) + ' ')
    fout.write('av'.center(10) + ' ')
    fout.write('scale'.center(10) + ' ')

    for par in list(t.columns.keys()) + list(additional.keys()):
        if par == 'MODEL_NAME':
            continue
        fout.write(par.lower().center(10) + ' ')

    fout.write('\n')

    fout.write('-' * (75 + 11 * (len(list(t.columns.keys()) + list(additional.keys())))))
    fout.write('\n')

    while True:

        # Read in next fit
        try:
            info = pickle.load(fin)
        except EOFError:
            break

        # Filter fits
        info.keep(select_format[0], select_format[1])

        # Get filtered and sorted table of parameters
        tsorted = info.filter_table(t)

        # Add additional parameter columns if necessary
        for parameter in additional:
            if parameter in tsorted.columns:
                raise Exception("Parameter {} already exists in file".format(parameter))
            tsorted.add_empty_column(parameter, dtype=float)
            for i, name in enumerate(tsorted['MODEL_NAME']):
                tsorted[parameter][i] = additional[parameter][name.strip()]

        fout.write("%30s " % info.source.name)
        fout.write("%10i " % info.source.n_data)
        fout.write("%10i " % info.n_fits)
        fout.write("\n")

        for fit_id in range(len(info.chi2)):
            fout.write('%10i ' % (fit_id + 1))
            fout.write('%30s ' % info.model_name[fit_id])
            fout.write('%10.3f ' % info.chi2[fit_id])
            fout.write('%10.3f ' % info.av[fit_id])
            fout.write('%10.3f ' % info.sc[fit_id])

            for par in tsorted.columns:
                if par == 'MODEL_NAME':
                    continue
                fout.write('%10.3e ' % (tsorted[par][fit_id]))

            fout.write('\n')

    # Close input and output files
    fin.close()
    fout.close()
    fout.close()
