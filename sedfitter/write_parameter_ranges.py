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


def write_parameter_ranges(input_file, output_file, select_format=("N", 1), additional={}):
    """
    Write out an ASCII file with ranges of paramters for each source.

    Parameters
    ----------
    input_file : str
        File containing the fit information
    output_file : str, optional
        The output ASCII file containing the parameter ranges
    select_format : tuple, optional
        Tuple specifying which fits should be output. See the documentation
        for a description of the tuple syntax.
    additional : dict, optional
        A dictionary giving additional parameters for each model. This should
        be a dictionary where each key is a parameter, and each value is a
        dictionary mapping the model names to the parameter values.
    """

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
    fout.write("%30s " % "")

    fout.write("%10s " % "")
    fout.write("%10s " % "")

    fout.write('chi2'.center(32) + ' ')
    fout.write('av'.center(32) + ' ')
    fout.write('scale'.center(32) + ' ')

    for par in list(t.columns.keys()) + list(additional.keys()):
        if par == 'MODEL_NAME':
            continue
        fout.write(par.lower().center(32) + ' ')

    fout.write('\n')

    # Second header line
    fout.write("source_name".center(30) + ' ')
    fout.write("n_data".center(10) + ' ')
    fout.write("n_fits".center(10) + ' ')

    fout.write("min".center(10) + " " + "best".center(10) + " " + "max".center(10) + " ")
    fout.write("min".center(10) + " " + "best".center(10) + " " + "max".center(10) + " ")
    fout.write("min".center(10) + " " + "best".center(10) + " " + "max".center(10) + " ")

    for par in list(t.columns.keys()) + list(additional.keys()):
        if par == 'MODEL_NAME':
            continue
        fout.write("min".center(10) + " " + "best".center(10) + " " + "max".center(10) + " ")

    fout.write('\n')

    # Third header line
    fout.write('-' * 30 + ' ')
    fout.write('-' * 10 + ' ')
    fout.write('-' * 10 + ' ')

    fout.write('-' * 32 + ' ')
    fout.write('-' * 32 + ' ')
    fout.write('-' * 32 + ' ')

    for par in list(t.columns.keys()) + list(additional.keys()):
        if par == 'MODEL_NAME':
            continue
        fout.write('-' * 32 + ' ')

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

        fout.write('%10.3f %10.3f %10.3f ' % (info.chi2.min(), info.chi2[0], info.chi2.max()))
        fout.write('%10.3f %10.3f %10.3f ' % (info.av.min(), info.av[0], info.av.max()))
        fout.write('%10.3f %10.3f %10.3f ' % (info.sc.min(), info.sc[0], info.sc.max()))

        for par in tsorted.columns:
            if par == 'MODEL_NAME':
                continue
            fout.write('%10.3e %10.3e %10.3e ' % (tsorted[par].min(), tsorted[par][0], tsorted[par].max()))

        fout.write('\n')

    # Close input and output files
    fin.close()
    fout.close()
    fout.close()
