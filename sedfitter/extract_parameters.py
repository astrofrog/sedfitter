from __future__ import print_function, division

import os
try:
    import cPickle as pickle
except ImportError:
    import pickle

from astropy.table import Table
import numpy as np

from .fit_info import FitInfo
from .extinction import Extinction


def extract_parameters(input=None, output_prefix=None, output_suffix=None,
                       parameters='all', select_format=("N", 1), header=True):

    fin = open(input, 'rb')

    # Read in header of output file
    model_dir = pickle.load(fin)
    filters = pickle.load(fin)
    extinction_law = pickle.load(fin)

    # Read in table of parameters for model grid
    if os.path.exists(model_dir + '/parameters.fits'):
        t = Table.read(model_dir + '/parameters.fits')
    elif os.path.exists(model_dir + '/parameters.fits.gz'):
        t = Table.read(model_dir + '/parameters.fits.gz')
    else:
        raise Exception("Parameter file not found in %s" % model_dir)

    format = {}
    for par in t.dtype.names:
        if par == 'MODEL_NAME':
            format[par] = "%11s"
        else:
            format[par] = "%11.3e"

    if parameters == 'all':
        parameters = t.dtype.names

    while True:

        # Read in next fit
        try:
            info = pickle.load(fin)
        except EOFError:
            break

        # Filter fits
        info.keep(select_format[0], select_format[1])

        output = ""
        if output_prefix:
            output += output_prefix
        output += info.source.name
        if output_suffix:
            output += output_suffix

        fout = open(output, 'w')

        if header:
            fout.write("%11s %11s %11s " % ("CHI2", "AV", "SC"))
            fout.write(" ".join([("%11s" % par) for par in parameters]))
            fout.write("\n")

        # Match good-fitting models to parameter list
        subset = np.in1d(t['MODEL_NAME'], info.model_name)

        tsub = t[subset]
        index = np.argsort(np.argsort(info.model_name))
        tsorted = tsub[index]
        if not np.all(info.model_name == tsorted['MODEL_NAME']):
            raise Exception("Parameter file sorting failed")

        for i in range(info.n_fits):

            row = tsorted[i]

            basic = "%11.3e %11.3e %11.3e " % (info.chi2[i], info.av[i], info.sc[i])
            pars = " ".join([(format[par] % row[par]) for par in parameters])

            fout.write(basic + pars + "\n")

        fout.close()

    # Close input and output files
    fin.close()
